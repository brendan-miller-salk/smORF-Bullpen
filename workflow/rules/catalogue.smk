# Catalogue: merge callers (build_catalogue.py; `catalogue_merge: exact` keeps every
# distinct ORF structure, `overlap` folds by reciprocal overlap), then, if
# `orf_collapse` is true, fold smORFs whose proteins MMseqs2 clusters together
# (>= 90% identity, >= 80% coverage). `orf_collapse: false` skips
# the collapse; the merged catalogue is then the final one. Provided mode (orf_gtf) replaces the caller merge and never
# collapses.

CAT = f"{OUT}/catalogue/cohort.catalogue"
MERGED = f"{OUT}/catalogue/merged/cohort.merged"


def _price_input(_):
    return [f"{OUT}/callers/price/allsamples.orfs.tsv"] if "price" in CALLERS else []


def _ribocode_input(_):
    return (expand(f"{OUT}/callers/ribocode/{{sample}}.txt", sample=SAMPLES)
            if "ribocode" in CALLERS else [])


def _ribotish_input(_):
    return (expand(f"{OUT}/callers/ribotish/{{tis}}.txt", tis=TIS_SAMPLES)
            if "ribotish" in CALLERS else [])


def _opt(flag, key):
    v = config.get(key)
    return f"{flag} {v}" if v is not None else ""


if PROVIDED:

    rule provided_catalogue:
        """Provided mode: the orf_gtf's ORFs, one per transcript, as the merged catalogue.
        Never MMseqs2-collapsed: every structure in the GTF stays a row."""
        input:
            gtf=config.get("orf_gtf", ""),
            genome=config["genome"],
            fai=f"{IDX}/genome.fa.fai",
        output:
            bed=f"{MERGED}.bed12",
            tsv=f"{MERGED}.tsv",
            consensus=f"{MERGED}.consensus.bed12",
            report=f"{OUT}/catalogue/catalogue_report.tsv",
        params:
            max_aa=config.get("smorf_max_aa", 150),
        conda:
            "../envs/score.yaml"
        log:
            f"{OUT}/logs/catalogue.log",
        shell:
            r"""
            python {SCRIPTS}/provided_catalogue.py --gtf {input.gtf} --genome {input.genome} \
                --smorf-max-aa {params.max_aa} --out-prefix {MERGED} --report {output.report} 2> {log}
            cat {log} >&2
            """

    rule provided_orf_ids:
        """The ORF list join_annotations joins onto, in place of the scores table."""
        input:
            f"{SC}/orf_candidates.bed12",
        output:
            f"{OUT}/catalogue/orf_ids.tsv",
        shell:
            "(echo orf_id; cut -f4 {input}) > {output}"

else:

    rule catalogue:
        """Harmonise seqnames, project RiboCode onto exons, merge callers and samples.
        Ribo-TISH calls are one 'sample' per TIS library."""
        input:
            price=_price_input,
            ribocode=_ribocode_input,
            ribotish=_ribotish_input,
            gtf=FULL_HYBRID,
            genome=config["genome"],
            fai=f"{IDX}/genome.fa.fai",
        output:
            bed=f"{MERGED}.bed12",
            tsv=f"{MERGED}.tsv",
            consensus=f"{MERGED}.consensus.bed12",
            report=f"{OUT}/catalogue/catalogue_report.tsv",
        params:
            prefix=MERGED,
            price=lambda wc, input: f"--price {input.price}" if input.price else "",
            ribocode=lambda wc, input: f"--ribocode {input.ribocode}" if input.ribocode else "",
            ribotish=lambda wc, input: f"--ribotish {input.ribotish}" if input.ribotish else "",
            opts=" ".join([
                _opt("--smorf-max-aa", "smorf_max_aa"),
                _opt("--min-callers", "orf_min_callers"),
                _opt("--min-samples", "orf_min_samples"),
                _opt("--price-max-p", "price_max_p"),
                _opt("--price-max-fdr", "price_max_fdr"),
                _opt("--ribocode-max-p", "ribocode_max_p"),
                _opt("--ribotish-max-tis-q", "ribotish_max_tis_q"),
                _opt("--merge", "catalogue_merge"),
            ]),
        resources:
            mem_mb=16000,
            runtime=60,
        conda:
            "../envs/score.yaml"
        log:
            f"{OUT}/logs/catalogue.log",
        shell:
            r"""
            python {SCRIPTS}/build_catalogue.py {params.price} {params.ribocode} {params.ribotish} \
                --gtf {input.gtf} --genome {input.genome} {params.opts} \
                --out-prefix {params.prefix} --report {output.report} 2> {log}
            cat {log} >&2
            """


rule orf_nt_fasta:
    """Every catalogue ORF, spliced and strand-aware; headers orf_id(+/-)."""
    input:
        bed=f"{MERGED}.bed12",
        genome=config["genome"],
    output:
        temp(f"{OUT}/catalogue/merged/cohort.nt.fa"),
    container:
        IMG["bedtools"]
    shell:
        "bedtools getfasta -split -s -nameOnly -fi {input.genome} -bed {input.bed} -fo {output}"


rule orf_aa_fasta:
    """Every catalogue ORF translated, stop trimmed; headers orf_id(+/-). Used by the
    MMseqs2 collapse and for orfs.tsv's aa_seq, so it is built either way."""
    input:
        f"{OUT}/catalogue/merged/cohort.nt.fa",
    output:
        f"{OUT}/catalogue/merged/cohort.aa.fa",
    threads: 4
    container:
        IMG["seqkit"]
    shell:
        "seqkit translate --threads {threads} --trim {input} > {output}"


if not PROVIDED and str(config.get("orf_collapse", True)).lower() not in ("false", "0", "no"):

    rule mmseqs_cluster:
        input:
            f"{OUT}/catalogue/merged/cohort.aa.fa",
        output:
            f"{OUT}/catalogue/merged/mmseqs/cohort_cluster.tsv",
        params:
            prefix=f"{OUT}/catalogue/merged/mmseqs/cohort",
            tmp=f"{OUT}/catalogue/merged/mmseqs/tmp",
        threads: 8
        resources:
            mem_mb=scaled(16000),
            runtime=60,
        container:
            IMG["mmseqs2"]
        log:
            f"{OUT}/logs/mmseqs.log",
        shell:
            r"""
            mmseqs easy-cluster {input} {params.prefix} {params.tmp}                 --min-seq-id 0.9 -c 0.8 --cov-mode 0 --threads {threads} > {log} 2>&1
            rm -rf {params.tmp}
            """

    rule collapse_catalogue:
        input:
            tsv=f"{MERGED}.tsv",
            bed=f"{MERGED}.bed12",
            clusters=f"{OUT}/catalogue/merged/mmseqs/cohort_cluster.tsv",
        output:
            bed=f"{CAT}.bed12",
            tsv=f"{CAT}.tsv",
            consensus=f"{CAT}.consensus.bed12",
        params:
            opts=" ".join([_opt("--smorf-max-aa", "smorf_max_aa"),
                           _opt("--min-callers", "orf_min_callers"),
                           _opt("--min-samples", "orf_min_samples")]),
        conda:
            "../envs/score.yaml"
        log:
            f"{OUT}/logs/collapse_catalogue.log",
        shell:
            r"""
            python {SCRIPTS}/collapse_catalogue.py --tsv {input.tsv} --bed12 {input.bed}                 --clusters {input.clusters} {params.opts} --out-prefix {CAT} 2> {log}
            cat {log} >&2
            """

else:

    rule catalogue_no_collapse:
        input:
            bed=f"{MERGED}.bed12",
            tsv=f"{MERGED}.tsv",
            consensus=f"{MERGED}.consensus.bed12",
        output:
            bed=f"{CAT}.bed12",
            tsv=f"{CAT}.tsv",
            consensus=f"{CAT}.consensus.bed12",
        shell:
            "cp {input.bed} {output.bed}; cp {input.tsv} {output.tsv}; cp {input.consensus} {output.consensus}"
