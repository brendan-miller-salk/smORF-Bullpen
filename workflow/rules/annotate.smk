# Join every enabled annotator's per-ORF table onto the scores, as
# scoring/orf_scores_annotated.tsv (provided mode: onto the bare ORF list, as
# orfs_annotated.tsv).

ANNOTATED = f"{OUT}/orfs_annotated.tsv" if PROVIDED else f"{OUT}/scoring/orf_scores_annotated.tsv"


# Unreviewed smORFs need pooled RPKM above min_unique_rpkm outside flanked canonical
# CDS (counts.tsv) to stay in Comprehensive / Primary: counting feeds back into orfs.csv.
COUNT_GATE = (not PROVIDED and "featurecounts" in COUNTERS
              and config.get("min_unique_rpkm") is not None)


rule join_annotations:
    input:
        scores=f"{OUT}/catalogue/orf_ids.tsv" if PROVIDED else f"{OUT}/scoring/orf_scores.tsv",
        annotations=expand(f"{OUT}/annotations/{{tool}}.tsv", tool=ANNOTATORS + CONSERVATION),
    output:
        ANNOTATED,
    params:
        exclude=(" ".join(["--exclude-types", *config.get("exclude_types", [])])
                 if config.get("exclude_types") and not PROVIDED else ""),
        unique=(" ".join(["--unique-types", *config.get("unique_evidence_types", []),
                          "--unique-min-inframe", str(config.get("unique_min_inframe_psites", 5))])
                if config.get("unique_evidence_types") and not PROVIDED else ""),
        own=("--require-unshared --exempt-classes " + " ".join(EXEMPT_CLASSES)
             if config.get("require_unshared_sequence") and not PROVIDED else ""),
        chroms=(" ".join(["--exclude-chroms", *config.get("exclude_chroms", [])])
                if config.get("exclude_chroms") and not PROVIDED else ""),
    conda:
        "../envs/score.yaml"
    log:
        f"{OUT}/logs/join_annotations.log",
    shell:
        r"""
        python {SCRIPTS}/join_annotations.py --scores {input.scores} \
            --annotations {input.annotations} {params.exclude} {params.unique} {params.own} \
            {params.chroms} --out {output} 2> {log}
        cat {log} >&2
        """


if not PROVIDED:

    rule pass_all_qc:
        """pass_all vs Chothani / GENCODE v2 2025: thresholds, Primary fraction, length and overlap
        biases. Warns by default; config pass_all_qc: fail makes a WARN fatal."""
        input:
            scores=ANNOTATED,
            thresholds=f"{OUT}/scoring/thresholds.tsv",
            cds=f"{OUT}/scoring/cds_calibration.tsv",
            catalogue=f"{OUT}/catalogue/cohort.catalogue.tsv",
        output:
            f"{OUT}/scoring/pass_all_qc.tsv",
        params:
            strict="--strict" if config.get("pass_all_qc", "warn") == "fail" else "",
        conda:
            "../envs/score.yaml"
        log:
            f"{OUT}/logs/pass_all_qc.log",
        shell:
            r"""
            python {SCRIPTS}/pass_all_qc.py --scores {input.scores} --thresholds {input.thresholds} \
                --cds {input.cds} --catalogue {input.catalogue} --out {output} {params.strict} 2> {log}
            cat {log} >&2
            """


rule final_outputs:
    """The three files a reader opens: orfs.csv, summary.csv, README.md at the top
    of the outdir. Everything else stays in its subfolder."""
    input:
        table=ANNOTATED,
        catalogue=f"{OUT}/catalogue/cohort.catalogue.tsv",
        thresholds=[] if PROVIDED else f"{OUT}/scoring/thresholds.tsv",
        qc=[] if PROVIDED else f"{OUT}/scoring/pass_all_qc.tsv",
        aa=f"{OUT}/catalogue/merged/cohort.aa.fa",
        counts=f"{OUT}/counts.tsv" if COUNT_GATE else [],
    output:
        orfs=f"{OUT}/orfs.csv",
        summary=f"{OUT}/summary.csv",
        readme=f"{OUT}/README.md",
    params:
        name=f"{'Provided' if PROVIDED else 'Ribo-seq'} ORFs: {os.path.basename(OUT)}",
        scoring="" if PROVIDED else (f"--thresholds {OUT}/scoring/thresholds.tsv "
                                     f"--qc {OUT}/scoring/pass_all_qc.tsv "
                                     f"--cds-median-fdr {config.get('cds_median_fdr', 0.05)}"),
        table=os.path.relpath(ANNOTATED, OUT),
        sets=(f'"reviewed_smorfs .gtf / .nuc / .pep: {"Comprehensive " if not PROVIDED else ""}smORFs '
              f'whose protein is exactly a Swiss-Prot sequence (blastoff sp_identical): ORF structure, CDS nt incl. stop, protein" '
              f'"unreviewed_smorfs .gtf / .nuc / .pep: the same for everything else '
              f'(novel, TrEMBL incl. exact, near-identical and other non-identical Swiss-Prot)"') if "blastoff" in ANNOTATORS else "",
        max_aa=config.get("smorf_max_aa", 150),
        counts=(f"--counts {OUT}/counts.tsv --min-unique-rpkm {config['min_unique_rpkm']}"
                if COUNT_GATE else ""),
        ucsc=" ".join(
            [f"--ucsc-genome {config.get('ucsc_genome', 'mm39')}"]
            + ([f"--ucsc-session {config['ucsc_session']}"] if config.get("ucsc_session") else [])
            + ([f"--ucsc-pad {config['ucsc_pad']}"] if config.get("ucsc_pad") is not None else [])),
    conda:
        "../envs/score.yaml"
    log:
        f"{OUT}/logs/final_outputs.log",
    shell:
        r"""
        python {SCRIPTS}/final_outputs.py --table {input.table} --catalogue {input.catalogue} \
            {params.scoring} --aa-fasta {input.aa} {params.ucsc} --outdir {OUT} \
            --smorf-max-aa {params.max_aa} {params.counts} \
            --name "{params.name}" \
            --details {params.sets} "subfolders: the full per-ORF table ({params.table}) and every step's files" 2> {log}
        cat {log} >&2
        """


# Comprehensive smORFs split by UniProt match, as GTF + nucleotide + peptide FASTA.
SMORF_SETS = [f"{OUT}/{s}.{x}" for s in ("reviewed_smorfs", "unreviewed_smorfs")
              for x in ("gtf", "nuc", "pep")] if "blastoff" in ANNOTATORS else []

if SMORF_SETS:

    rule smorf_sets:
        """reviewed_smorfs (blastoff sp_identical) and unreviewed_smorfs (the rest): Comprehensive
        smORFs (provided mode: all), each as .gtf, .nuc (CDS + stop) and .pep (aa_seq)."""
        input:
            orfs=f"{OUT}/orfs.csv",
            bed12=f"{OUT}/scoring/orf_candidates.bed12",
            genome=config["genome"],
        output:
            SMORF_SETS,
        conda:
            "../envs/score.yaml"
        log:
            f"{OUT}/logs/smorf_sets.log",
        shell:
            r"""
            python {SCRIPTS}/smorf_sets.py --orfs {input.orfs} --bed12 {input.bed12} \
                --genome {input.genome} --outdir {OUT} 2> {log}
            cat {log} >&2
            """
