# diamond: conservation of each ORF's protein in other species' transcripts.
# DIAMOND has no tblastn, so the search runs the other way round: the ORF proteins
# are the database and each species' Ensembl cDNA + ncRNA are the blastx queries.
# Same hits as tblastn of the ORFs against the transcripts. Only forward frames
# count (Ensembl transcripts are sense). Per ORF and species the best hit is kept;
# the ORF is conserved there if it has E <= diamond_evalue and covers at least
# diamond_min_cov of the protein. E is recomputed from the bitscore as tblastn would
# report it (ORF as query), so it does not depend on the run's ORF count; DIAMOND's
# own E-value would. Two databases: the whole transcriptome, and (with
# diamond_orthologs) the transcripts of the host gene's Ensembl orthologs, where a
# short ORF can still reach significance. ORFs that even a perfect match could not
# make significant are called untestable, not not_conserved.
# Needs config diamond_species (written by refs/fetch_conservation.sh).
# Contract: <OUT>/annotations/diamond.tsv

DM = f"{OUT}/annotations/diamond"
DM_SPECIES = config.get("diamond_species") or {}
DM_ORTH = config.get("diamond_orthologs") or {}
DIAMOND_FMT = "6 qseqid sseqid pident length qstart qend qlen qframe sstart send slen evalue bitscore"

if not DM_SPECIES:
    raise WorkflowError("tools.conservation has 'diamond' but config has no diamond_species "
                        "(name: transcript FASTA; refs/fetch_conservation.sh writes them).")
_bad = [s for s in DM_SPECIES if not re.fullmatch(r"[a-z][a-z0-9_]*", s)]
if _bad:
    raise WorkflowError(f"diamond_species names {_bad} must be lowercase [a-z0-9_]: they become columns.")
_extra = sorted(set(DM_ORTH) - set(DM_SPECIES))
if _extra:
    raise WorkflowError(f"diamond_orthologs names {_extra} are not in diamond_species.")
if "diamond" not in IMG:
    raise WorkflowError("tools.conservation has 'diamond' but containers has no 'diamond' "
                        "(template: docker://quay.io/biocontainers/diamond:2.2.8--he361c42_0).")


wildcard_constraints:
    species="|".join(map(re.escape, DM_SPECIES)),


rule diamond_query:
    """The ORF proteins, translated exactly as blastoff does (internal stops skipped)."""
    input:
        bed=f"{OUT}/scoring/orf_candidates.bed12",
        genome=config["genome"],
    output:
        f"{DM}/orfs.faa",
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/diamond/query.log",
    shell:
        r"""
        python {SCRIPTS}/blastoff/blastoff.py query --bed12 {input.bed} --genome {input.genome} \
            --short-aa 0 --out-short /dev/null --out-long {output} 2> {log}
        cat {log} >&2
        """


rule diamond_db:
    input:
        f"{DM}/orfs.faa",
    output:
        f"{DM}/db/orfs.dmnd",
    params:
        prefix=f"{DM}/db/orfs",
    container:
        IMG["diamond"]
    log:
        f"{OUT}/logs/diamond/db.log",
    shell:
        "diamond makedb --in {input} -d {params.prefix} > {log} 2>&1"


rule diamond_blastx:
    """--min-orf 1 and --masking 0: smORFs are short, and DIAMOND's defaults drop
    short translated stretches and mask low-complexity ones. -k 0 keeps every ORF a
    transcript hits (a CDS matches its isoforms and truncations too), so the output is
    large: MEASURED ~75 hits per human transcript vs muscle's 60k ORFs, gzipped here."""
    input:
        query=lambda wc: DM_SPECIES[wc.species],
        db=f"{DM}/db/orfs.dmnd",
    output:
        f"{DM}/{{species}}.tsv.gz",
    params:
        fmt=DIAMOND_FMT,
        bits=config.get("diamond_search_bits", 30),
    threads: 16
    resources:
        mem_mb=scaled(32000),
        # 6 h, doubled on the profile's retry. MEASURED: 3,000 human transcripts in 137 s
        # on 4 cpu, i.e. ~1 h per species on 16.
        runtime=lambda wc, attempt: 360 * attempt,
    container:
        IMG["diamond"]
    log:
        f"{OUT}/logs/diamond/{{species}}.log",
    shell:
        r"""
        diamond blastx -q {input.query} -d {input.db} -o {output} --outfmt {params.fmt} \
            --ultra-sensitive --min-orf 1 --masking 0 --max-target-seqs 0 --compress 1 \
            --min-score {params.bits} --threads {threads} > {log} 2>&1
        """


rule diamond:
    input:
        hits=expand(f"{DM}/{{species}}.tsv.gz", species=list(DM_SPECIES)),
        transcripts=list(DM_SPECIES.values()),
        orthologs=list(DM_ORTH.values()),
        proteins=f"{DM}/orfs.faa",
        catalogue=f"{OUT}/catalogue/cohort.catalogue.tsv",
        annotator=[f"{OUT}/annotations/annotator.tsv"] if "annotator" in ANNOTATORS else [],
        bed=f"{OUT}/scoring/orf_candidates.bed12",
    output:
        f"{OUT}/annotations/diamond.tsv",
    params:
        hits=" ".join(f"{s}={DM}/{s}.tsv.gz" for s in DM_SPECIES),
        transcripts=" ".join(f"{s}={p}" for s, p in DM_SPECIES.items()),
        orthologs=" ".join(f"{s}={p}" for s, p in DM_ORTH.items()),
        annotator=f"--annotator {OUT}/annotations/annotator.tsv" if "annotator" in ANNOTATORS else "",
        evalue=config.get("diamond_evalue", 1e-3),
        cov=config.get("diamond_min_cov", 0.5),
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/diamond/classify.log",
    shell:
        r"""
        python {SCRIPTS}/diamond_classify.py --hits {params.hits} \
            --transcripts {params.transcripts} --orthologs {params.orthologs} \
            --proteins {input.proteins} --catalogue {input.catalogue} {params.annotator} --bed12 {input.bed} \
            --evalue {params.evalue} --min-cov {params.cov} --out {output} 2> {log}
        cat {log} >&2
        """
