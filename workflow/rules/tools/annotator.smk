# Annotator: smORF typing (uORF, uaORF, dORF, ...) against the Ensembl GTF.
# The Annotator is part of this repo (workflow/scripts/annotator/). Requires an
# ENSEMBL GTF: GENCODE silently types every ORF 'UA'.
# Contract: <OUT>/annotations/annotator.tsv

ANN = f"{OUT}/annotations/annotator"
ANN_SCRIPTS = f"{SCRIPTS}/annotator"


rule annotator_gtf:
    """Scored ORFs -> CDS-feature GTF with chrN:start-end embedded in gene_id."""
    input:
        bed=f"{OUT}/scoring/orf_candidates.bed12",
    output:
        gtf=f"{ANN}/smorfs.gtf",
        skipped=f"{ANN}/skipped_scaffolds.tsv",
        idmap=f"{ANN}/id_map.tsv",
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/annotator/gtf.log",
    shell:
        r"""
        python {ANN_SCRIPTS}/bed12_to_annotator_gtf.py \
            --bed12 {input.bed} \
            --out-gtf {output.gtf} \
            --out-skipped {output.skipped} \
            --out-idmap {output.idmap} 2> {log}
        cat {log} >&2
        """


rule annotator_intersect:
    input:
        smorf=f"{ANN}/smorfs.gtf",
        ref=config["gtf"],
    output:
        intersect=f"{ANN}/intersect.gtf",
        nonintersect=f"{ANN}/nonintersect.gtf",
    resources:
        mem_mb=32000,
        runtime=120,
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/annotator/intersect.log",
    shell:
        r"""
        bedtools intersect -wo -s -a {input.smorf} -b {input.ref} > {output.intersect} 2> {log}
        bedtools intersect -v  -s -a {input.smorf} -b {input.ref} > {output.nonintersect} 2>> {log}
        wc -l {output.intersect} {output.nonintersect} >&2
        """


rule annotator:
    input:
        orfs=f"{OUT}/scoring/orf_candidates.bed12",
        intersect=f"{ANN}/intersect.gtf",
        nonintersect=f"{ANN}/nonintersect.gtf",
        ref=config["gtf"],
        idmap=f"{ANN}/id_map.tsv",
        skipped=f"{ANN}/skipped_scaffolds.tsv",
    output:
        csv=f"{ANN}/smORF_annotation.csv",
        tsv=f"{OUT}/annotations/annotator.tsv",
    params:
        # A UCSC assembly name for the links, not a FASTA. The script defaults to hg38.
        genome=config.get("ucsc_genome", "mm39"),
    resources:
        mem_mb=16000,
        runtime=120,
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/annotator/annotate.log",
    shell:
        r"""
        python {ANN_SCRIPTS}/smorf_annotator.py \
            --intersect_output {input.intersect} \
            --non_intersect_output {input.nonintersect} \
            --ensembl_gtf {input.ref} \
            --genome {params.genome} \
            --output_file {output.csv} > {log} 2>&1
        python {ANN_SCRIPTS}/annotator_table.py \
            --orfs {input.orfs} --annotation {output.csv} \
            --idmap {input.idmap} --skipped {input.skipped} \
            --out {output.tsv} 2>> {log}
        tail -20 {log} >&2
        """
