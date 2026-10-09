# Settle the catalogue's stop-codon convention before anything reads it. Shared
# by both modes: scoring and every annotator read scoring/orf_candidates.bed12.

SC = f"{OUT}/scoring"


rule normalize:
    """Settle the stop-codon convention and collapse duplicate block structures.
    Aborts (exit 2) on a catalogue that mixes conventions."""
    input:
        bed=f"{OUT}/catalogue/cohort.catalogue.bed12",
        genome=config["genome"],
    output:
        bed=f"{SC}/orf_candidates.bed12",
        report=f"{SC}/orf_validation.tsv",
    params:
        # NB: 'extend' alone is a Snakemake-reserved params name (>=8).
        extend_arg=f"--extend-stop {config.get('extend_stop', 'auto')}",
        dedup="" if config.get("dedup", True) else "--no-dedup",
    resources:
        mem_mb=16000,
        runtime=60,
    conda:
        "../envs/score.yaml"
    log:
        f"{OUT}/logs/normalize.log",
    shell:
        r"""
        python {SCRIPTS}/normalize_orf_bed12.py \
            --bed12 {input.bed} --genome {input.genome} \
            {params.extend_arg} {params.dedup} \
            --report {output.report} \
            --out-bed12 {output.bed} 2> {log}
        cat {log} >&2
        """
