# Chothani translation-signature scoring and smORF typing of the catalogue,
# from this workflow's catalogue and P-site tracks.

PSITE_FW = expand(f"{OUT}/psites/{{sample}}_fw.bedgraph", sample=SAMPLES)
PSITE_RC = expand(f"{OUT}/psites/{{sample}}_rc.bedgraph", sample=SAMPLES)
PSITE_FW_UNIQUE = expand(f"{OUT}/psites/unique/{{sample}}_fw.bedgraph", sample=SAMPLES)
PSITE_RC_UNIQUE = expand(f"{OUT}/psites/unique/{{sample}}_rc.bedgraph", sample=SAMPLES)


rule score:
    """PIF / Uniformity / Drop-off, calibrated on annotated CDS, then thresholded."""
    input:
        orfs=f"{SC}/orf_candidates.bed12",
        gtf=config["gtf"],
        fwd=PSITE_FW,
        rev=PSITE_RC,
        fwd_unique=PSITE_FW_UNIQUE,
        rev_unique=PSITE_RC_UNIQUE,
    output:
        orfs=f"{SC}/orf_scores.tsv",
        cds=f"{SC}/cds_calibration.tsv",
        thresholds=f"{SC}/thresholds.tsv",
    params:
        tail=config["threshold_tail"],
        pct=config["threshold_percentile"],
        zero=config["uniformity_zero_codons"],
        method=config.get("threshold_method", "empirical"),
        min_ps=config.get("calibration_min_psites", 100),
        require=" ".join(config["require"]),
        min_ppm=(f"--min-ppm {config['min_ppm']}" if config.get("min_ppm") is not None else ""),
        unshared=config.get("unshared_evidence", "report"),
        unshared_min=config.get("unshared_min_inframe_psites", 1),
        unshared_frame=config.get("unshared_frame", "aware"),
    resources:
        mem_mb=32000,
        runtime=120,
    conda:
        "../envs/score.yaml"
    log:
        f"{OUT}/logs/score.log",
    shell:
        r"""
        python {SCRIPTS}/score_orfs.py \
            --orf-bed12 {input.orfs} \
            --gtf {input.gtf} \
            --psite-forward {input.fwd} \
            --psite-reverse {input.rev} \
            --psite-forward-unique {input.fwd_unique} \
            --psite-reverse-unique {input.rev_unique} \
            --out-orfs {output.orfs} \
            --out-cds {output.cds} \
            --out-thresholds {output.thresholds} \
            --threshold-tail {params.tail} \
            --threshold-percentile {params.pct} \
            --uniformity-zero-codons {params.zero} \
            --threshold-method {params.method} \
            --calibration-min-psites {params.min_ps} \
            --require {params.require} \
            --unshared-evidence {params.unshared} \
            --unshared-min-inframe {params.unshared_min} \
            --unshared-frame {params.unshared_frame} \
            {params.min_ppm} 2> {log}
        cat {log} >&2
        """


rule primary_set:
    """orfs.csv's Primary ORFs, as BED12 for browser loading. Read from orfs.csv, not the
    scores table, so the count gate (min_unique_rpkm, applied in final_outputs) holds."""
    input:
        orfs=f"{OUT}/orfs.csv",
        bed=f"{SC}/orf_candidates.bed12",
    output:
        f"{SC}/primary_set.bed12",
    conda:
        "../envs/score.yaml"
    shell:
        r"""
        python -c 'import pandas as pd, sys
o = pd.read_csv(sys.argv[1], usecols=["orf_id", "orf_set"], dtype=str)
print("\n".join(sorted(set(o.loc[o.orf_set == "primary", "orf_id"]))))' {input.orfs} > {output}.ids
        awk -F'\t' 'NR==FNR {{ keep[$1]; next }} $4 in keep' {output}.ids {input.bed} > {output}
        rm -f {output}.ids
        """


rule plots:
    input:
        orfs=f"{SC}/orf_scores.tsv",
        cds=f"{SC}/cds_calibration.tsv",
        thresholds=f"{SC}/thresholds.tsv",
    output:
        f"{SC}/score_distributions.png",
    conda:
        "../envs/score.yaml"
    script:
        "../scripts/plot_scores.py"
