# ShortStop: classifies smORFs as SAMs (sam_intracellular / sam_secreted,
# Swiss-Prot-like) or PRISMs, with its bundled model (trained on human).
# Installed from GitHub at a pinned commit (envs/shortstop.yaml), not vendored:
# its license is non-commercial academic.
# Contract: <OUT>/annotations/shortstop.tsv

SS = f"{OUT}/annotations/shortstop"
SS_SCRIPTS = f"{SCRIPTS}/shortstop"


rule shortstop_input:
    input:
        bed=f"{OUT}/scoring/orf_candidates.bed12",
        genome=config["genome"],
    output:
        gtf=f"{SS}/input/smorfs.gtf",
        skipped=f"{SS}/input/skipped.tsv",
        expected=f"{SS}/input/expected_protein.tsv",
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/shortstop/input.log",
    shell:
        r"""
        python {SS_SCRIPTS}/shortstop_input.py --bed12 {input.bed} --genome {input.genome} \
            --out-gtf {output.gtf} --out-skipped {output.skipped} \
            --out-expected {output.expected} 2> {log}
        cat {log} >&2
        """


rule shortstop_predict:
    """ShortStop creates --outdir with a one-level mkdir, and reads features from a
    CWD-relative default unless --orfs_to_be_predicted is given, so both are
    pinned here and it runs inside its own directory."""
    input:
        gtf=f"{SS}/input/smorfs.gtf",
        genome=config["genome"],
    output:
        pred=f"{SS}/run/predictions/shortstop_classifications.csv",
        seqs=f"{SS}/run/sequences/unknown_sequences.csv",
    params:
        run=os.path.abspath(f"{SS}/run"),
        log=os.path.abspath(f"{OUT}/logs/shortstop/predict.log"),
    resources:
        mem_mb=scaled(24000),
        runtime=480,
    conda:
        "../../envs/shortstop.yaml"
    log:
        f"{OUT}/logs/shortstop/predict.log",
    shell:
        r"""
        gtf=$(realpath {input.gtf}); genome=$(realpath {input.genome})
        log={params.log:q}; run={params.run:q}
        rm -rf "$run"; mkdir -p "$(dirname "$run")"
        cd "$(dirname "$run")"
        shortstop predict --genome "$genome" --putative_smorfs_gtf "$gtf" \
            --outdir "$run" \
            --orfs_to_be_predicted "$run/features/extracted_features_of_smorfs.csv" > "$log" 2>&1
        """


rule shortstop:
    input:
        orfs=f"{OUT}/scoring/orf_candidates.bed12",
        pred=f"{SS}/run/predictions/shortstop_classifications.csv",
        seqs=f"{SS}/run/sequences/unknown_sequences.csv",
        expected=f"{SS}/input/expected_protein.tsv",
        skipped=f"{SS}/input/skipped.tsv",
    output:
        tsv=f"{OUT}/annotations/shortstop.tsv",
        report=f"{SS}/protein_check.tsv",
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/shortstop/table.log",
    shell:
        r"""
        python {SS_SCRIPTS}/shortstop_table.py --orfs {input.orfs} \
            --predictions {input.pred} --sequences {input.seqs} \
            --expected {input.expected} --skipped {input.skipped} \
            --out {output.tsv} --report {output.report} 2> {log}
        cat {log} >&2
        """
