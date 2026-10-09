# Caller: RiboCode, per sample on transcriptome BAMs.
# Contract: <OUT>/callers/ribocode/<sample>.txt (empty if no periodic read lengths)


rule ribocode_annotation:
    input:
        fasta=config["genome"],
        gtf=HYBRID,
    output:
        gtf=f"{IDX}/ribocode/hybrid_updated.gtf",
        annot=directory(f"{IDX}/ribocode/annotation"),
    threads: 1
    resources:
        mem_mb=scaled(8000),
        runtime=60,
    container:
        IMG["ribocode"]
    log:
        f"{IDX}/logs/ribocode_annotation.log",
    shell:
        r"""
        GTFupdate {input.gtf} > {output.gtf} 2> {log}
        prepare_transcripts -g {output.gtf} -f {input.fasta} -o {output.annot} >> {log} 2>&1
        python -c "from RiboCode.prepare_transcripts import GenomeSeq; GenomeSeq('{output.annot}/transcripts_sequence.fa')" >> {log} 2>&1
        """


rule ribocode_metaplots:
    input:
        annot=f"{IDX}/ribocode/annotation",
        # The hybrid alignment's transcriptome BAM straight from STAR.
        bam=f"{OUT}/align/hybrid/{{sample}}/Aligned.toTranscriptome.out.bam",
        strand=f"{OUT}/qc/strand/{{sample}}.tsv",   # gate: this library is sense-stranded
    output:
        f"{OUT}/callers/ribocode/{{sample}}_pre_config.txt",
    params:
        prefix=f"{OUT}/callers/ribocode/{{sample}}",
    resources:
        mem_mb=scaled(8000),
        runtime=60,
    container:
        IMG["ribocode"]
    log:
        f"{OUT}/logs/ribocode/{{sample}}_metaplots.log",
    shell:
        r"""
        # The config records the BAM path as given; RiboCode re-reads it from there.
        metaplots -a {input.annot} -r $(realpath {input.bam}) -o {params.prefix} > {log} 2>&1
        """


rule ribocode:
    """An empty metaplots config means no read length showed 3-nt periodicity;
    the sample then contributes an empty table rather than failing the run.

    RiboCode caches P-sites as <CWD>/<config sample name>_psites.hd5 and reuses
    that file if it exists. The sample name comes from the BAM basename, which
    is `Aligned.toTranscriptome.out` for every sample, so a shared CWD makes every
    sample silently reuse the first one's P-sites. Each run gets its own CWD."""
    input:
        annot=f"{IDX}/ribocode/annotation",
        config=f"{OUT}/callers/ribocode/{{sample}}_pre_config.txt",
    output:
        f"{OUT}/callers/ribocode/{{sample}}.txt",
    params:
        prefix=lambda wc: os.path.abspath(f"{OUT}/callers/ribocode/{wc.sample}"),
        workdir=lambda wc: os.path.abspath(f"{OUT}/callers/ribocode/work_{wc.sample}"),
        annot=lambda wc, input: os.path.abspath(input.annot),
        config=lambda wc, input: os.path.abspath(input.config),
        log=lambda wc: os.path.abspath(f"{OUT}/logs/ribocode/{wc.sample}.log"),
    resources:
        mem_mb=scaled(16000),
        runtime=120,
    container:
        IMG["ribocode"]
    log:
        f"{OUT}/logs/ribocode/{{sample}}.log",
    shell:
        r"""
        if ! grep -qvE '^(#|[[:space:]]*$)' {input.config}; then
            echo "WARNING: {wildcards.sample}: no periodic read lengths; RiboCode skipped" | tee {log} >&2
            : > {output}
            exit 0
        fi
        rm -rf {params.workdir}; mkdir -p {params.workdir}; cd {params.workdir}
        RiboCode -a {params.annot:q} -c {params.config:q} -o {params.prefix:q} 2>&1 | tee {params.log:q}
        if grep -q Error {params.log:q}; then exit 1; fi
        rm -rf {params.workdir}
        """
