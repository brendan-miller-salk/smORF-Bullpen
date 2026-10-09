# Caller: PRICE (GEDI). One job on the pooled cohort; memory grows with
# library count. Contract: <OUT>/callers/price/allsamples.orfs.tsv


rule gedi_index:
    """PRICE's genome index. The .oml holds absolute paths, so it is built in
    place and never moved."""
    input:
        fasta=config["genome"],
        gtf=FULL_HYBRID,
    output:
        oml=f"{IDX}/gedi/full_hybrid_reference.oml",
    params:
        dir=f"{IDX}/gedi",
        # Absolute paths from Python: container realpath may lack -m.
        oml=lambda wc, output: os.path.abspath(output.oml),
    threads: 4
    resources:
        mem_mb=scaled(12000),
        runtime=60,
    container:
        IMG["gedi_index"]
    log:
        f"{IDX}/logs/gedi_index.log",
    shell:
        r"""
        mkdir -p {params.dir}
        # IndexGenome keeps an existing index: after a GTF change it rewrote only the
        # .oml and PRICE ran on the old transcripts. Start from an empty index.
        rm -f {params.dir}/full_hybrid_reference.* {params.dir}/genome.fi
        ln -sf $(realpath {input.fasta}) {params.dir}/genome.fa
        gedi -e IndexGenome -s {params.dir}/genome.fa -a $(realpath {input.gtf}) \
            -n full_hybrid_reference -f {params.dir} -o {params.oml:q} -nomapping -p > {log} 2>&1
        """


rule price:
    input:
        bams=expand(f"{OUT}/align/{{sample}}.genome.sorted.bam", sample=SAMPLES),
        bais=expand(f"{OUT}/align/{{sample}}.genome.sorted.bam.bai", sample=SAMPLES),
        oml=f"{IDX}/gedi/full_hybrid_reference.oml",
        strand=f"{OUT}/qc/strand_report.tsv",   # gate: every library is sense-stranded
    output:
        tsv=f"{OUT}/callers/price/allsamples.orfs.tsv",
    params:
        dir=f"{OUT}/callers/price",
    threads: 16
    resources:
        # 7.5 GB / 56 min on one library; memory grows with library count.
        mem_mb=lambda wc, attempt: max(16000, 6000 * len(SAMPLES)) * attempt,
        runtime=600,
    container:
        IMG["gedi_price"]
    log:
        f"{OUT}/logs/price.log",
    shell:
        r"""
        mkdir -p {params.dir}
        printf '%s\n' $(realpath {input.bams}) > {params.dir}/input.bamlist
        # GEDI refuses to overwrite .cit files, so a rerun fails on the last run's;
        # a stale allsamples.* would also carry the previous run's sample set.
        rm -f {params.dir}/input.bamlist.cit {params.dir}/allsamples.* $(printf '%s.cit ' $(realpath {input.bams}))
        bamlist2cit -n {threads} -p {params.dir}/input.bamlist > {log} 2>&1
        gedi -e Price -reads {params.dir}/input.bamlist.cit -genomic {input.oml} \
            -prefix {params.dir}/allsamples -nthreads {threads} >> {log} 2>&1
        """
