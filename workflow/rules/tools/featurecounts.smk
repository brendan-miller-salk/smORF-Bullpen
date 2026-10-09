# Counter: featureCounts on canonical CDS plus the unique parts of novel smORFs.
# Novel = Comprehensive with blastoff_class 'novel' (no Swiss-Prot/TrEMBL match),
# at most smorf_max_aa long.
# Only libraries called forward in qc/strand_report.tsv are counted.

COUNT = f"{OUT}/counting"
FLANK_NT = config.get("counting_flank_nt", "auto")


def forward_samples(report):
    """Shell snippet: the samples strand_report.tsv calls forward."""
    return f"awk -F'\\t' 'NR>1 && $3==\"forward\" {{ print $1 }}' {report}"


rule count_flank_windows:
    """Canonical transcripts with long UTRs and no neighbouring CDS, as positions around start and stop."""
    input:
        canonical=f"{IDX}/reference.canonical.gtf",
        gtf=REF_GTF,
    output:
        map=f"{COUNT}/flank/windows.tsv",
        plus=f"{COUNT}/flank/windows_plus.bed",
        minus=f"{COUNT}/flank/windows_minus.bed",
    resources:
        mem_mb=scaled(16000),
        runtime=60,
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/counting/flank_windows.log",
    shell:
        r"""
        python {SCRIPTS}/count_flank.py windows --canonical-gtf {input.canonical} --gtf {input.gtf} \
            --out-map {output.map} --out-bed-plus {output.plus} --out-bed-minus {output.minus} 2> {log}
        """


rule count_flank_depth:
    """Whole-read sense coverage of uniquely mapped reads (MAPQ 255) in those windows."""
    input:
        bam=f"{OUT}/align/{{sample}}.genome.sorted.bam",
        plus=f"{COUNT}/flank/windows_plus.bed",
        minus=f"{COUNT}/flank/windows_minus.bed",
    output:
        plus=f"{COUNT}/flank/{{sample}}_plus.depth",
        minus=f"{COUNT}/flank/{{sample}}_minus.depth",
    resources:
        mem_mb=scaled(4000),
        runtime=120,
    container:
        IMG["star"]
    shell:
        r"""
        samtools view -u -q 255 -F 16 {input.bam} | samtools depth -a -b {input.plus} - > {output.plus}
        samtools view -u -q 255 -f 16 {input.bam} | samtools depth -a -b {input.minus} - > {output.minus}
        """


rule count_flank_fit:
    """5' and 3' counting flanks from the forward libraries' pooled metagene."""
    input:
        map=f"{COUNT}/flank/windows.tsv",
        plus=expand(f"{COUNT}/flank/{{sample}}_plus.depth", sample=SAMPLES),
        minus=expand(f"{COUNT}/flank/{{sample}}_minus.depth", sample=SAMPLES),
        strand=f"{OUT}/qc/strand_report.tsv",
    output:
        fit=f"{COUNT}/flank_fit.tsv",
        profile=f"{COUNT}/flank_profile.tsv",
        png=f"{COUNT}/flank_metagene.png",
    params:
        dir=f"{COUNT}/flank",
        fwd=forward_samples(f"{OUT}/qc/strand_report.tsv"),
        frac=config.get("counting_flank_excess_frac", 0.05),
        baseline=" ".join(str(x) for x in config.get("counting_flank_baseline", [40, 100])),
    resources:
        mem_mb=scaled(32000),
        runtime=60,
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/counting/flank_fit.log",
    shell:
        r"""
        fwd=$({params.fwd})
        python {SCRIPTS}/count_flank.py fit --map {input.map} \
            --depth-plus $(for s in $fwd; do printf '%s ' {params.dir}/${{s}}_plus.depth; done) \
            --depth-minus $(for s in $fwd; do printf '%s ' {params.dir}/${{s}}_minus.depth; done) \
            --excess-frac {params.frac} --baseline {params.baseline} \
            --out {output.fit} --out-profile {output.profile} --plot {output.png} 2> {log}
        cat {log} >&2
        """


rule count_gtf:
    """Canonical CDS + novel smORF CDS (stop trimmed) minus flanked canonical CDS."""
    input:
        scores=ANNOTATED,
        bed12=f"{OUT}/scoring/orf_candidates.bed12",
        gtf=REF_GTF,
        fit=f"{COUNT}/flank_fit.tsv",
    output:
        gtf=f"{COUNT}/count_annotation.gtf",
        orfs=f"{COUNT}/orf_unique_nt.tsv",
    params:
        override="" if str(FLANK_NT) == "auto" else f"--flank-5p {int(FLANK_NT)} --flank-3p {int(FLANK_NT)}",
        max_aa=config.get("smorf_max_aa", 150),
        exempt=" ".join(EXEMPT_CLASSES),
    resources:
        mem_mb=scaled(8000),
        runtime=30,
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/counting/count_gtf.log",
    shell:
        r"""
        python {SCRIPTS}/build_count_gtf.py --scores {input.scores} --bed12 {input.bed12} \
            --gtf {input.gtf} --flank {input.fit} {params.override} \
            --max-aa {params.max_aa} --exempt-classes {params.exempt} \
            --out-gtf {output.gtf} --out-orfs {output.orfs} 2> {log}
        cat {log} >&2
        """


rule featurecounts:
    """One pass over every forward library: -O (a read on two features counts for both),
    unique reads only (no -M), sense strand (-s 1, guaranteed by the strand check)."""
    input:
        gtf=f"{COUNT}/count_annotation.gtf",
        bams=expand(f"{OUT}/align/{{sample}}.genome.sorted.bam", sample=SAMPLES),
        strand=f"{OUT}/qc/strand_report.tsv",
    output:
        f"{COUNT}/featurecounts.txt",
    params:
        fwd=forward_samples(f"{OUT}/qc/strand_report.tsv"),
    threads: 8
    resources:
        mem_mb=scaled(8000),
        runtime=120,
    container:
        IMG["subread"]
    log:
        f"{OUT}/logs/counting/featurecounts.log",
    shell:
        r"""
        fwd=$({params.fwd})
        [ -n "$fwd" ] || {{ echo "featurecounts: no forward-stranded library to count" >&2; exit 1; }}
        featureCounts -T {threads} -O -t CDS -g gene_id -s 1 -a {input.gtf} -o {output} \
            $(for s in $fwd; do printf '%s ' {OUT}/align/$s.genome.sorted.bam; done) > {log} 2>&1
        """


rule counts_table:
    """counts.tsv: canonical genes and novel smORFs by sample."""
    input:
        counts=f"{COUNT}/featurecounts.txt",
        orfs=f"{COUNT}/orf_unique_nt.tsv",
    output:
        f"{OUT}/counts.tsv",
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/counting/counts_table.log",
    shell:
        r"""
        python {SCRIPTS}/counts_table.py --counts {input.counts} --orfs {input.orfs} \
            --out {output} 2> {log}
        cat {log} >&2
        """
