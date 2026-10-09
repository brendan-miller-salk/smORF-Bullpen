# Strandedness gate, per library: a library's P-sites and RiboCode wait only for
# its own check; PRICE, which pools the cohort, waits for the whole report.


rule strand_counts:
    """Uniquely mapped reads on annotated CDS, counted on the sense (-s 1) and antisense (-s 2) strand."""
    input:
        bam=f"{OUT}/align/{{sample}}.genome.sorted.bam",
        gtf=REF_GTF,
    output:
        sense=f"{OUT}/qc/strand/{{sample}}_s1.txt.summary",
        anti=f"{OUT}/qc/strand/{{sample}}_s2.txt.summary",
    params:
        prefix=f"{OUT}/qc/strand/{{sample}}",
    threads: 4
    resources:
        mem_mb=scaled(8000),
        runtime=60,
    container:
        IMG["subread"]
    log:
        f"{OUT}/logs/qc/{{sample}}_strand_counts.log",
    shell:
        r"""
        : > {log}
        for s in 1 2; do
            featureCounts -T {threads} -t CDS -g gene_id -s $s -a {input.gtf} \
                -o {params.prefix}_s$s.txt {input.bam} >> {log} 2>&1
        done
        """


rule strand_check:
    """One row for the library. Under strand_check: fail, a library that isn't
    sense-stranded fails here, so nothing downstream of it runs."""
    input:
        sense=f"{OUT}/qc/strand/{{sample}}_s1.txt.summary",
        anti=f"{OUT}/qc/strand/{{sample}}_s2.txt.summary",
    output:
        f"{OUT}/qc/strand/{{sample}}.tsv",
    params:
        min_sense=config.get("strand_min_sense", 0.9),
        min_assigned=config.get("strand_min_assigned", 100000),
        strict="--strict" if config.get("strand_check", "fail") == "fail" else "",
    conda:
        "../envs/score.yaml"
    log:
        f"{OUT}/logs/qc/{{sample}}_strand_check.log",
    shell:
        r"""
        python {SCRIPTS}/strand_check.py --sense {input.sense} --antisense {input.anti} \
            --min-sense {params.min_sense} --min-assigned {params.min_assigned} \
            --out {output} {params.strict} 2> {log} || {{ cat {log} >&2; exit 1; }}
        cat {log} >&2
        """


rule strand_report:
    """One table for the cohort."""
    input:
        expand(f"{OUT}/qc/strand/{{sample}}.tsv", sample=SAMPLES),
    output:
        f"{OUT}/qc/strand_report.tsv",
    shell:
        r"""
        awk 'FNR==1 && NR!=1 {{ next }} {{ print }}' {input} > {output}
        column -t -s $'\t' {output} >&2
        """
