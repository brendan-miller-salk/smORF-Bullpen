ADAPTER = str(config.get("adapter", "auto")).strip()


rule trim:
    """`adapter: auto` lets Trim Galore detect Illumina / Nextera / small RNA;
    a sequence is passed as --adapter. adapter_check then verifies either one."""
    input:
        lambda wc: SAMPLE_FASTQ[wc.sample],
    output:
        fq=temp(f"{OUT}/preprocess/trim/{{sample}}_trimmed.fq.gz"),
        report=f"{OUT}/preprocess/trim/{{sample}}_trimming_report.txt",
        json=f"{OUT}/preprocess/trim/{{sample}}_trimming_report.json",
        stdout=f"{OUT}/preprocess/trim/{{sample}}.trim_galore.log",
    params:
        dir=f"{OUT}/preprocess/trim",
        adapter="" if ADAPTER.lower() == "auto" else f"--adapter {ADAPTER}",
        minlen=config.get("min_read_length", 20),
    threads: 8
    resources:
        mem_mb=scaled(8000),
        runtime=240,
    container:
        IMG["trimgalore"]
    shell:
        r"""
        trim_galore {params.adapter} --length {params.minlen} \
            --cores $(( {threads} > 3 ? {threads} - 3 : 1 )) --gzip \
            --basename {wildcards.sample} -o {params.dir} {input} > {output.stdout} 2>&1
        mv {params.dir}/$(basename {input})_trimming_report.txt {output.report}
        mv {params.dir}/$(basename {input})_trimming_report.json {output.json}
        """


rule adapter_check:
    """Fail the sample unless the adapter used is in >= adapter_min_fraction of reads."""
    input:
        fastq=lambda wc: SAMPLE_FASTQ[wc.sample],
        json=f"{OUT}/preprocess/trim/{{sample}}_trimming_report.json",
        log=f"{OUT}/preprocess/trim/{{sample}}.trim_galore.log",
    output:
        f"{OUT}/preprocess/adapter_check/{{sample}}.tsv",
    params:
        mode=ADAPTER,
        min_frac=config.get("adapter_min_fraction", 0.5),
    conda:
        "../envs/score.yaml"
    log:
        f"{OUT}/logs/adapter_check/{{sample}}.log",
    shell:
        r"""
        python {SCRIPTS}/adapter_check.py --sample {wildcards.sample} --fastq {input.fastq} \
            --json {input.json} --log {input.log} --mode {params.mode} \
            --min-fraction {params.min_frac} --out {output} 2> {log}
        cat {log} >&2
        """


rule adapter_report:
    """One table for the cohort; fails if samples were trimmed with different adapters."""
    input:
        expand(f"{OUT}/preprocess/adapter_check/{{sample}}.tsv", sample=SAMPLES),
    output:
        f"{OUT}/preprocess/adapter_report.tsv",
    shell:
        r"""
        awk 'FNR==1 && NR!=1 {{ next }} {{ print }}' {input} > {output}
        n=$(awk -F'\t' 'NR>1 {{ print $3 }}' {output} | sort -u | wc -l)
        column -t -s $'\t' {output} >&2
        if [ "$n" -gt 1 ]; then
            echo "ERROR: samples were trimmed with $n different adapters; set adapter: explicitly." >&2
            exit 1
        fi
        """


rule sortmerna:
    """Remove rRNA. Each sample gets its own workdir: SortMeRNA refuses to start
    on an existing kvdb."""
    input:
        fq=f"{OUT}/preprocess/trim/{{sample}}_trimmed.fq.gz",
        adapter_ok=f"{OUT}/preprocess/adapter_check/{{sample}}.tsv",
        idx=f"{IDX}/sortmerna/idx",
        dbs=expand(f"{IDX}/sortmerna/db/{{db}}.fasta", db=config["sortmerna_dbs"]),
    output:
        fq=f"{OUT}/preprocess/sortmerna/{{sample}}.non_rRNA.fq.gz",
        log=f"{OUT}/preprocess/sortmerna/{{sample}}.sortmerna.log",
    params:
        refs=lambda wc, input: " ".join(f"--ref {f}" for f in input.dbs),
        wd=f"{OUT}/preprocess/sortmerna/work_{{sample}}",
    threads: 32
    resources:
        # Measured on 8 quadriceps libraries: peak 17-58 GB;
        # 1h47-2h24 at cohort concurrency, but one library took 7 h.
        mem_mb=scaled(64000),
        runtime=600,
    container:
        IMG["sortmerna"]
    log:
        f"{OUT}/logs/sortmerna/{{sample}}.log",
    shell:
        r"""
        rm -rf {params.wd}; mkdir -p {params.wd}
        sortmerna {params.refs} --idx-dir {input.idx} --index 0 \
            --reads {input.fq} --workdir {params.wd} --threads {threads} \
            --aligned {params.wd}/rRNA_reads --other {params.wd}/non_rRNA_reads \
            --fastx --num_alignments 1 -v > {log} 2>&1
        other=$(ls {params.wd}/non_rRNA_reads*.f*q* | head -1)
        case "$other" in
            *.gz) mv "$other" {output.fq} ;;
            *)    gzip -c "$other" > {output.fq} ;;
        esac
        mv {params.wd}/rRNA_reads.log {output.log}
        rm -rf {params.wd}
        """
