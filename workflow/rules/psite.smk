# plastid P-site offsets, a check of the fit, and genome-coordinate P-site tracks.
# psite runs with --aggregate: without it the median metagene was empty on every
# library tested, so every length silently got --default (see psite_check.py).


rule plastid_metagene:
    input:
        f"{IDX}/reference.canonical.gtf",
    output:
        f"{IDX}/plastid/reference.canonical_rois.txt",
    params:
        prefix=f"{IDX}/plastid/reference.canonical",
    resources:
        mem_mb=scaled(8000),
        runtime=60,
    container:
        IMG["plastid"]
    log:
        f"{IDX}/logs/plastid_metagene.log",
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        metagene generate {params.prefix} --annotation_files {input} \
            --landmark cds_start > {log} 2>&1
        """


rule plastid_psite:
    input:
        rois=f"{IDX}/plastid/reference.canonical_rois.txt",
        bam=f"{OUT}/align/{{sample}}.genome.sorted.bam",
        bai=f"{OUT}/align/{{sample}}.genome.sorted.bam.bai",
        strand=f"{OUT}/qc/strand/{{sample}}.tsv",   # gate: this library is sense-stranded
    output:
        offsets=f"{OUT}/psites/offsets/{{sample}}_p_offsets.txt",
        profiles=f"{OUT}/psites/offsets/{{sample}}_metagene_profiles.txt",
    params:
        prefix=f"{OUT}/psites/offsets/{{sample}}",
        minlen=config.get("psite_min_length", 20),
        maxlen=config.get("psite_max_length", 50),
        default=config.get("psite_default_offset", 12),
    resources:
        mem_mb=scaled(32000),   # measured peak 0.8-6.6 GB, but 29.3 GB on a 616M-read library
        runtime=60,
    container:
        IMG["plastid"]
    log:
        f"{OUT}/logs/plastid/{{sample}}_psite.log",
    shell:
        r"""
        psite {input.rois} {params.prefix} --count_files {input.bam} \
            --min_length {params.minlen} --max_length {params.maxlen} \
            --default {params.default} --require_upstream --aggregate > {log} 2>&1
        """


rule psite_check:
    """Keep only read lengths with a fitted, in-frame offset; summarize the sample."""
    input:
        offsets=f"{OUT}/psites/offsets/{{sample}}_p_offsets.txt",
        profiles=f"{OUT}/psites/offsets/{{sample}}_metagene_profiles.txt",
    output:
        offsets=f"{OUT}/psites/offsets/{{sample}}_p_offsets.kept.txt",
        lengths=f"{OUT}/psites/check/{{sample}}_lengths.tsv",
        summary=f"{OUT}/psites/check/{{sample}}.tsv",
    params:
        min_share=config.get("psite_min_share", 0.015),
        min_frame0=config.get("psite_min_frame0", 0.45),
        offset_range=" ".join(str(x) for x in config.get("psite_offset_range", [10, 16])),
        min_reads=config.get("psite_min_start_reads", 20000),
        min_kept=config.get("psite_min_kept_share", 0.5),
    conda:
        "../envs/score.yaml"
    log:
        f"{OUT}/logs/plastid/{{sample}}_psite_check.log",
    shell:
        r"""
        python {SCRIPTS}/psite_check.py --sample {wildcards.sample} \
            --offsets {input.offsets} --profiles {input.profiles} \
            --min-share {params.min_share} --min-frame0 {params.min_frame0} \
            --offset-range {params.offset_range} --min-start-reads {params.min_reads} \
            --min-kept-share {params.min_kept} --out-offsets {output.offsets} \
            --out-lengths {output.lengths} --out {output.summary} 2> {log}
        cat {log} >&2
        """


rule psite_report:
    """One table for the cohort. Fails on any FAIL sample (psite_check: warn only
    reports them), and always on a sample with no usable read length."""
    input:
        expand(f"{OUT}/psites/check/{{sample}}.tsv", sample=SAMPLES),
    output:
        f"{OUT}/psites/psite_report.tsv",
    params:
        mode=config.get("psite_check", "fail"),
    shell:
        r"""
        awk 'FNR==1 && NR!=1 {{ next }} {{ print }}' {input} > {output}
        column -t -s $'\t' {output} >&2
        bad=$(awk -F'\t' 'NR>1 && $2=="FAIL" {{ print $1 }}' {output} | paste -sd' ')
        empty=$(awk -F'\t' 'NR>1 && $5=="" {{ print $1 }}' {output} | paste -sd' ')
        if [ -n "$empty" ]; then
            echo "ERROR: no read length has a usable P-site offset in: $empty. Drop these libraries from samples.tsv." >&2
            exit 1
        fi
        if [ -n "$bad" ] && [ "{params.mode}" = "fail" ]; then
            echo "ERROR: P-site check failed for: $bad (reasons in {output}). Drop them from samples.tsv, or set psite_check: warn to score them on their kept lengths." >&2
            exit 1
        fi
        """


rule plastid_tracks:
    input:
        bam=f"{OUT}/align/{{sample}}.genome.sorted.bam",
        bai=f"{OUT}/align/{{sample}}.genome.sorted.bam.bai",
        offsets=f"{OUT}/psites/offsets/{{sample}}_p_offsets.kept.txt",
        report=f"{OUT}/psites/psite_report.tsv",   # gate: no tracks until the cohort passes
    output:
        fw=f"{OUT}/psites/{{sample}}_fw.bedgraph",
        rc=f"{OUT}/psites/{{sample}}_rc.bedgraph",
    params:
        prefix=f"{OUT}/psites/{{sample}}",
        minlen=config.get("psite_min_length", 20),
    resources:
        mem_mb=scaled(48000),   # make_wiggle measured peak 2.4-17.4 GB, 44.4 GB on a 616M-read library
        runtime=120,
    container:
        IMG["plastid"]
    log:
        f"{OUT}/logs/plastid/{{sample}}_tracks.log",
    shell:
        r"""
        make_wiggle --count_files {input.bam} --offset {input.offsets} \
            --fiveprime_variable -o {params.prefix} \
            --min_length {params.minlen} --output_format bedgraph > {log} 2>&1
        mv {params.prefix}_fw.wig {output.fw}
        mv {params.prefix}_rc.wig {output.rc}
        """


rule unique_bam:
    """Uniquely mapped reads only (STAR MAPQ 255). STAR keeps up to 20 loci per read and
    plastid counts a read at each, so a pseudogene nearly identical to its parent gets
    the parent's reads in full; the unique tracks exclude them."""
    input:
        bam=f"{OUT}/align/{{sample}}.genome.sorted.bam",
    output:
        bam=temp(f"{OUT}/psites/unique/{{sample}}.unique.bam"),
        bai=temp(f"{OUT}/psites/unique/{{sample}}.unique.bam.bai"),
    threads: 4
    resources:
        mem_mb=scaled(4000),
        runtime=60,
    container:
        IMG["star"]
    shell:
        r"""
        samtools view -@ {threads} -b -q 255 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """


rule plastid_tracks_unique:
    """P-site tracks from uniquely mapped reads, with the same fitted offsets."""
    input:
        bam=f"{OUT}/psites/unique/{{sample}}.unique.bam",
        bai=f"{OUT}/psites/unique/{{sample}}.unique.bam.bai",
        offsets=f"{OUT}/psites/offsets/{{sample}}_p_offsets.kept.txt",
    output:
        fw=f"{OUT}/psites/unique/{{sample}}_fw.bedgraph",
        rc=f"{OUT}/psites/unique/{{sample}}_rc.bedgraph",
    params:
        prefix=f"{OUT}/psites/unique/{{sample}}",
        minlen=config.get("psite_min_length", 20),
    resources:
        mem_mb=scaled(48000),
        runtime=120,
    container:
        IMG["plastid"]
    log:
        f"{OUT}/logs/plastid/{{sample}}_tracks_unique.log",
    shell:
        r"""
        make_wiggle --count_files {input.bam} --offset {input.offsets} \
            --fiveprime_variable -o {params.prefix} \
            --min_length {params.minlen} --output_format bedgraph > {log} 2>&1
        mv {params.prefix}_fw.wig {output.fw}
        mv {params.prefix}_rc.wig {output.rc}
        """
