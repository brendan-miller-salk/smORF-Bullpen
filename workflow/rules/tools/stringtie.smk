# Assembly: StringTie 3 on bulk RNA-seq builds the run's novel transcript models, which
# rules/reference.smk then treats as novel_gtf (gffcompare -> novel_class_codes -> hybrid GTFs).
#   fastp per lane -> STAR, reference index, per-sample 2-pass -> strandedness from ReadsPerGene (auto)
#   -> stringtie -G per sample -> --merge -G -> stringtie -e per sample on the merge
#   -> novel (not a reference transcript_id) with TPM >= assembly_min_tpm in
#      >= assembly_min_samples samples -> assembly/novel.gtf
# RNA-seq samples use the wildcard `rna`: `sample` is constrained to the Ribo-seq libraries.

ASM = f"{OUT}/assembly"
RNA_STRAND = str(config.get("rnaseq_strandedness", "auto"))
if RNA_STRAND not in ("auto", "unstranded", "rf", "fr"):
    raise WorkflowError("rnaseq_strandedness must be auto, unstranded, rf (dUTP) or fr.")
if "rnaseq_samplesheet" not in config:
    raise WorkflowError("tools.assembly needs rnaseq_samplesheet (columns sample, fastq_1, fastq_2).")
for _img, _eg in (("stringtie", "stringtie:3.0.3--h29c0135_0"), ("fastp", "fastp:1.0.1--heae3180_0")):
    if _img not in IMG:
        raise WorkflowError(f"tools.assembly needs containers.{_img}, e.g. docker://quay.io/biocontainers/{_eg}")


def _read_rnaseq(path):
    """{sample: [(fastq_1, fastq_2 or ''), ...]}: one row per lane, a sample's rows concatenated."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or not {"sample", "fastq_1", "fastq_2"} <= set(rows[0]):
        raise WorkflowError(f"{path}: need tab-separated columns 'sample', 'fastq_1', 'fastq_2'.")
    lanes = {}
    for r in rows:
        s, fq1, fq2 = (str(r[k] or "").strip() for k in ("sample", "fastq_1", "fastq_2"))
        if not fq1 or any(f and not os.path.isabs(f) for f in (fq1, fq2)):
            raise WorkflowError(f"{path}: {s} needs an absolute fastq_1 (and fastq_2 if paired).")
        lanes.setdefault(s, []).append((fq1, fq2))
    for s, ls in lanes.items():
        if len({bool(fq2) for _, fq2 in ls}) > 1:
            raise WorkflowError(f"{path}: {s} mixes single- and paired-end rows.")
    return lanes


RNA_LANES = _read_rnaseq(config["rnaseq_samplesheet"])
RNA_SAMPLES = sorted(RNA_LANES)
# StringTie's flag for the strandedness rnaseq_strand settled on.
STRAND_FLAG = r"$(sed 's/^unstranded$//; s/^rf$/--rf/; s/^fr$/--fr/' {input.strand})"


wildcard_constraints:
    rna="|".join(map(re.escape, RNA_SAMPLES)),
    lane=r"\d+",


# The NSCLC isoform catalogue's trimming (isoform_catalog stage 02): both TruSeq adapters
# given (fastp would otherwise copy R1's adapter to R2), poly-G (NextSeq / NovaSeq), a
# gentle 3' quality trim, and pairs with a mate < 36 nt dropped.
FASTP_ARGS = ("--trim_poly_g --poly_g_min_len 10 --cut_tail --cut_tail_window_size 4 "
              "--cut_tail_mean_quality 20 --length_required 36")
# Adapters: auto = fastp detects them per lane (--detect_adapter_for_pe; single-end detects
# by default), recorded in assembly/adapter_report.tsv. Or give both sequences: with only
# R1's, fastp would copy it to R2.
ADAPTER_R1 = str(config.get("rnaseq_adapter_r1", "auto"))
ADAPTER_R2 = str(config.get("rnaseq_adapter_r2", "auto"))
if (ADAPTER_R1 == "auto") != (ADAPTER_R2 == "auto"):
    raise WorkflowError("rnaseq_adapter_r1 and rnaseq_adapter_r2 must both be auto or both sequences.")
ADAPT_PE = ("--detect_adapter_for_pe" if ADAPTER_R1 == "auto" else
            f"--adapter_sequence {ADAPTER_R1} --adapter_sequence_r2 {ADAPTER_R2}")
ADAPT_SE = "" if ADAPTER_R1 == "auto" else f"--adapter_sequence {ADAPTER_R1}"


def trimmed(wc, mate):
    """A sample's trimmed lanes for STAR: R1 / R2 when paired, SE otherwise."""
    lanes = RNA_LANES[wc.rna]
    paired = bool(lanes[0][1])
    if mate == 2 and not paired:
        return []
    name = ("R1" if mate == 1 else "R2") if paired else "SE"
    return [f"{ASM}/trim/{wc.rna}_L{i}_{name}.fq.gz" for i in range(len(lanes))]


# One job per lane, each reading its .gz directly: fastp stops on a corrupt gzip, but
# from a pipe (several lanes zcat'ed together) it reads nothing and still exits 0.
rule rnaseq_fastp_pe:
    input:
        r1=lambda wc: RNA_LANES[wc.rna][int(wc.lane)][0],
        r2=lambda wc: RNA_LANES[wc.rna][int(wc.lane)][1],
    output:
        r1=temp(f"{ASM}/trim/{{rna}}_L{{lane}}_R1.fq.gz"),
        r2=temp(f"{ASM}/trim/{{rna}}_L{{lane}}_R2.fq.gz"),
        json=f"{ASM}/trim/{{rna}}_L{{lane}}.fastp.json",
        html=f"{ASM}/trim/{{rna}}_L{{lane}}.fastp.html",
    params:
        args=FASTP_ARGS,
    threads: 8
    resources:
        mem_mb=scaled(8000),
        runtime=240,
    container:
        IMG["fastp"]
    log:
        f"{OUT}/logs/assembly/fastp_{{rna}}_L{{lane}}.log",
    shell:
        r"""
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            {ADAPT_PE} {params.args} \
            --thread {threads} --json {output.json} --html {output.html} \
            --report_title {wildcards.rna}_L{wildcards.lane} > {log} 2>&1
        """


rule rnaseq_fastp_se:
    input:
        lambda wc: RNA_LANES[wc.rna][int(wc.lane)][0],
    output:
        fq=temp(f"{ASM}/trim/{{rna}}_L{{lane}}_SE.fq.gz"),
        json=f"{ASM}/trim/{{rna}}_L{{lane}}_SE.fastp.json",
        html=f"{ASM}/trim/{{rna}}_L{{lane}}_SE.fastp.html",
    params:
        args=FASTP_ARGS,
    threads: 8
    resources:
        mem_mb=scaled(8000),
        runtime=240,
    container:
        IMG["fastp"]
    log:
        f"{OUT}/logs/assembly/fastp_{{rna}}_L{{lane}}_SE.log",
    shell:
        r"""
        fastp -i {input} -o {output.fq} {ADAPT_SE} {params.args} \
            --thread {threads} --json {output.json} --html {output.html} \
            --report_title {wildcards.rna}_L{wildcards.lane} > {log} 2>&1
        """


def fastp_jsons(wc=None):
    out = []
    for s, lanes in RNA_LANES.items():
        se = "" if lanes[0][1] else "_SE"
        out += [f"{ASM}/trim/{s}_L{i}{se}.fastp.json" for i in range(len(lanes))]
    return out


rule rnaseq_adapter_report:
    """The adapters fastp used per lane (detected or given) and what it trimmed."""
    input:
        fastp_jsons(),
    output:
        f"{ASM}/adapter_report.tsv",
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/assembly/adapter_report.log",
    shell:
        r"""
        python {SCRIPTS}/fastp_report.py --json {input} --out {output} 2> {log}
        cat {log} >&2
        """


rule rnaseq_star:
    """ENCODE long-RNA settings. intronMotif gives unstranded spliced reads the XS tag
    StringTie reads their strand from; alignSJDBoverhangMin 3 because pass-1 junctions
    count as annotated in pass 2."""
    input:
        r1=lambda wc: trimmed(wc, 1),
        r2=lambda wc: trimmed(wc, 2),
        idx=f"{IDX}/star_reference",
    output:
        bam=temp(f"{ASM}/align/{{rna}}/Aligned.out.bam"),
        counts=f"{ASM}/align/{{rna}}/ReadsPerGene.out.tab",
        log=f"{ASM}/align/{{rna}}/Log.final.out",
    params:
        prefix=f"{ASM}/align/{{rna}}/",
        reads=lambda wc, input: " ".join(",".join(m) for m in (input.r1, input.r2) if m),
        extra=config.get("rnaseq_star_extra_args", ""),
    threads: 32
    resources:
        mem_mb=scaled(56000),
        runtime=360,
    container:
        IMG["star"]
    log:
        f"{OUT}/logs/assembly/star_{{rna}}.log",
    shell:
        r"""
        STAR --genomeDir {input.idx} --readFilesIn {params.reads} --readFilesCommand zcat \
            --runThreadN {threads} --outFileNamePrefix {params.prefix} \
            --outSAMattrRGline ID:{wildcards.rna} SM:{wildcards.rna} \
            --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --outFilterMultimapNmax 20 --outFilterMismatchNoverReadLmax 0.04 \
            --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 \
            --outSAMtype BAM Unsorted --quantMode GeneCounts \
            --runRNGseed 0 --twopassMode Basic {params.extra} > {log} 2>&1
        """


rule rnaseq_sort:
    input:
        f"{ASM}/align/{{rna}}/Aligned.out.bam",
    output:
        bam=f"{ASM}/align/{{rna}}.sorted.bam",
        bai=f"{ASM}/align/{{rna}}.sorted.bam.bai",
    threads: 6
    resources:
        mem_mb=scaled(16000),
        runtime=240,
    container:
        IMG["star"]
    shell:
        r"""
        samtools sort -T {output.bam}.tmp --threads {threads} -m 1500M -o {output.bam} {input}
        samtools index -@ {threads} {output.bam}
        """


rule rnaseq_strand:
    """One strandedness for the cohort (rnaseq_strandedness: auto), else the configured one."""
    input:
        expand(f"{ASM}/align/{{rna}}/ReadsPerGene.out.tab", rna=RNA_SAMPLES),
    output:
        report=f"{ASM}/strand_report.tsv",
        flag=f"{ASM}/strandedness.txt",
    params:
        mode=RNA_STRAND,
        # 'x' in the class codes: every library must measure stranded (rules/common.smk)
        require="--require-stranded" if NEED_STRANDED else "",
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/assembly/rnaseq_strand.log",
    shell:
        r"""
        python {SCRIPTS}/rnaseq_strand.py --counts {input} --mode {params.mode} {params.require} \
            --out {output.report} --out-flag {output.flag} 2> {log}
        cat {log} >&2
        """


rule stringtie_assemble:
    """Reference-guided, per sample. exclude_chroms (chrM) is left out of the assembly."""
    input:
        bam=f"{ASM}/align/{{rna}}.sorted.bam",
        bai=f"{ASM}/align/{{rna}}.sorted.bam.bai",
        gtf=REF_GTF,
        strand=f"{ASM}/strandedness.txt",
    output:
        f"{ASM}/stringtie/{{rna}}.gtf",
    params:
        exclude=(f"-x {','.join(config['exclude_chroms'])}" if config.get("exclude_chroms") else ""),
    threads: 8
    resources:
        mem_mb=scaled(16000),
        runtime=240,
    container:
        IMG["stringtie"]
    log:
        f"{OUT}/logs/assembly/stringtie_{{rna}}.log",
    shell:
        r"""
        stringtie {input.bam} -G {input.gtf} """ + STRAND_FLAG + r""" {params.exclude} \
            -p {threads} -l {wildcards.rna} -o {output} > {log} 2>&1
        """


rule stringtie_merge:
    """The reference enters as guides only: inputs need TPM/coverage values."""
    input:
        gtfs=expand(f"{ASM}/stringtie/{{rna}}.gtf", rna=RNA_SAMPLES),
        ref=REF_GTF,
    output:
        f"{ASM}/merged.gtf",
    params:
        inputs=f"{ASM}/merge_inputs.txt",
        args=config.get("assembly_merge_args", "-F 0 -T 1"),
    threads: 8
    resources:
        mem_mb=scaled(16000),
        runtime=120,
    container:
        IMG["stringtie"]
    log:
        f"{OUT}/logs/assembly/stringtie_merge.log",
    shell:
        r"""
        printf '%s\n' {input.gtfs} > {params.inputs}
        stringtie --merge -G {input.ref} {params.args} -p {threads} -o {output} {params.inputs} > {log} 2>&1
        """


rule stringtie_quant:
    """Each sample's TPM for every merged transcript (-e: no new models)."""
    input:
        bam=f"{ASM}/align/{{rna}}.sorted.bam",
        bai=f"{ASM}/align/{{rna}}.sorted.bam.bai",
        gtf=f"{ASM}/merged.gtf",
        strand=f"{ASM}/strandedness.txt",
    output:
        f"{ASM}/expr/{{rna}}.gtf",
    threads: 8
    resources:
        mem_mb=scaled(16000),
        runtime=240,
    container:
        IMG["stringtie"]
    log:
        f"{OUT}/logs/assembly/stringtie_quant_{{rna}}.log",
    shell:
        r"""
        stringtie {input.bam} -e -G {input.gtf} """ + STRAND_FLAG + r""" \
            -p {threads} -o {output} > {log} 2>&1
        """


rule assembly_filter:
    input:
        merged=f"{ASM}/merged.gtf",
        ref=REF_GTF,
        expr=expand(f"{ASM}/expr/{{rna}}.gtf", rna=RNA_SAMPLES),
    output:
        gtf=f"{ASM}/novel.gtf",
        tpm=f"{ASM}/novel_tpm.tsv",
        report=f"{ASM}/assembly_report.tsv",
    params:
        min_tpm=config.get("assembly_min_tpm", 1),
        min_samples=config.get("assembly_min_samples", 2),
    resources:
        mem_mb=scaled(16000),
        runtime=60,
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/assembly/assembly_filter.log",
    shell:
        r"""
        python {SCRIPTS}/assembly_filter.py --merged {input.merged} --reference {input.ref} \
            --expr {input.expr} --min-tpm {params.min_tpm} --min-samples {params.min_samples} \
            --out-gtf {output.gtf} --out-tpm {output.tpm} --out-report {output.report} 2> {log}
        cat {log} >&2
        """
