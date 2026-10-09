# STAR: the same reads align twice, with identical
# arguments, to the `reference` index (filtered Ensembl) and the `hybrid` index
# (backbone + novel transcripts).
#   reference genome BAM, sorted  -> plastid P-sites, PRICE
#   hybrid transcriptome BAM, as STAR wrote it (unsorted) -> RiboCode


rule star_align:
    input:
        fq=f"{OUT}/preprocess/sortmerna/{{sample}}.non_rRNA.fq.gz",
        idx=f"{IDX}/star_{{ref}}",
        gtf=lambda wc: REF_GTF if wc.ref == "reference" else HYBRID,
    output:
        genome=temp(f"{OUT}/align/{{ref}}/{{sample}}/Aligned.out.bam"),
        tx=f"{OUT}/align/{{ref}}/{{sample}}/Aligned.toTranscriptome.out.bam",
        log=f"{OUT}/align/{{ref}}/{{sample}}/Log.final.out",
    wildcard_constraints:
        ref="reference|hybrid",
    params:
        prefix=f"{OUT}/align/{{ref}}/{{sample}}/",
        extra=config.get("star_extra_args", ""),
    threads: 32
    resources:
        mem_mb=scaled(56000),
        runtime=120,
    container:
        IMG["star"]
    log:
        f"{OUT}/logs/star/{{ref}}_{{sample}}.log",
    shell:
        r"""
        STAR --genomeDir {input.idx} --readFilesIn {input.fq} --readFilesCommand zcat \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} --outFileNamePrefix {params.prefix} \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
            --alignSJDBoverhangMin 1 --alignEndsType EndToEnd \
            --outFilterMultimapNmax 20 --outSAMattributes All \
            --outSAMstrandField intronMotif --outSAMtype BAM Unsorted \
            --quantMode TranscriptomeSAM --quantTranscriptomeSAMoutput BanSingleEnd \
            --runRNGseed 0 --twopassMode Basic \
            {params.extra} > {log} 2>&1
        """


rule sort_genome_bam:
    """The reference alignment's genome BAM, sorted and indexed."""
    input:
        f"{OUT}/align/reference/{{sample}}/Aligned.out.bam",
    output:
        bam=f"{OUT}/align/{{sample}}.genome.sorted.bam",
        bai=f"{OUT}/align/{{sample}}.genome.sorted.bam.bai",
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
