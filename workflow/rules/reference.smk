# References and indexes. Everything here depends only on genome + GTFs, so it
# lands in IDX, which several runs may share.


rule genome_fai:
    input:
        config["genome"],
    output:
        f"{IDX}/genome.fa.fai",
    container:
        IMG["star"]
    shell:
        "samtools faidx {input} --fai-idx {output}"


# The GTF set (all built from the filtered reference):
#   REF_GTF      filtered Ensembl       -> STAR "reference" index: plastid, PRICE BAMs
#   HYBRID       backbone + novel       -> STAR "hybrid" index, RiboCode
#   FULL_HYBRID  reference + novel      -> GEDI/PRICE index, catalogue projection
# backbone = reference.canonical.gtf (default) or REF_GTF (ribocode_annotation: full).
REF_GTF = f"{IDX}/reference.filtered.gtf"
HYBRID = f"{IDX}/hybrid_reference.gtf"
FULL_HYBRID = f"{IDX}/full_hybrid_reference.gtf"
RIBOCODE_ANNOTATION = config.get("ribocode_annotation", "canonical")
if RIBOCODE_ANNOTATION not in ("canonical", "full"):
    raise WorkflowError("ribocode_annotation must be 'canonical' or 'full'.")
BACKBONE = f"{IDX}/reference.canonical.gtf" if RIBOCODE_ANNOTATION == "canonical" else REF_GTF

if NOVEL_GTF:

    rule gffcompare_novel:
        input:
            ref=REF_GTF,
            novel=NOVEL_GTF,
        output:
            f"{IDX}/gffcompare/novel.annotated.gtf",
        params:
            prefix=f"{IDX}/gffcompare/novel",
        container:
            IMG["gffcompare"]
        log:
            f"{IDX}/logs/gffcompare.log",
        shell:
            "gffcompare -r {input.ref} -o {params.prefix} {input.novel} > {log} 2>&1"


rule hybrid_gtf:
    input:
        ref=REF_GTF,
        backbone=BACKBONE,
        annotated=([f"{IDX}/gffcompare/novel.annotated.gtf"] if NOVEL_GTF else []),
    output:
        hybrid=HYBRID,
        full=FULL_HYBRID,
        report=f"{IDX}/hybrid_gffcompare.tsv",
    resources:
        mem_mb=scaled(16000),
        runtime=60,
    params:
        codes=",".join(config.get("novel_class_codes", ["u", "i", "j"])),
        ribo_codes=",".join(config.get("ribocode_novel_class_codes")
                            or config.get("novel_class_codes", ["u", "i", "j"])),
        annotated=lambda wc, input: f"--annotated {input.annotated}" if input.annotated else "",
    conda:
        "../envs/score.yaml"
    log:
        f"{IDX}/logs/hybrid_gtf.log",
    shell:
        r"""
        python {SCRIPTS}/hybrid_gtf.py --reference {input.ref} --backbone {input.backbone}             {params.annotated} --class-codes {params.codes}             --ribocode-class-codes {params.ribo_codes}             --out-hybrid {output.hybrid} --out-full-hybrid {output.full}             --report {output.report} 2> {log}
        cat {log} >&2
        """


rule gtf_filter:
    """Keep lines on seqnames in the FASTA that carry a transcript_id."""
    input:
        gtf=config["gtf"],
        fai=f"{IDX}/genome.fa.fai",
    output:
        f"{IDX}/reference.filtered.gtf",
    conda:
        "../envs/score.yaml"
    shell:
        "python {SCRIPTS}/gtf_filter.py --gtf {input.gtf} --fai {input.fai} --out {output}"


rule canonical_gtf:
    """One transcript per gene (reference.canonical.gtf): AGAT's
    longest isoform (longest CDS for coding genes, else longest exons), then
    gffread back to GTF. Used by plastid's metagene."""
    input:
        f"{IDX}/reference.filtered.gtf",
    output:
        gff=f"{IDX}/canonical/reference.canonical_agat.longest.gff",
    params:
        dir=f"{IDX}/canonical",
        # Absolute paths from Python: the container's busybox realpath has no -m.
        gtf=lambda wc, input: os.path.abspath(input[0]),
        out=lambda wc, output: os.path.abspath(output.gff),
        log=os.path.abspath(f"{IDX}/logs/agat.log"),
    threads: 4
    resources:
        mem_mb=scaled(32000),
        runtime=240,
    container:
        IMG["agat"]
    log:
        f"{IDX}/logs/agat.log",
    shell:
        r"""
        mkdir -p {params.dir} && cd {params.dir}
        agat_sp_keep_longest_isoform.pl --gff {params.gtf:q} --cpu {threads} \
            --out {params.out:q} > {params.log:q} 2>&1
        """


rule canonical_gtf_convert:
    input:
        f"{IDX}/canonical/reference.canonical_agat.longest.gff",
    output:
        f"{IDX}/reference.canonical.gtf",
    container:
        IMG["gffread"]
    shell:
        "gffread {input} --keep-exon-attrs -F -T -o {output}"


rule star_index:
    """Two indexes: `reference` (filtered Ensembl) and `hybrid`."""
    input:
        fasta=config["genome"],
        gtf=lambda wc: REF_GTF if wc.ref == "reference" else HYBRID,
    output:
        directory(f"{IDX}/star_{{ref}}"),
    wildcard_constraints:
        ref="reference|hybrid",
    threads: 32
    resources:
        mem_mb=scaled(96000),
        runtime=240,   # 13 min alone; >2 h when sharing I/O with a cohort's SortMeRNA
    container:
        IMG["star"]
    log:
        f"{IDX}/logs/star_index_{{ref}}.log",
    shell:
        r"""
        mkdir -p {output}
        STAR --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --genomeSAindexNbases 14 \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} \
            --limitGenomeGenerateRAM $(( {resources.mem_mb} * 900000 )) \
            --outFileNamePrefix {output}/ > {log} 2>&1
        """


SORTMERNA_URL = "https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases"


# Downloads run on the head node: compute nodes on many clusters have no internet.
localrules:
    sortmerna_db,


rule sortmerna_db:
    output:
        f"{IDX}/sortmerna/db/{{db}}.fasta",
    retries: 2
    shell:
        "curl -fsSL -o {output} {SORTMERNA_URL}/{wildcards.db}.fasta"


rule sortmerna_index:
    input:
        expand(f"{IDX}/sortmerna/db/{{db}}.fasta", db=config.get("sortmerna_dbs", [])),
    output:
        directory(f"{IDX}/sortmerna/idx"),
    params:
        refs=lambda wc, input: " ".join(f"--ref {f}" for f in input),
        workdir=f"{IDX}/sortmerna/work_index",
    threads: 16
    resources:
        mem_mb=scaled(32000),
        runtime=120,
    container:
        IMG["sortmerna"]
    log:
        f"{IDX}/logs/sortmerna_index.log",
    shell:
        r"""
        rm -rf {params.workdir}; mkdir -p {params.workdir} {output}
        sortmerna {params.refs} --idx-dir {output} --workdir {params.workdir} \
            --threads {threads} --index 1 > {log} 2>&1
        rm -rf {params.workdir}
        """
