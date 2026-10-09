# Visualizer: ORF-relative P-site frame bigWigs for a genome browser.
# Pools every library's P-site bedGraph, splits each nucleotide by its frame from
# the start codon of the ORF it lies in (canonical CDS first, then ORF candidates),
# and scales to P-sites per million in-ORF P-sites (see scripts/frame_tracks.py).
# Hosting and the UCSC hub file are a separate step (the project's /hf-track-upload):
# one stacked multiWig row per strand, blue frame 0, red 1, gold 2.

TRK = f"{OUT}/tracks"
TRACKS_GENOME = config.get("ucsc_genome", "mm39")
TRACKS_LABEL = config.get("tracks_label") or os.path.basename(os.getcwd())
UCSC_BIGZIPS = f"https://hgdownload.soe.ucsc.edu/goldenPath/{TRACKS_GENOME}/bigZips"
TRACK_FILES = expand(f"{TRK}/{TRACKS_GENOME}/{TRACKS_LABEL}_{{strand}}_f{{frame}}.bw",
                     strand=["plus", "minus"], frame=[0, 1, 2])

wildcard_constraints:
    trk_strand="plus|minus",
    trk_frame="[012]",


# Downloads run on the head node: compute nodes on many clusters have no internet.
localrules:
    ucsc_chrom_files,


rule ucsc_chrom_files:
    """UCSC chromosome sizes and the alias table that renames contigs (GL456210.1 -> chr1_GL456210v1_random)."""
    output:
        alias=f"{IDX}/ucsc/{TRACKS_GENOME}.chromAlias.txt",
        sizes=f"{IDX}/ucsc/{TRACKS_GENOME}.chrom.sizes",
    retries: 2
    shell:
        r"""
        curl -fsSL -o {output.alias} {UCSC_BIGZIPS}/{TRACKS_GENOME}.chromAlias.txt
        curl -fsSL -o {output.sizes} {UCSC_BIGZIPS}/{TRACKS_GENOME}.chrom.sizes
        """


rule frame_tracks:
    """Pooled P-sites split by ORF-relative frame, as UCSC-named bedGraphs."""
    input:
        fw=expand(f"{OUT}/psites/{{sample}}_fw.bedgraph", sample=SAMPLES),
        rc=expand(f"{OUT}/psites/{{sample}}_rc.bedgraph", sample=SAMPLES),
        canonical=f"{IDX}/reference.canonical.gtf",
        orfs=f"{SC}/orf_candidates.bed12",
        fai=f"{IDX}/genome.fa.fai",
        alias=f"{IDX}/ucsc/{TRACKS_GENOME}.chromAlias.txt",
    output:
        bedgraphs=temp(expand(f"{TRK}/bedgraph/{TRACKS_LABEL}_{{strand}}_f{{frame}}.bedgraph",
                              strand=["plus", "minus"], frame=[0, 1, 2])),
        summary=f"{TRK}/frame_summary.tsv",
    params:
        label=TRACKS_LABEL,
    resources:
        mem_mb=scaled(16000),   # 0.9 GB on heart (160M P-sites); arrays are per chromosome
        runtime=60,
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/tracks/frame_tracks.log",
    shell:
        r"""
        python {SCRIPTS}/frame_tracks.py --fw {input.fw} --rc {input.rc} \
            --canonical-gtf {input.canonical} --orfs {input.orfs} --fai {input.fai} \
            --alias {input.alias} --label {params.label} --bedgraph-dir {TRK}/bedgraph \
            --summary {output.summary} > {log} 2>&1
        """


rule frame_bigwig:
    input:
        bedgraph=f"{TRK}/bedgraph/{TRACKS_LABEL}_{{trk_strand}}_f{{trk_frame}}.bedgraph",
        sizes=f"{IDX}/ucsc/{TRACKS_GENOME}.chrom.sizes",
    output:
        f"{TRK}/{TRACKS_GENOME}/{TRACKS_LABEL}_{{trk_strand}}_f{{trk_frame}}.bw",
    resources:
        mem_mb=scaled(4000),
        runtime=30,
    container:
        IMG.get("bedgraphtobigwig", "docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:482--hdc0a859_0")
    shell:
        "bedGraphToBigWig {input.bedgraph} {input.sizes} {output}"
