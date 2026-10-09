# smORF-Bullpen configuration guide

Every option is documented inline in `config/config.template.yaml`. Full usage guide:
<https://github.com/brendan-miller-salk/smORF-Bullpen/blob/main/docs/usage.md>.

## Reference checks

Run these before using references the user supplies:

```bash
# Ensembl GTF, not GENCODE: Ensembl has five_prime_utr features and gene_biotype.
zcat -f genes.gtf | awk '$3=="five_prime_utr"' | head -1
zcat -f genes.gtf | grep -m1 -o 'gene_biotype "[^"]*"'
# UCSC-style chromosome names in both files, and the same names.
grep '^>' genome.fa | head -3
zcat -f genes.gtf | grep -v '^#' | cut -f1 | sort -u | head
```

Requirements: Ensembl GTF (GENCODE silently types every ORF `UA`); chromosomes named `chr1` to
`chrY` and `chrM`; GTF and FASTA names match. `refs/fetch_references.sh` builds a set that
passes. Other species: only `fetch_references.sh` and `ucsc_genome` are mouse-specific, but
chromosome names must be `chr` plus digits, `X`, `Y` or `M`, and ShortStop's model is
human-trained.

## Required keys

| Key | Value |
|---|---|
| `samplesheet` | `config/samples.tsv` (Ribo-seq mode) |
| `orf_gtf` | Path to the ORF GTF (provided mode only; `null` otherwise) |
| `outdir` | Output directory, relative to the run directory |
| `genome`, `gtf` | Reference FASTA and Ensembl GTF |
| `index_dir` | Where STAR, GEDI and SortMeRNA indexes are built; share one between runs on the same references (but not with assembly runs) |
| `blastoff_fasta` | UniProt FASTA from `fetch_references.sh` (needed by `blastoff`) |
| `ucsc_genome` | UCSC assembly for links and tracks (`mm39` for mouse GRCm39) |

## Choosing modules

Ask these questions, then set `tools:`:

| Question | If yes |
|---|---|
| Do you have bulk RNA-seq of the same cells and want ORFs on unannotated transcripts? | `assembly: [stringtie]`, `novel_gtf: null`, and an RNA-seq samplesheet (`sample`, `fastq_1`, `fastq_2`; one row per lane; `fastq_2` empty for single-end) at `rnaseq_samplesheet` |
| Already have a StringTie/gffcompare assembly? | `novel_gtf: <path>` (and `novel_gtf_stranded: true` if stranded) instead of assembly |
| Want both ORF callers? (recommended) | `callers: [price, ribocode]` |
| Need smORF types (uORF, iORF, ...)? | `annotator` |
| Need to know which ORFs are known proteins, isoforms or novel? | `blastoff` |
| Interested in secreted vs intracellular microproteins? | `shortstop` |
| Need conservation in other species? | `conservation: [diamond]` and the `diamond_species` / `diamond_orthologs` block from `fetch_conservation.sh` |
| Need per-sample counts or the expression filter? | `counters: [featurecounts]` (requires `blastoff`) |
| Want genome browser tracks by frame? | `visualizers: [frame_tracks]` |

Dependencies: `featurecounts` needs `blastoff`; provided mode allows only `annotators` and
`conservation` (`callers`, `counters`, `visualizers`, `assembly` must be `[]` or absent);
`novel_gtf` and `assembly` are mutually exclusive.

## Options users commonly change

| Option | Template | When to change |
|---|---|---|
| `adapter` | `auto` | Custom 3' linker (e.g. `CTGTAGGCACCATCAAT`): give the sequence, or the adapter check fails |
| `smorf_max_aa` | 150 | Different smORF length cut-off (flag only) |
| `min_unique_rpkm` | 5 | Expression filter for unreviewed smORFs; remove to disable |
| `exclude_types` | `[Trunc-Iso]` | Keep or drop whole smORF types |
| `strand_check`, `psite_check` | `fail` | `warn` to continue past a failing library (explain the risk) |
| `ucsc_session` | `null` | A saved UCSC session (`user/name`) for the links in `orfs.csv` |
| `tracks_label` | `null` | Track name prefix; defaults to the run directory name |

## Profiles

| Profile | Launch from the run directory |
|---|---|
| `profiles/slurm` | `sbatch -A <account> -p <partition> $BULLPEN_REPO/profiles/slurm/run_head.sbatch config/bullpen.yaml` |
| `profiles/<site>` | `sbatch $BULLPEN_REPO/profiles/<site>/run_head.sbatch config/bullpen.yaml` |
| `profiles/local` | `snakemake -s $BULLPEN_REPO/workflow/Snakefile --configfile config/bullpen.yaml --profile $BULLPEN_REPO/profiles/local --cores 16 --resources mem_mb=120000` (needs about 96 GB for the STAR index) |

Append `-n` for a dry run, `stage1` to stop after the catalogue and P-site tracks, or
`--forcerun <rule>` to redo a step. Tool images and environments are built on the first run,
which needs internet access; `conda-prefix` and `apptainer-prefix` in the profile share them
between runs.
