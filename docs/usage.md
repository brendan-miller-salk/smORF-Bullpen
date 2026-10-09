# smORF-Bullpen: Usage

## Contents

- [Installation](#installation)
- [References](#references)
- [Setting up a run](#setting-up-a-run)
  - [Ribo-seq mode](#ribo-seq-mode)
  - [Provided mode](#provided-mode)
- [Running the pipeline](#running-the-pipeline)
- [Choosing tools](#choosing-tools)
- [Novel transcripts](#novel-transcripts)
- [Profiles](#profiles)
- [Resources](#resources)
- [Adding a tool](#adding-a-tool)
- [Tests](#tests)

## Installation

The pipeline needs conda (or mamba) and Apptainer (or Singularity). Snakemake and its Slurm
executor run in a small conda environment; every other tool runs in a pinned container image or
a conda environment that Snakemake builds on the first run.

```bash
git clone https://github.com/brendan-miller-salk/smORF-Bullpen.git
mamba env create -n bullpen-snakemake -f smORF-Bullpen/workflow/envs/snakemake.yaml
mamba activate bullpen-snakemake
```

## References

```bash
bash smORF-Bullpen/refs/fetch_references.sh /path/to/refs
```

This downloads and builds about 3 GB of references: the GRCm39 primary assembly, the Ensembl 110
GTF with `chr`-prefixed chromosome names, and the mouse UniProt reference proteome (a pinned
release) for `blastoff`. The script exits with an error if the GTF lacks the features the
Annotator needs.

References must satisfy two conditions:

- **The GTF must be from Ensembl, not GENCODE.** The Annotator relies on `five_prime_utr` and
  `gene_biotype`; without them every ORF is typed `UA`.
- **Chromosomes must use UCSC-style names** (`chr1` to `chrY`, `chrM`).

For the `diamond` conservation tool, `refs/fetch_conservation.sh` downloads Ensembl cDNA and
ncRNA for human, rat, cow, chicken and zebrafish, together with Ensembl Compara orthologs of
mouse genes, and prints the matching config block. The `SOURCE`, `RELEASE` and `MART`
variables change the run species and Ensembl release.

## Setting up a run

Each analysis lives in its own run directory outside the repository. Outputs, logs and
Snakemake's state are written there, and all commands below are run from inside it.

### Ribo-seq mode

```bash
mkdir -p ~/runs/brain/config ~/runs/brain/logs && cd ~/runs/brain
cp $BULLPEN_REPO/config/config.template.yaml config/bullpen.yaml
cp $BULLPEN_REPO/config/samples.template.tsv config/samples.tsv
```

`config/samples.tsv` is tab-separated, with one gzipped single-end library per row and absolute
paths:

```tsv
sample	fastq
brain_rep1	/data/riboseq/brain_rep1.fastq.gz
brain_rep2	/data/riboseq/brain_rep2.fastq.gz
```

In `config/bullpen.yaml`, set at least `genome`, `gtf`, `outdir`, `blastoff_fasta` and the
`tools:` block. Every option is documented in the template.

> [!NOTE]
> Libraries must be sense-stranded. UMIs are not handled. If the library uses a custom 3'
> linker, give its sequence as `adapter`, because the run checks that the adapter is present in
> at least half of each sample's first million reads.

### Provided mode

Provided mode annotates a set of ORFs you already have. Set `orf_gtf` to turn it on:

```bash
mkdir -p ~/runs/my_smorfs/config ~/runs/my_smorfs/logs && cd ~/runs/my_smorfs
cp $BULLPEN_REPO/config/config.template.yaml config/bullpen.yaml
```

```yaml
orf_gtf: /path/to/my_smorfs.gtf
outdir: results_my_smorfs
genome: /path/to/refs/GRCm39.primary_assembly.genome.fa
gtf: /path/to/refs/ensembl110.GRCm39.chr.gtf   # the Annotator types ORFs against this
index_dir: results_my_smorfs/index
tools:
  callers: []          # required: provided mode has no callers
  annotators: [annotator, shortstop, blastoff]
  counters: []         # required: there is nothing to count
```

No samplesheet is needed. Trimming, alignment, ORF calling, scoring and counting are not run.

The GTF needs a `transcript` row and its `CDS` rows for each ORF, one ORF per transcript.
`stop_codon` rows are folded in if present, and a CDS that ends one codon short of a stop is
extended to include it. Chromosomes may be bare (`7`) or `chr`-prefixed. ORFs with identical
structures are merged into one row. The run stops and lists the offending ORFs if any have CDS
rows without a transcript row, a length that is not a multiple of 3, no stop codon, or a
chromosome that is not in the genome.

## Running the pipeline

### Dry run

```bash
export BULLPEN_REPO=/path/to/smORF-Bullpen
sbatch -A <account> -p <partition> $BULLPEN_REPO/profiles/slurm/run_head.sbatch config/bullpen.yaml -n
```

Once the job has finished, check the job table in `logs/bullpen_head_<jobid>.out`. Provided mode
should list no `trim`, `star_align`, `price`, `ribocode` or `score` jobs.

### Launch

Drop the `-n`:

```bash
sbatch -A <account> -p <partition> $BULLPEN_REPO/profiles/slurm/run_head.sbatch config/bullpen.yaml
```

This submits one small head job (2 CPUs, 8 GB, up to 2 days) that runs Snakemake, which in turn
submits each step as its own Slurm job. Follow progress with `squeue -u $USER` and
`tail -f logs/bullpen_head_<jobid>.out`. Per-step logs are written to `<outdir>/logs/`.

### Resuming

Re-running the same command resumes the run and skips finished steps. Snakemake does not detect
edits to scripts, so use `--forcerun <rule>` to redo a step.

### Targets

Targets are passed after the config file, in the same way as `-n`:

| Target | Runs |
|---|---|
| `all` (default) | The full pipeline |
| `stage1` | Up to the ORF catalogue and P-site tracks (Ribo-seq mode only) |

## Choosing tools

The `tools:` block of the config decides which modules run. Each tool is one file in
[`workflow/rules/tools/`](../workflow/rules/tools/), and its columns appear in `orfs.csv` only if
it ran.

| Tool | Kind | Adds |
|---|---|---|
| `stringtie` | assembly | Novel transcripts from bulk RNA-seq |
| `price` | caller | ORF calls, from all libraries pooled |
| `ribocode` | caller | ORF calls, per library |
| `annotator` | annotator | `smorf_type`, `general_smorf_type`, `host_gene_name`, `host_gene_id` |
| `shortstop` | annotator | `shortstop_class`, `shortstop_probability` |
| `blastoff` | annotator | Match to UniProt (`blastoff_class`, `blastoff_acc`, `blastoff_gene`, ...) |
| `diamond` | conservation | `diamond_*` columns: orthologous hits in other species |
| `featurecounts` | counter | `counts.tsv`; `unique_reads` and `unique_rpkm` in `orfs.csv` |
| `frame_tracks` | visualizer | Frame-resolved P-site bigWig tracks |

How each tool works is described in the [methods documentation](methods.md#tools).

### Common configurations

| Goal | `annotators:` | `counters:` |
|---|---|---|
| Everything (template default) | `[annotator, shortstop, blastoff]` | `[featurecounts]` |
| smORF type only | `[annotator]` | `[]` |
| ShortStop only | `[shortstop]` | `[]` |
| Type and ShortStop, without UniProt | `[annotator, shortstop]` | `[]` |
| Coordinates and protein sequence only | `[]` | `[]` |

### Dependencies between tools

- `featurecounts` needs `blastoff`, which decides which smORFs are counted. Requesting it without
  `blastoff` stops the run before any job is submitted.
- With both `annotator` and `blastoff`, UniProt matches refine the Annotator's `smorf_type`. The
  original call is kept in `smorf_type_annotator`.
- `stringtie` needs `rnaseq_samplesheet` and replaces `novel_gtf`; setting both is an error.
- `visualizers:`, `conservation:` and `assembly:` are off when absent from the config.
- Provided mode never calls, counts, assembles or draws tracks.

## Novel transcripts

ORFs on transcripts that Ensembl lacks can be found by supplying novel transcript models in one
of two ways:

- **Supply an assembly:** set `novel_gtf` to a StringTie assembly annotated with gffcompare
  class codes. Set `novel_gtf_stranded: true` if it was assembled from stranded RNA-seq.
- **Build one:** set `tools.assembly: [stringtie]` and point `rnaseq_samplesheet` at bulk
  RNA-seq of the same cells. The samplesheet has the columns `sample`, `fastq_1` and `fastq_2`
  (leave `fastq_2` empty for single-end), with one row per lane.

`novel_class_codes` picks which gffcompare classes PRICE and the catalogue see.
`ribocode_novel_class_codes` picks the subset that RiboCode and its STAR index see. Classes that
repeat a reference intron chain (`k`, `m`, `n`) make transcriptome reads multi-map, and RiboCode
then loses periodicity, so the template excludes them on the RiboCode side.

The antisense class `x` needs stranded RNA-seq. With `tools.assembly`, every library must
measure as `rf` or `fr` in `assembly/strand_report.tsv`, or the run stops after alignment.

> [!WARNING]
> With `tools.assembly`, the hybrid GTFs and STAR index in `index_dir` depend on that run's
> RNA-seq. Do not share `index_dir` with a run that assembles different data.

## Profiles

| Profile | Use | Launch, from the run directory |
|---|---|---|
| `profiles/slurm` | Any Slurm cluster | `sbatch -A <account> -p <partition> $BULLPEN_REPO/profiles/slurm/run_head.sbatch config/bullpen.yaml` |
| `profiles/local` | One machine, no scheduler | `snakemake -s $BULLPEN_REPO/workflow/Snakefile --configfile config/bullpen.yaml --profile $BULLPEN_REPO/profiles/local --cores 16 --resources mem_mb=120000` |
| `profiles/<site>` | Your cluster, configured once | `sbatch $BULLPEN_REPO/profiles/<site>/run_head.sbatch config/bullpen.yaml` |

To set up a site profile, copy `profiles/example/` to `profiles/<site>/` and replace every value
marked `EDIT`: account, partition, bind paths, environment and image locations, `module load`
lines and `SITE`. A site profile can fix the account, partition and repository path, so they no
longer need to be given at launch.

Site profiles are excluded by `.gitignore` so that account names and paths stay out of the
repository. Keep a separate copy of yours.

## Resources

Per-rule CPU, memory and time requests are set in [`workflow/rules/`](../workflow/rules/) and
can be overridden with `set-resources` in the profile. The two heaviest steps are:

- **STAR genome index:** 96 GB, rising to 192 GB on retry.
- **PRICE:** memory scales with the number of libraries.

Tool images and conda environments are built once. Set `conda-prefix` and `apptainer-prefix` in
the profile to share them between runs, and point `index_dir` at a shared directory to build the
STAR, GEDI and SortMeRNA indexes once.

The SortMeRNA databases and UCSC chromosome aliases are downloaded on the head node, so compute
nodes do not need internet access for those steps.

## Adding a tool

1. Write `workflow/rules/tools/<name>.smk`.
2. Add `<name>` to `KNOWN_TOOLS` in [`workflow/rules/common.smk`](../workflow/rules/common.smk).
3. Pin its software in `workflow/envs/` or the config's `containers:` block.

A caller's output is read by `workflow/scripts/build_catalogue.py`. An annotator or conservation
tool must write `<outdir>/annotations/<name>.tsv`, with `orf_id` as the first column and one row
per ORF. The `join_annotations` rule then merges it automatically.

## Tests

Each test script prints PASS or FAIL and exits with status 1 on failure. Run them with a Python
that has numpy, pandas and pyfaidx, such as the `workflow/envs/score.yaml` environment:

```bash
for t in tests/test_*.py; do python "$t" > /dev/null && echo "ok   $t" || echo "FAIL $t"; done
```
