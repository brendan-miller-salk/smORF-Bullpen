---
name: setup
description: Set up and launch a smORF-Bullpen run. Use when the user wants to run smORF-Bullpen on their Ribo-seq FASTQs or annotate an existing ORF GTF, needs a run directory, samplesheet or config, wants help choosing modules, or wants to dry-run or submit the pipeline.
argument-hint: "[folder of FASTQs or ORF GTF]"
---

# Set up a smORF-Bullpen run

Guide the user from data to a launched run. Work step by step, confirm choices before writing
files, and never submit a job without showing the exact command and getting approval.

Option details, module choices and reference checks are in [config-guide.md](config-guide.md).

## 1. Find the pipeline and the environment

- Locate the smORF-Bullpen repository: `$BULLPEN_REPO` if set, otherwise ask. If the user has
  none, offer `git clone https://github.com/brendan-miller-salk/smORF-Bullpen.git` into a
  location they choose (not inside the run directory).
- Check the head-job environment: `snakemake --version` (8 or newer) and `apptainer --version`
  (or `singularity`). If Snakemake is missing, offer
  `mamba env create -n bullpen-snakemake -f $BULLPEN_REPO/workflow/envs/snakemake.yaml`.
  On clusters with modules, `module spider apptainer` finds the container runtime.
- The head job inherits the shell environment when `sbatch` is called, so the environment must
  be active, and Apptainer loaded, in the shell that submits. A site profile can do this
  instead (see step 7).

## 2. Choose the mode

- **Ribo-seq mode**: the user has single-end Ribo-seq FASTQs (gzipped).
- **Provided mode**: the user has ORFs as a GTF and wants them annotated. Set `orf_gtf`; no
  samplesheet is needed.

If `$ARGUMENTS` points to a folder of FASTQs or a GTF, infer the mode and confirm it.

## 3. References

Ask where references live, or offer to build them with
`bash $BULLPEN_REPO/refs/fetch_references.sh /path/to/refs` (about 3 GB of downloads; mouse
GRCm39, Ensembl 110). Before using user-supplied references, run the checks in
[config-guide.md](config-guide.md#reference-checks). The most common silent failure is a
GENCODE GTF, which makes every ORF type `UA`.

If the user wants the `diamond` conservation module, the species files come from
`bash $BULLPEN_REPO/refs/fetch_conservation.sh /path/to/refs`, which prints the config block to
paste.

## 4. Create the run directory

Ask where the run should live (outside the repository, on storage with room for large outputs)
and what to call it. Then:

```bash
mkdir -p <run>/config <run>/logs
cp $BULLPEN_REPO/config/config.template.yaml <run>/config/bullpen.yaml
```

## 5. Samplesheet (Ribo-seq mode)

From the FASTQ folder, list `*.fastq.gz` / `*.fq.gz` and propose sample names (file name
without extension and lane/read suffixes). Show the table and let the user edit names before
writing `config/samples.tsv`:

```
sample<TAB>fastq<TAB>type<TAB>condition
brain_rep1<TAB>/absolute/path/brain_rep1.fastq.gz<TAB>elongation<TAB>brain
brain_harr1<TAB>/absolute/path/brain_harr1.fastq.gz<TAB>tis_harr<TAB>brain
```

Rules: tab-separated, header `sample`, `fastq` (then optional `type`, `condition`), unique
names, absolute paths, one single-end library per row. If you see R1/R2 pairs, ask: Ribo-seq is
single-end, and R2 files usually should not be included. Re-sequenced runs of one library need
merging first.

Ask whether any library is initiation-enriched (harringtonine or lactimidomycin; names with
`harr`, `ltm`, `lact` or `tis` are a hint). If none, the `type` and `condition` columns can be
left out. If some are, set `type` to `tis_harr` or `tis_ltm` for them (`elongation` for the
rest), set `condition` so each TIS library shares it with at least one elongation library, and
add `ribotish` to `callers` in step 6. See
[config-guide.md](config-guide.md#tis-libraries-and-ribo-tish).

## 6. Configure

Edit `config/bullpen.yaml` and show the diff of what you changed. Always set `genome`, `gtf`,
`outdir`, `index_dir` and, if `blastoff` is on, `blastoff_fasta`. Then choose modules with the
user using the questions in [config-guide.md](config-guide.md#choosing-modules), and set
`tools:` accordingly. Check the dependencies listed there (for example `featurecounts` needs
`blastoff`; assembly needs an RNA-seq samplesheet; `ribotish` needs TIS libraries).

Leave scoring and filter options at their template values unless the user asks; explain any
change you make.

## 7. Cluster settings

- Account: `sacctmgr show assoc user=$USER format=account%30,partition%20 -P`. Partition and
  limits: `sinfo -s`. Ask the user which to use; the head job needs 2 days on that partition.
- If the user will run several analyses, offer to create a site profile: copy
  `$BULLPEN_REPO/profiles/example/` to `$BULLPEN_REPO/profiles/<site>/` and fill in every value
  marked `EDIT` (account, partition, module loads, conda path, temp and cache directories). Site
  profiles are not tracked by git, so suggest the user keeps a backup.
- Without a scheduler, use the local profile (see [config-guide.md](config-guide.md#profiles)).

## 8. Dry run, then launch

From the run directory, show the dry-run command and submit only after approval:

```bash
export BULLPEN_REPO=/path/to/smORF-Bullpen
sbatch -A <account> -p <partition> $BULLPEN_REPO/profiles/slurm/run_head.sbatch config/bullpen.yaml -n
```

When it finishes (`squeue -u $USER` no longer lists `bullpen_head`), read
`logs/bullpen_head_<jobid>.out`. Check for `WorkflowError`, and that the job table matches the
chosen modules (provided mode shows no `trim`, `star_align`, `price`, `ribocode` or `score`;
with TIS libraries, one `ribotish_predict` per TIS library).
Fix any errors with the user before continuing.

Then show the real launch (the same command without `-n`) and submit only after approval. Tell
the user the run takes hours to days, that progress can be checked with the `status` skill, and
that re-running the same command resumes a stopped run.

## Guardrails

- Never write inside the pipeline repository except when creating a site profile the user asked
  for.
- Never submit jobs, delete files or overwrite an existing config without explicit approval.
- Do not run the pipeline or heavy tools on a login node; everything runs through `sbatch`.
