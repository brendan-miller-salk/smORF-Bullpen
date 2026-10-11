# smORF-Bullpen failure reference

## Config errors

Raised while Snakemake reads the config, before any job runs. The head log shows `WorkflowError`.

| Message contains | Cause and fix |
|---|---|
| `Config is missing` | Required keys absent: `outdir`, `genome`, `gtf`, `index_dir`, `containers`, plus `samplesheet` and `sortmerna_dbs` in Ribo-seq mode. Start from `config/config.template.yaml`. |
| `Unknown <kind>` | A `tools:` entry names a tool that does not exist; the message lists valid ones. |
| `orf_gtf is set ... must be [] or absent` | Provided mode with callers, counters, visualizers or assembly listed. Set those to `[]`. |
| `tools.callers is empty` | Ribo-seq mode needs at least one of `price`, `ribocode`, `ribotish`. |
| `Set novel_gtf or tools.assembly, not both` | Choose one source of novel transcripts. |
| `needs a strand-specific novel_gtf` | `x` in `novel_class_codes` with a supplied `novel_gtf`; set `novel_gtf_stranded: true` only if it came from stranded RNA-seq, otherwise drop `x`. |
| `rnaseq_strandedness is unstranded` | `x` cannot be used with unstranded RNA-seq. |
| `tools.counters needs blastoff` | Add `blastoff` to `annotators` or set `counters: []`. |
| `need tab-separated columns 'sample' and 'fastq'` | Samplesheet header or separator wrong (spaces instead of tabs is common). |
| `duplicate sample` | Sample names must be unique. |
| `must be an absolute path` | FASTQ paths in the samplesheet must be absolute. |
| `type of <sample> is '<t>'` | `type` must be `elongation`, `tis_harr` or `tis_ltm` (blank means `elongation`). |
| `no elongation libraries` | Every library is TIS. Add ordinary Ribo-seq libraries (`type` `elongation` or blank). |
| `are only used by Ribo-TISH` | TIS libraries without `ribotish` in `tools.callers`. Add it, or drop those rows. |
| `has ribotish but the samplesheet has no TIS libraries` | Mark the TIS libraries `tis_harr` / `tis_ltm`, or remove `ribotish`. |
| `look like TIS libraries by name ... but have no type` | Ask the user whether they are harringtonine / lactimidomycin libraries: set `type` to `tis_harr` / `tis_ltm` (plus a shared `condition`), or `elongation` if not. |
| `have no elongation library with the same condition` | Fix `condition` so each TIS library shares it with an elongation library (blank means `all`). |
| `containers has no 'ribotish'` | Add `ribotish: docker://quay.io/biocontainers/ribotish:0.2.8--pyhdfd78af_0` under `containers:`. |
| `ribocode_annotation` | Must be `canonical` or `full`. |
| blastoff / diamond / stringtie messages | Missing `blastoff_fasta`, `diamond_species` (names lowercase `[a-z0-9_]`; `diamond_orthologs` keys must be a subset), `rnaseq_samplesheet`, or a container entry. |

## Quality gates

Deliberate stops that protect the results from bad input.

| Rule | Message | Meaning | Fix |
|---|---|---|---|
| `adapter_check` | `only X% of reads contain ...` | The 3' adapter is not in at least `adapter_min_fraction` (0.5) of the first million reads; usually a custom linker Trim Galore cannot detect | Set `adapter:` to the actual sequence. No warn mode. |
| `adapter_report` | `trimmed with N different adapters` | Auto-detection picked different adapters per sample | Set `adapter:` explicitly. |
| `strand_check` / `strand_report` | `not sense-stranded or too few CDS reads` | A library is below `strand_min_sense` (0.9) sense reads on CDS, or has fewer than `strand_min_assigned` (100,000) CDS reads. See `qc/strand_report.tsv`. | Check the library prep; antisense libraries are not supported. `strand_check: warn` continues, but P-sites and calls will be wrong for a non-sense library. |
| `psite_report` | `P-site check failed for: ...` | A library has too few start-codon reads (`psite_min_start_reads`), its dominant read length has no usable offset, or kept lengths cover too little of it. See `psites/psite_report.tsv`. | `psite_check: warn` scores the library on its kept lengths. Consider dropping the library. |
| `psite_report` | `look like TIS (harringtonine / lactimidomycin) libraries` | An elongation library's start-codon peak is >= `tis_peak_max` x its CDS body (`start_peak_ratio` in `psites/psite_report.tsv`): almost always an unlabelled TIS library. | Fails even with `warn`. Ask the user; set its `type` to `tis_harr` / `tis_ltm`. Raise `tis_peak_max` only if they confirm it is elongation. |
| `psite_report` | `no read length has a usable P-site offset` | No periodicity at all in a library | Fails even with `warn`; remove the library from the samplesheet. |
| `normalize` | `callers disagree on the stop-codon convention` / `not consistent` | Caller outputs mix stop-inclusive and stop-exclusive coordinates | Set `extend_stop: always` or `never` after checking `catalogue/catalogue_report.tsv`. |
| `normalize` | empty catalogue | No ORFs were called | Check `price.log`, `ribocode/*.log` and `ribotish/*.log`, and the P-site report. |
| `ribocode` | log contains `Error` | RiboCode failed for a library, often no periodic read length | Read `<outdir>/logs/ribocode/<sample>.log`. |
| `ribotish_quality` | `ribotish quality: no read length passed for <sample>` | The library has no usable P-site offset for Ribo-TISH | Check `callers/ribotish/quality/<sample>_qual.pdf`; drop the library or fix its `type`. Log: `logs/ribotish/quality_<sample>.log`. |
| `ribotish_predict` | `AssertionError: group argument must be None for now` | Ribo-TISH 0.2.7/0.2.8 multiprocessing on Python >= 3.8; the rule's shim should prevent it | Report as a bug with `logs/ribotish/predict_<tis>.log`. |
| `shortstop` | `protein check failed for N of M submitted ORFs` | ShortStop returned a different protein than expected. See `annotations/shortstop/protein_check.tsv`. | Report as a bug with that file. |
| `smorf_sets` | `ORFs not in` / `genome translation disagrees with aa_seq` | Internal consistency check between the table and the genome | Report as a bug with `logs/smorf_sets.log`. |
| `pass_all_qc` | WARN lines | Results fall outside GENCODE smORFs v2 ranges (`scoring/pass_all_qc.tsv`) | Warns only; stops only with `pass_all_qc: fail`. |
| `provided_catalogue` | CDS without transcript row, off-chromosome, not a multiple of 3, no stop, seqnames not in genome | Problems in a provided-mode GTF | Fix the GTF; the message lists the ORFs. |
| `rnaseq_strand` | `cannot infer one strandedness` | RNA-seq samples disagree on strandedness | Set `rnaseq_strandedness` explicitly, or remove the odd sample. |

## Resources

Memory requests double on the single retry (`retries: 1` in the profiles).

| Rule | Request |
|---|---|
| `star_index` | 32 threads, 96 GB (192 GB on retry), 4 h |
| `star_align` | 32 threads, 56 GB, 2 h |
| `sortmerna` | 32 threads, 64 GB, 10 h |
| `price` | 16 threads, max(16 GB, 6 GB x libraries), 10 h |
| `ribotish_quality` | 4 threads, 16 GB, 4 h |
| `ribotish_predict` | 16 threads, 64 GB, 24 h |
| `plastid_tracks`, `plastid_tracks_unique` | 48 GB, 2 h |
| `canonical_gtf` (AGAT) | 32 GB, 4 h |
| `shortstop_predict` | 24 GB, 8 h |
| `diamond_blastx` | 16 threads, 32 GB, 6 h |
| `rnaseq_star` | 32 threads, 56 GB, 6 h |

`OUT_OF_MEMORY` or `TIMEOUT` in `sacct` after the retry: raise the request in the profile, e.g.

```yaml
set-resources:
  price:
    mem_mb: 192000
    runtime: 1440
```

## Other

| Symptom | Cause and fix |
|---|---|
| Every ORF typed `UA` | GTF is GENCODE, not Ensembl. Rebuild references with `refs/fetch_references.sh`. |
| `sortmerna_db` or `ucsc_chrom_files` fails | These download on the head node (GitHub, UCSC). Check internet access from the node running the head job. |
| `LockException` | Previous head job died. If `squeue` shows no head job for this run, rerun once with `--unlock`, then launch normally. |
| `snakemake not found` / no apptainer | `run_head.sbatch` exits if either is missing from the head job's environment; check the `module load` lines or conda activation in the site profile. |
| `BULLPEN_REPO` unset | Export it before `sbatch`, or use a site profile that sets it. |

Head job log: `logs/bullpen_head_<jobid>.out` (stdout and stderr). Per-step logs:
`<outdir>/logs/`; reference steps: `<index_dir>/logs/`; Slurm logs:
`.snakemake/slurm_logs/rule_<name>/`; Snakemake logs: `.snakemake/log/`.
