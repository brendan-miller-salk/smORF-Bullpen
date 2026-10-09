---
name: results
description: Explain and query the results of a finished smORF-Bullpen run (orfs.csv, summary.csv, counts.tsv, smORF sets, QC reports). Use when the user asks what a run found, how many smORFs passed, why an ORF was excluded, or wants to filter ORFs by type, conservation, UniProt match or expression.
argument-hint: "[run directory or output directory]"
---

# Explore smORF-Bullpen results

Help the user understand and query a finished run. Answer from the files; never guess numbers.

## 1. Find the output directory

- If `$ARGUMENTS` is given, use it. A run directory contains `config/bullpen.yaml`; its `outdir`
  key (relative to the run directory) is the output directory. An output directory contains
  `orfs.csv`.
- Otherwise look in the current directory for `config/bullpen.yaml` or `orfs.csv`, and ask if
  neither is there.
- Read `config/bullpen.yaml` for the settings that shape the results: `orf_gtf` (set = provided
  mode), the `tools:` block, `smorf_max_aa`, and the filters listed in
  [reference.md](reference.md#filters).

## 2. Start with the run's own summary

Read `<outdir>/README.md` first. The pipeline writes it for each run: headline counts, the QC
comparison against GENCODE smORFs v2, and a definition of every column present. Then read
`summary.csv` (columns `section, item, primary, comprehensive, excluded, all, note`):

| `section` | Rows |
|---|---|
| `orfs` | all ORFs, smORFs, smORFs above CDS median |
| `type (all ORFs)`, `type (smORFs)` | counts per `smorf_type` |
| `general type (smORFs)` | Upstream, Internal, Isoform, Downstream, Non-coding, Reviewed, TrEMBL |
| `threshold` | fitted PIF, uniformity and drop-off cut-offs and how many CDS they came from |
| `qc vs Chothani / GENCODE v2 2025` | OK or WARN per check, observed vs expected |

Give the user the headline first (Primary and Comprehensive smORFs), then anything flagged WARN.
Provided-mode runs have no sets, scores or thresholds; report counts per type instead.

## 3. Answer questions with small, explicit queries

`orfs.csv` is one row per ORF and is usually a few MB, so load it with pandas and show the
code you ran. Useful filters:

```python
import pandas as pd
orfs = pd.read_csv(f"{outdir}/orfs.csv")
smorfs = orfs[orfs.is_smorf & (orfs.orf_set != "excluded")]        # Comprehensive smORFs
primary = smorfs[smorfs.orf_set == "primary"]
novel = smorfs[smorfs.blastoff_class == "novel"]
conserved_human = smorfs[smorfs.diamond_species.fillna("").str.contains("human")]
excluded = orfs[orfs.orf_set == "excluded"]                        # reason in `failed`
```

- **"Why was this ORF excluded?"** Read its `failed` value; the text after `excluded:` names the
  rule. Map it to the config option with [reference.md](reference.md#filters).
- **Conservation.** `diamond_n_species` counts species with a conserved hit. A species listed in
  `diamond_n_untestable` could not have produced a significant hit, so absence there is not
  evidence of non-conservation. Per-species calls are in `scoring/orf_scores_annotated.tsv`
  (`diamond_<species>_call`).
- **Expression.** `unique_rpkm` is pooled across libraries; per-sample counts and RPKM are in
  `counts.tsv` (`feature_type == "unreviewed_smorf"`, joined on `feature_id == orf_id`).
- **Viewing an ORF.** The `ucsc` column is a ready-made browser link. Frame tracks, if built, are
  in `tracks/<ucsc_genome>/`.
- **Columns not in `orfs.csv`** (`smorf_type_annotator`, `host_gene_id`, other `blastoff_*`
  fields, per-species DIAMOND results) are in `scoring/orf_scores_annotated.tsv`
  (`orfs_annotated.tsv` in provided mode), keyed by `orf_id`.

Column meanings and class definitions are in [reference.md](reference.md).

## 4. Report

- Lead with the answer and the count, then the filter used, so the user can reproduce it.
- Name files and columns exactly.
- If the user wants a table or list saved, write it inside the run directory (for example
  `<outdir>/../analysis/`), never inside the pipeline repository, and say where it went.
- Do not modify anything under the output directory.
