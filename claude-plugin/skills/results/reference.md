# smORF-Bullpen results reference

Full documentation: <https://github.com/brendan-miller-salk/smORF-Bullpen/blob/main/docs/output.md>
and <https://github.com/brendan-miller-salk/smORF-Bullpen/blob/main/docs/methods.md>.

## Output files

| Path (under `outdir`) | Contents |
|---|---|
| `orfs.csv` | One row per ORF; the main table |
| `summary.csv` | Counts per set and type, fitted thresholds, QC vs GENCODE smORFs v2 |
| `README.md` | Run description and column definitions |
| `counts.tsv` | `feature_id`, `feature_type` (`canonical` or `unreviewed_smorf`), `length_nt`, one count column per sample, `rpkm_<sample>` |
| `reviewed_smorfs.{gtf,nuc,pep}` | Comprehensive smORFs with `blastoff_class == sp_identical` |
| `unreviewed_smorfs.{gtf,nuc,pep}` | All other Comprehensive smORFs |
| `scoring/orf_scores_annotated.tsv` | Full table with every annotator column (Ribo-seq mode) |
| `orfs_annotated.tsv` | Full table (provided mode) |
| `scoring/thresholds.tsv`, `scoring/cds_calibration.tsv` | Fitted thresholds and the CDS they came from |
| `scoring/pass_all_qc.tsv` | QC checks against GENCODE smORFs v2 ranges |
| `qc/strand_report.tsv`, `psites/psite_report.tsv`, `preprocess/adapter_report.tsv` | Per-library QC gates |
| `counting/flank_fit.tsv`, `counting/flank_metagene.png` | Fitted CDS flanks used for counting |
| `tracks/<ucsc_genome>/*.bw`, `tracks/frame_summary.tsv` | Frame-resolved P-site bigWigs |

## orf_set (Ribo-seq mode)

| Value | Meaning |
|---|---|
| `primary` | Called by PRICE or RiboCode, passes PIF, uniformity and drop-off thresholds |
| `comprehensive` | Called, not excluded, fails at least one score |
| `excluded` | Removed by a filter; reason after `excluded:` in `failed` |

"Comprehensive smORFs" in summaries means `primary` + `comprehensive`. Provided mode has no
`orf_set`.

## Key orfs.csv columns

| Column | Meaning |
|---|---|
| `orf_id`, `chrom`, `start`, `end`, `strand` | BED coordinates, 0-based half-open |
| `isoform_group` | Stop-codon position; start/splice variants of one ORF share it |
| `is_smorf` | At most `smorf_max_aa` (default 150) amino acids |
| `callers` | `PRICE`, `RiboCode` or `PRICE+RiboCode` |
| `ppm` | P-sites per million (column sums to 1e6) |
| `pif`, `uniformity`, `dropoff` | Translation signature scores; thresholds in `summary.csv` |
| `unshared_nt`, `unshared_inframe_psites` | ORF sequence not read in frame by any annotated CDS, and in-frame P-sites there |
| `*_unique` | Same measures from uniquely mapped reads (MAPQ 255) |
| `unique_reads`, `unique_rpkm` | Pooled counts outside flanked canonical CDS (unreviewed smORFs only) |
| `above_cds_median` | Primary smORFs with ribosome density above the dataset's median CDS (BH FDR < `cds_median_fdr`) |
| `failed` | Failed scores, `dropoff n/a`, and `excluded: <reason>` |
| `smorf_type`, `general_smorf_type` | Annotator type, refined by UniProt match; architecture group |
| `blastoff_class`, `blastoff_acc`, `blastoff_gene` | UniProt match |
| `shortstop_class`, `shortstop_probability` | `sam_secreted`, `sam_intracellular`, `prisms`, or `not_scored_*` |
| `diamond_n_species`, `diamond_species`, `diamond_n_untestable` | Conservation summary |
| `ucsc` | Genome browser link |

## blastoff_class

Prefix `sp_` (Swiss-Prot) or `tr_` (TrEMBL), then:

| Class | Meaning |
|---|---|
| `identical` | Protein is exactly the UniProt sequence (string match) |
| `near_identical` | >= 95% identity end to end, not exact |
| `insertion` / `deletion` | Same ends; extra / skipped internal segment |
| `N_truncated`, `C_truncated`, `internal` | Lies within the UniProt protein |
| `N_extended`, `C_extended`, `NC_extended` | Extends it at one or both ends |
| `partial` | Shares >= 10 aa in frame, then diverges |
| `homolog` | Similar, not in frame with a known protein |

`novel` has no prefix: no match.

## smorf_type and general_smorf_type

| general_smorf_type | smorf_type values |
|---|---|
| Upstream | uORF, uoORF, uaORF, uaoORF, udORF |
| Internal | iORF, riORF |
| Isoform | N-Iso, C-Iso, NC-Iso, Ins-Iso, Del-Iso, Trunc-Iso (from Swiss-Prot matches) |
| Downstream | dORF, doORF, daORF, daoORF |
| Non-coding | lncRNA, intergenic, psORF, aORF, eORF, UA, not_annotated_* |
| Reviewed | annotated_CDS (exact Swiss-Prot) |
| TrEMBL | TrEMBL (exact TrEMBL) |

Every ORF typed `UA` usually means the GTF was GENCODE rather than Ensembl.

## Filters

Each sets `orf_set` to `excluded` with the reason in `failed`.

| Config option (template value) | Excludes |
|---|---|
| `unshared_evidence: drop`, `unshared_min_inframe_psites` (5) | ORFs with own sequence but too few in-frame P-sites in it |
| `unique_evidence_types` ([psORF]), `unique_min_inframe_psites` (5) | Those types without enough uniquely mapped in-frame P-sites |
| `exclude_types` ([Trunc-Iso]) | Listed `smorf_type`s |
| `require_unshared_sequence` (true) | ORFs with no own sequence, unless exact Swiss-Prot (or exact TrEMBL with `trembl_identical_exempt: true`) |
| `exclude_chroms` ([chrM]) | ORFs on listed chromosomes |
| `min_unique_rpkm` (5) | Unreviewed smORFs at or below this pooled RPKM outside flanked CDS |

ORFs excluded at the final RPKM step keep their `counts.tsv` row; earlier exclusions are not
counted.
