# smORF-Bullpen: Output

## Contents

- [Main results](#main-results)
- [ORF sets](#orf-sets)
- [orfs.csv columns](#orfscsv-columns)
- [smORF GTF and FASTA files](#smorf-gtf-and-fasta-files)
- [Read counts](#read-counts)
- [Frame tracks](#frame-tracks)
- [Quality control reports](#quality-control-reports)
- [Intermediate files](#intermediate-files)

## Main results

All paths are relative to `outdir`. The files at the top level are the ones to open first.

| File | Contents |
|---|---|
| `orfs.csv` | One row per ORF: coordinates, protein, scores, annotations, counts and set membership |
| `summary.csv` | smORF counts per set and per `general_smorf_type` |
| `README.md` | A description of the run and of every output column |
| `counts.tsv` | Read counts per sample for canonical genes and unreviewed smORFs (`featurecounts`) |
| `reviewed_smorfs.{gtf,nuc,pep}` | smORFs whose protein is exactly a Swiss-Prot sequence (`blastoff`) |
| `unreviewed_smorfs.{gtf,nuc,pep}` | All other smORFs (`blastoff`) |

## ORF sets

In Ribo-seq mode, the `orf_set` column places every ORF in exactly one set:

| `orf_set` | Meaning |
|---|---|
| `primary` | Passes every translation signature threshold in `require` (PIF, uniformity, drop-off) |
| `comprehensive` | Called by a caller and not excluded, but below at least one threshold |
| `excluded` | Removed by one of the filters below; the reason is recorded in `failed` |

Excluded ORFs stay in `orfs.csv`, so every call can be traced. Each filter runs at the first
step that has the information it needs:

| Step | Filter (config option) | Outcome |
|---|---|---|
| Scoring | PIF, uniformity or drop-off below the fitted threshold (`require`) | Comprehensive, not Primary |
| Scoring | Has sequence of its own but fewer than `unshared_min_inframe_psites` in-frame P-sites there (`unshared_evidence: drop`) | Excluded |
| Join | A `unique_evidence_types` ORF (e.g. psORF) with fewer than `unique_min_inframe_psites` P-sites from uniquely mapped reads in its own sequence | Excluded |
| Join | `smorf_type` listed in `exclude_types` (e.g. Trunc-Iso) | Excluded |
| Join | No sequence of its own and not a reviewed protein (`require_unshared_sequence`) | Excluded |
| Join | On a chromosome listed in `exclude_chroms` (e.g. chrM) | Excluded |
| Final | Unreviewed smORF with pooled RPKM at or below `min_unique_rpkm` outside flanked canonical CDS | Excluded |

ORFs excluded at scoring or the join are not counted and have no row in `counts.tsv`. ORFs
excluded at the final step keep their row.

Primary smORFs also get an `above_cds_median` flag. It is `True` when the smORF's in-frame
ribosome density is significantly above the dataset's median annotated CDS (one-sided Poisson
test, Benjamini-Hochberg FDR below `cds_median_fdr`). The flag does not change `orf_set`.

In provided mode there are no sets: ORFs were not called from Ribo-seq, so `orfs.csv` has no
`orf_set`, `callers` or score columns.

## orfs.csv columns

| Column | Contents | Present |
|---|---|---|
| `orf_id` | ORF identifier | Always |
| `source_id` | The ORF's ID in the provided GTF (`orf_id` attribute, else `transcript_id`) | Provided mode |
| `ucsc` | UCSC Genome Browser link to the ORF, optionally within a saved session (`ucsc_session`) | Always |
| `chrom`, `start`, `end`, `strand` | Coordinates, 0-based half-open (BED) | Always |
| `isoform_group` | Position of the ORF's stop codon; start and splice variants of one ORF share it | Always |
| `host_gene_name`, `gene_id` | Host gene name and Ensembl gene ID | Always |
| `orf_set` | `primary`, `comprehensive` or `excluded` (see [ORF sets](#orf-sets)) | Ribo-seq mode |
| `smorf_type` | ORF type relative to the host gene, refined by the UniProt match | `annotator` |
| `general_smorf_type` | `smorf_type` grouped by architecture: Upstream, Internal, Isoform, Downstream, Non-coding, Reviewed, TrEMBL | `annotator` |
| `is_smorf` | `True` if at most `smorf_max_aa` amino acids (150 by default); a flag, not a filter | Always |
| `aa_length` | Protein length in amino acids | Always |
| `callers` | Which callers found the ORF: `PRICE`, `RiboCode`, `Ribo-TISH`, joined with `+` | Ribo-seq mode |
| `price_p` | PRICE p-value | Ribo-seq mode |
| `ppm` | P-sites per million, scaled so the column sums to 10<sup>6</sup> | Ribo-seq mode |
| `pif`, `uniformity`, `dropoff` | Translation signature scores (see [Scoring](methods.md#scoring)) | Ribo-seq mode |
| `unshared_nt`, `unshared_psites`, `unshared_inframe_psites` | Length of the ORF's own sequence (not read in frame by any annotated CDS), and the P-sites in it | Ribo-seq mode |
| `psites_unique`, `ppm_unique`, `unshared_inframe_psites_unique` | The same P-site measures from uniquely mapped reads only | Ribo-seq mode |
| `unique_reads`, `unique_rpkm` | Pooled reads and RPKM outside flanked canonical CDS (unreviewed smORFs only) | `featurecounts` |
| `above_cds_median`, `cds_median_density_ratio`, `cds_median_fdr` | Ribosome density compared with the median annotated CDS (Primary smORFs only) | Ribo-seq mode |
| `failed` | Scores below threshold, `dropoff n/a` if too few reads to measure it, and `excluded: <reason>` | Ribo-seq mode |
| `blastoff_class`, `blastoff_acc`, `blastoff_gene` | UniProt match class, accession and gene (see [blastoff](methods.md#blastoff)) | `blastoff` |
| `shortstop_class`, `shortstop_probability` | ShortStop class and probability | `shortstop` |
| `shortstop_sam_consensus`, `shortstop_sam_consensus_percentage` | highest SAM-Consensus level (`SAM-Consensus-100`/`-90`/`-50`/`none`) and % of 100 models voting SAM | `shortstop` |
| `diamond_n_species`, `diamond_species`, `diamond_n_untestable` | Species with a conserved ortholog, and species where the test had no power | `diamond` |
| `aa_seq` | Protein sequence | Always |

The `README.md` written to each output directory defines every column for that run. The full
table, with every column from every annotator (such as `smorf_type_annotator`, `host_gene_id`,
the remaining `blastoff_*` fields and per-species `diamond_<species>_*` calls), is
`scoring/orf_scores_annotated.tsv`, or `orfs_annotated.tsv` in provided mode.

The catalogue (`catalogue/cohort.catalogue.tsv`) has `called_by_<caller>` (1 or 0) and
`pval_<caller>` for `price` (PRICE p-value) and `ribocode` (combined p-value), plus
`called_by_ribotish` and `pval_ribotish` (Fisher p-value) when Ribo-TISH ran. `n_samples` counts
the RiboCode libraries and Ribo-TISH TIS libraries that called the ORF; PRICE runs pooled.

## smORF GTF and FASTA files

With `blastoff` on, the Comprehensive smORFs (Primary and comprehensive, at most `smorf_max_aa`
amino acids; in provided mode, all smORFs) are split by their match to UniProt:

| Set | `blastoff_class` |
|---|---|
| `reviewed_smorfs` | `sp_identical`: the protein is exactly a Swiss-Prot sequence |
| `unreviewed_smorfs` | Everything else: novel ORFs, all TrEMBL matches, and inexact Swiss-Prot matches such as isoforms |

Each set is written as three files:

| File | Contents |
|---|---|
| `.gtf` | Per ORF: `transcript`, `exon`, `CDS` (with phase), `start_codon` and `stop_codon` rows, in transcription order. `gene_id` and `transcript_id` are both the `orf_id`. |
| `.nuc` | CDS nucleotides from the genome, stop codon included |
| `.pep` | Protein sequence (`aa_seq`): literal first codon, stop removed |

The run fails if a genome translation disagrees with `aa_seq` or lacks a stop codon.
`scoring/primary_set.bed12` holds the Primary ORFs in BED12 format.

## Read counts

`counts.tsv` holds unique-read counts per sample for canonical genes and for unreviewed smORFs,
with `feature_type` distinguishing the two. smORFs are counted only on the part that no
canonical CDS reaches, after each CDS is widened by a 5' and 3' flank fitted to the dataset.

| Column | Contents |
|---|---|
| `length_nt` | Counted length; a smORF that lies entirely inside flanked CDS has length 0 and zero counts |
| `<sample>` | Unique reads assigned per library |
| `rpkm_<sample>` | Reads per kb of `length_nt` per million reads assigned in the table; blank when `length_nt` is 0 |

The flank fit is reported in `counting/flank_fit.tsv` and plotted in
`counting/flank_metagene.png`. See [Counting](methods.md#featurecounts) for the method.

## Frame tracks

With `frame_tracks` listed under `visualizers:`, pooled P-sites are written as bigWig files
split by strand and reading frame:

```
tracks/<ucsc_genome>/<label>_{plus,minus}_f{0,1,2}.bw
tracks/frame_summary.tsv
```

`<label>` is `tracks_label`, or the run directory name if unset. Tracks are scaled per million
P-sites that fall inside ORFs. No track hub file is written. In the UCSC Genome Browser, show each
strand as a single stacked multiWig track with frame 0 in blue (`0,0,255`), frame 1 in red
(`255,0,0`) and frame 2 in gold (`230,170,0`).

## Quality control reports

By default each check stops the run when a library fails. Most can be set to `warn` to report
and continue.

| Report | Checks | Config |
|---|---|---|
| `preprocess/adapter_report.tsv` | The adapter is present in at least `adapter_min_fraction` of each sample's first million reads | `adapter`, `adapter_min_fraction` |
| `qc/strand_report.tsv` | At least `strand_min_sense` of each library's unique CDS reads are sense-stranded | `strand_check` |
| `psites/psite_report.tsv` | Each library has read lengths with a fitted, in-frame P-site offset | `psite_check` |
| `scoring/pass_all_qc.tsv` | Thresholds, Primary fraction and length and overlap biases compared with Chothani et al. (warns by default) | `pass_all_qc` |
| `callers/ribotish/quality/<sample>_qual.pdf` | Ribo-TISH found at least one read length with a usable P-site offset (fails otherwise) | `ribotish` |
| `catalogue/catalogue_report.tsv` | Stop-codon convention of each caller, measured against the genome | |
| `assembly/strand_report.tsv`, `assembly/adapter_report.tsv`, `assembly/assembly_report.tsv` | RNA-seq strandedness, adapters and per-step transcript counts (`stringtie`) | |

## Intermediate files

| Directory | Contents |
|---|---|
| `preprocess/` | Trimmed and rRNA-depleted reads, trimming reports |
| `align/` | STAR genome and transcriptome alignments, against the reference and the hybrid index |
| `qc/` | Strandedness counts per library |
| `psites/` | P-site offsets and bedGraph tracks; `psites/unique/` holds the same from uniquely mapped reads |
| `callers/` | Raw PRICE, RiboCode and Ribo-TISH output; `callers/ribotish/` holds `<tis_sample>.txt` (calls) and `<tis_sample>.tisBackground.txt`, and `callers/ribotish/quality/` the P-site offsets per library (`<sample>.para.py`, `<sample>_qual.{txt,pdf}`) |
| `catalogue/` | The merged cohort catalogue and MMseqs2 clusters |
| `scoring/` | CDS calibration, fitted thresholds, ORF scores and the annotated score table |
| `annotations/` | One table per annotator, keyed by `orf_id` |
| `counting/` | Flank fit and counting GTF |
| `assembly/` | StringTie assembly and reports (`stringtie`) |
| `index/` | STAR, GEDI and SortMeRNA indexes (or `index_dir`) |
| `logs/` | One log per job |
