# smORF-Bullpen: Methods

## Contents

- [Workflow](#workflow)
- [Preprocessing and quality control](#preprocessing-and-quality-control)
- [ORF calling](#orf-calling)
- [ORF catalogue](#orf-catalogue)
- [Scoring](#scoring)
- [Evidence filters](#evidence-filters)
- [Tools](#tools)
  - [Annotator](#annotator)
  - [blastoff](#blastoff)
  - [ShortStop](#shortstop)
  - [featureCounts](#featurecounts)
  - [DIAMOND](#diamond)
  - [Frame tracks](#frame-tracks)
  - [StringTie](#stringtie)
- [Conventions](#conventions)

## Workflow

### Ribo-seq mode

```
[optional: bulk RNA-seq -> fastp -> STAR -> StringTie per sample -> merge -> TPM filter -> novel transcripts]

Trim Galore -> SortMeRNA
  -> STAR: reference index (P-sites, PRICE, Ribo-TISH) and hybrid index (RiboCode)
  -> strandedness check per library
  -> P-site offsets and fit check per library
  -> ORF calling: PRICE (pooled), RiboCode (per library), Ribo-TISH (per TIS library)
  -> cohort ORF catalogue
  -> translation signature scoring
  -> annotators and conservation, in parallel
  -> join: merge annotations, refine smorf_type, apply exclusions
       |-> pass_all QC against Chothani et al.
       |-> counting: unreviewed smORFs outside flanked canonical CDS -> counts.tsv
  -> final: count gate, set assignment -> orfs.csv, summary.csv, README.md
  -> reviewed and unreviewed smORF GTF and FASTA files

Frame tracks branch off the P-sites and do not depend on ORF calls.
```

The workflow never loops back. Counting branches off the join because it counts only
unreviewed, non-excluded smORFs, which are known only after the join. The final step waits for
counting because counts decide set membership.

### Provided mode

```
ORF GTF -> catalogue (stop codons folded in, chromosome names matched to the genome,
           identical structures merged)
        -> annotators -> orfs_annotated.tsv -> orfs.csv, summary.csv, README.md
```

## Preprocessing and quality control

Reads are trimmed with Trim Galore and depleted of rRNA with SortMeRNA against its eight default
databases (Rfam 5S and 5.8S, and SILVA archaeal, bacterial and eukaryotic small and large
subunit rRNA). Before rRNA removal, each sample is checked for its adapter: it must be found in
at least `adapter_min_fraction` of the first million reads.

Reads are aligned with STAR twice: to the reference genome, whose alignments give the P-sites
used for scoring and PRICE, and to a hybrid index that includes novel transcripts, whose
transcriptome alignments RiboCode uses. `ribocode_annotation` selects the backbone RiboCode
sees: one transcript per gene (`canonical`, the longest isoform chosen by AGAT) or every Ensembl
isoform (`full`). PRICE and the catalogue always use every isoform plus novel transcripts.

Two checks gate each library before it is used:

- **Strandedness.** featureCounts counts the library's unique reads on annotated CDS on each
  strand. At least `strand_min_sense` (0.9) must be sense, from at least `strand_min_assigned`
  reads. A library's P-sites and RiboCode run wait for its own check, and PRICE waits for all
  of them.
- **P-site offsets.** plastid fits an offset per read length (`psite --aggregate`). A read
  length is kept only if plastid fitted it (rather than falling back to the default), it makes
  up at least `psite_min_share` of start-codon reads, its offset lies in `psite_offset_range`,
  and at least `psite_min_frame0` of its P-sites are in frame. A library fails if it has too few
  start-codon reads, if its dominant read length is not kept, or if kept lengths cover less than
  `psite_min_kept_share` of its reads.
- **TIS guard.** `start_peak_ratio` is P-sites per codon on codons -1..+1 of annotated start
  codons over codons 11-19, pooled over read lengths with an offset in range. An elongation
  library at or above `tis_peak_max` (default 12) stops the run as a likely unlabelled TIS
  library; PRICE and RiboCode wait for this check.

P-site tracks from several libraries are pooled by summing overlapping bedGraph intervals
exactly.

## ORF calling

- **PRICE** runs once on all libraries pooled. PRICE reports every candidate it tested, not only
  its calls, so the table is filtered to a Benjamini-Hochberg FDR below `price_max_fdr` (0.1,
  PRICE's own default). `price_max_p` adds an optional cap on the raw p-value. Gene IDs come
  from the transcript PRICE names.
- **RiboCode** runs per library at its default cut-off (combined p ≤ 0.05).
- **Ribo-TISH** runs `predict` once per TIS library (`type: tis_harr` or `tis_ltm`), with that
  library as TIS data (`-t`, plus `--harr` for harringtonine) and the elongation libraries of
  its `condition` as ordinary Ribo-seq (`-b`). Each call combines a TIS test (start-codon peak)
  and a frame test (elongation body); calls are Ribo-TISH's own output at Fisher q ≤
  `ribotish_max_q` (0.05, its `--fsqth`), and `pval_ribotish` is the Fisher p-value. The Fisher
  q passes any in-frame start inside a translated CDS on frame evidence alone, so the catalogue
  also requires an initiation peak at the start codon, TIS q ≤ `ribotish_max_tis_q` (0.05; on
  chr22 of one LTM library, 35,926 calls became 1,680, 86% AUG), and drops calls that run off
  the transcript without a stop codon. Candidates
  are every Ensembl isoform plus novel transcripts, as for PRICE; the TIS background and the
  known-CDS labels (`-a`) use the canonical GTF, one transcript per gene. `ribotish_alt: true` tests AUG
  and near-cognate starts (one base from AUG, `--alt`); `false` tests AUG only. P-site offsets
  come from `ribotish quality` on each library's reference genome alignments against the
  canonical GTF (`-t` for TIS libraries); the step fails if no read length passes Ribo-TISH's
  quality threshold. Each TIS library counts as one sample in `n_samples`. Ribo-TISH's TisType
  maps to the catalogue class: Annotated, Truncated and Extended to `canonical_cds`, 5'UTR to
  `uORF`, 3'UTR to `dORF`, Internal to `intORF` and Novel to `novel_u`, ignoring `:Known` and
  `:CDSFrameOverlap` suffixes. Its `--blocks` give the genomic structure directly (0-based
  half-open, stop included). Ribo-TISH 0.2.7 and 0.2.8 crash with more than one process on
  Python 3.8 or later (`group argument must be None for now`), so `predict` runs through a small
  inline shim that drops the pool context passed as `group` and is otherwise unchanged.

TIS libraries are used only by Ribo-TISH. Their reads pile up at start codons, which would
distort PRICE's model, RiboCode's periodicity and the PIF, uniformity and drop-off scores, so
they are excluded from P-site offset fitting and its check, scoring, counting and frame tracks.

Chromosome names from all callers are harmonised to the genome before merging, and each
caller's stop-codon convention is measured against the genome and reported in
`catalogue/catalogue_report.tsv`.

## ORF catalogue

Calls from all callers and libraries are merged into one cohort catalogue.

- `catalogue_merge: exact` (template default) keeps every distinct ORF structure as its own
  row, retaining all isoforms as in Chothani et al. Calls of the same structure from several
  callers, or from several transcripts, become one row credited to each caller.
- `catalogue_merge: overlap` folds ORFs with at least 80% reciprocal overlap and identical
  splice junctions where they overlap into one representative, chosen by class and then length.
  Calls that differ only in their start codon are merged, so the called start can be lost.
- `orf_collapse: true` additionally folds smORFs whose proteins MMseqs2 clusters together (at
  least 90% identity and 80% coverage), within one locus only. Copies elsewhere, such as
  processed pseudogenes, keep their own rows. Surviving rows list their members in
  `collapsed_orf_ids`. The template default is `false`.

## Scoring

Every ORF is scored on three translation signatures, as defined for GENCODE smORFs v2
(Chothani et al. 2026; Menon et al. 2026):

| Score | Measures |
|---|---|
| PIF | Proportion of P-sites in frame |
| Uniformity | Share of codons whose P-sites are mostly in frame |
| Drop-off | In-frame P-site density before the stop codon relative to after it |

Thresholds are fitted to each dataset from complete annotated CDS, never hardcoded. A CDS is
used for calibration if it has no `cds_start_NF` or `cds_end_NF` tag, a 3-nt stop codon, a first
phase of 0, a length that is a multiple of 3, and at least `calibration_min_psites` (100)
P-sites. The stop codon is appended as it is for catalogue ORFs, so the drop-off window sits in
the same place. With the template defaults, each threshold is set where 95% of annotated CDS
score above it (`threshold_tail: lower`, `threshold_percentile: 0.95`, `threshold_method:
empirical`).

ORFs that pass every score in `require` are Primary; the rest are Comprehensive. Codons without
P-sites are ignored when computing uniformity (`uniformity_zero_codons: ignore`): counting them
as failures pushes the threshold below zero on a single library, so that every ORF passes.

P-sites per million (`ppm`) is reported but is not a criterion by default, since gating on
expression penalises lowly expressed smORFs. Setting `min_ppm` to a number enables it.

`scoring/pass_all_qc.tsv` compares the run's thresholds, Primary fraction and length and overlap
biases with the ranges reported by Chothani et al.

## Evidence filters

An ORF can pass the translation signatures on P-sites it shares with a highly translated
annotated CDS. These filters require each ORF to have evidence of its own.

**Unshared sequence.** Every ORF gets `unshared_nt`, `unshared_psites` and
`unshared_inframe_psites`: the nucleotides that no annotated CDS (from any Ensembl transcript,
stop codons included) reads in the same frame, and the pooled P-sites there. With
`unshared_frame: aware` (default), sequence that a CDS covers only in another frame counts as
the ORF's own, since upstream, internal and downstream overlapping ORFs are alternative-frame
overlaps. With `blind`, any CDS overlap counts as shared.

| Option (template value) | Effect |
|---|---|
| `unshared_evidence: drop`, `unshared_min_inframe_psites: 5` | An ORF with unshared sequence but fewer than 5 in-frame P-sites there is excluded |
| `require_unshared_sequence: true` | An ORF with no unshared sequence is excluded, unless its protein is exactly a Swiss-Prot sequence |
| `trembl_identical_exempt: true` | Exact TrEMBL matches are also exempt from the previous rule and from counting; they must still pass the scores |
| `unique_evidence_types: [psORF]`, `unique_min_inframe_psites: 5` | These types need 5 in-frame P-sites from uniquely mapped reads in their unshared sequence; otherwise a pseudogene almost identical to its parent gene picks up the parent's reads |
| `exclude_types: [Trunc-Iso]` | Truncated isoforms lie entirely inside a known protein, so all their P-sites belong to the host CDS |
| `exclude_chroms: [chrM]` | The callers use the standard genetic code, so mitochondrial ORF structures and proteins are wrong |
| `min_unique_rpkm: 5` | An unreviewed smORF needs a pooled RPKM above 5 outside flanked canonical CDS (see [featureCounts](#featurecounts)) |

Uniquely mapped reads are those with STAR MAPQ 255. STAR keeps up to 20 loci per read and
plastid counts each, so `psites/unique/` holds a second set of P-site tracks, with the same
offsets, built from uniquely mapped reads only.

Primary smORFs are also tested for ribosome density above the dataset's median annotated CDS,
using uniquely mapped reads over the codons between start and stop (one-sided Poisson test,
BH FDR below `cds_median_fdr`). For ORFs that overlap an annotated CDS and are not reviewed, only
codons that no CDS touches are used. The result is the `above_cds_median` flag and does not
change `orf_set`.

## Tools

### Annotator

The Annotator assigns each ORF a `smorf_type` from its position relative to Ensembl transcripts,
for example uORF, uoORF, iORF, dORF, lncRNA, intergenic or psORF. Scaffold ORFs (`GL*`, `JH*`,
`MU*`) are typed `not_annotated_scaffold` and kept. `general_smorf_type` groups the types by
architecture:

| Group | Types |
|---|---|
| Upstream | uORF, uoORF, uaORF, uaoORF, udORF |
| Internal | iORF, riORF |
| Isoform | N-Iso, C-Iso, NC-Iso, Ins-Iso, Del-Iso, Trunc-Iso |
| Downstream | dORF, doORF, daORF, daoORF |
| Non-coding | lncRNA, intergenic, psORF, aORF, eORF, UA, not_annotated_* |
| Reviewed | annotated_CDS |
| TrEMBL | TrEMBL |

### blastoff

blastoff classifies each ORF's protein against the mouse UniProt reference proteome, including
isoforms, using BLASTP (`blastp-short` for proteins of at most 30 amino acids).

An ORF is `identical` when its protein, read from the literal first codon, is exactly a UniProt
sequence. This is a string comparison, not a BLAST result, and Swiss-Prot canonical sequences
take precedence over isoforms, which take precedence over TrEMBL. Otherwise, alignments are
classified by how the two proteins' ends compare:

| Class | Meaning |
|---|---|
| `near_identical` | At least `blastoff_min_pident` (95%) identity end to end, but not exact |
| `insertion` | Same ends; the ORF is longer, with an extra internal exon or segment |
| `deletion` | Same ends; the ORF is shorter, skipping an exon |
| `N_truncated`, `C_truncated`, `internal` | The ORF lies within the UniProt protein |
| `N_extended`, `C_extended`, `NC_extended` | The ORF extends the UniProt protein at one or both ends |
| `partial` | Shares an in-frame stretch of at least 10 amino acids, then diverges |
| `homolog` | Similar, but not in frame with a known protein |
| `novel` | No match |

Each class is prefixed `sp_` (Swiss-Prot) or `tr_` (TrEMBL). Every in-frame class except
`identical` needs an alignment of at least `blastoff_min_shared_aa` (10) amino acids, so a 5-9
amino acid ORF is either identical, a homolog or novel. `blastoff_offset_aa` gives the amino
acids gained or lost.

For long queries, up to 10 HSPs per protein are kept, and HSPs in the same order on both
sequences are chained into one hit, so an exon inserted or skipped in the middle of an alignment
still gives a single match. Identity is computed over ungapped columns.

With the Annotator also enabled, UniProt matches refine `smorf_type`:

| Match | `smorf_type` |
|---|---|
| Exact Swiss-Prot | `annotated_CDS` |
| Exact TrEMBL | `TrEMBL` |
| Swiss-Prot truncated or internal | `Trunc-Iso` (shorter, no residues of its own) |
| Swiss-Prot extended or partial | `N-Iso`, `C-Iso` or `NC-Iso`, by the end where the ORF has residues of its own |
| Swiss-Prot insertion or deletion | `Ins-Iso` or `Del-Iso` |
| Any other match | The Annotator's type is kept |

Isoform relabelling applies only to internal and overlapping types (`iORF`, `uoORF`, `doORF`,
`uaoORF`, `daoORF`, `uaORF`, `daORF`). The original call is kept in `smorf_type_annotator`.

blastoff can also be run on its own: `workflow/scripts/blastoff/blastoff.py query|classify`.

### ShortStop

[ShortStop](https://github.com/brendan-miller-salk/ShortStop) classifies smORF proteins as
secreted (`sam_secreted`), intracellular (`sam_intracellular`) or `prisms`, with a probability.
ORFs it cannot score are labelled `not_scored_non_AUG`, `not_scored_length` or
`not_scored_sequence`. Every submitted ORF must come back with the expected protein
(`protein_check.tsv`), or the step fails. ShortStop is installed at a pinned commit (release
v1.1.0), and its bundled model is trained on human data.

ShortStop v1.1 also scores each ORF with 100 consensus models, each voting SAM when
P(intracellular) + P(secreted) > 0.5. `shortstop_sam_consensus_percentage` is the share voting
SAM, and `shortstop_sam_consensus` the highest level reached: `SAM-Consensus-100` (all models),
`-90` (>= 90%), `-50` (> 50%) or `none`. The levels are nested; `-90` or `-100` is the most
robust set. `shortstop_class` stays the standard model's call. The model bundle is downloaded
and checksummed when the ShortStop env is built.

### featureCounts

featureCounts counts unique reads (no `-M`) on canonical CDS and on unreviewed smORFs:
Comprehensive smORFs of at most `smorf_max_aa` amino acids whose protein is not exactly a
Swiss-Prot sequence. Only the part of each smORF that no canonical CDS reaches is counted, after
each canonical CDS is widened by a 5' and a 3' flank. The subtraction ignores frame: reads in
another frame of a CDS are attributed to the CDS.

The flanks are fitted per dataset from read spill-over past canonical start and stop codons. A
flank ends where coverage falls to within `counting_flank_excess_frac` (5%) of the way from the
adjacent UTR level (the median over `counting_flank_baseline`, 40-100 nt) to the CDS level.
`counting_flank_nt` replaces both fitted flanks with a fixed length. Counting uses one
featureCounts pass with `-O -s 1`, so only forward-stranded libraries are supported.

With `min_unique_rpkm` set, a counted smORF whose pooled RPKM is not above it moves to
`excluded`. Pooled RPKM is the read count summed over libraries, per kb of counted length, per
million reads assigned in `counts.tsv`.

### DIAMOND

DIAMOND tests whether each ORF's protein is encoded in other species. Because DIAMOND has no
tblastn mode, each species' Ensembl cDNA and ncRNA are used as blastx queries against a database
of ORF proteins (`--ultra-sensitive --min-orf 1 --masking 0 -k 0`, minimum score
`diamond_search_bits`). Reverse-frame hits are discarded, since transcripts are sense.

E-values are recomputed from the bitscore as tblastn would report them, so they do not depend on
the run's ORF set. They are computed against two databases:

- The species' whole transcriptome (`diamond_<species>_evalue`).
- The transcripts of orthologs of the ORF's host gene (`diamond_<species>_ortholog_evalue`). This
  much smaller search space lets short ORFs reach significance.

Each species gets a call:

| Call | Meaning |
|---|---|
| `conserved` | Best forward hit at E ≤ `diamond_evalue` (10<sup>-3</sup>), covering at least `diamond_min_cov` (50%) of the protein |
| `untestable` | Even a perfect match could not reach significance in either database |
| `not_conserved` | Testable, with no qualifying hit |

On mouse muscle compared with human, 0.8% of ORFs of at most 15 amino acids were testable
against the transcriptome alone, and 39% with the ortholog database; for 16-20 amino acids, 61%
and 89%; every ORF over 25 amino acids was testable. Intergenic ORFs have no host gene and so no
ortholog database.

### Frame tracks

Pooled P-sites are assigned a frame from the start codon of the ORF they fall in: canonical CDS
first (honouring the first CDS phase), then catalogue ORFs where there is no CDS. A nucleotide
that two ORFs of the same tier place in different frames is dropped, as are P-sites outside
ORFs; the shares are reported in `tracks/frame_summary.tsv`. Tracks are scaled per million
in-ORF P-sites, because mitochondrial rRNA survives depletion. Contigs are renamed with UCSC's
chromAlias table.

### StringTie

Novel transcripts are assembled from bulk RNA-seq of the same cells:

1. **fastp**, per lane: adapters detected automatically or given (`rnaseq_adapter_r1`,
   `rnaseq_adapter_r2`), poly-G trimming, 3' quality trimming (window 4, Q20), and removal of
   pairs with a mate shorter than 36 nt. Each lane is read directly from its `.gz` file, so a
   corrupt archive fails the run.
2. **STAR**, 2-pass, against the reference index with ENCODE long-RNA settings. Strandedness is
   inferred from STAR gene counts (`rnaseq_strandedness: auto`), and all samples must agree.
3. **StringTie** 3, reference-guided, per sample; chromosomes in `exclude_chroms` are skipped.
4. **StringTie `--merge`** across samples (`assembly_merge_args`, default `-F 0 -T 1`).
5. **StringTie `-e`** to estimate TPM per sample on the merged assembly.

A merged transcript absent from the reference is kept if its TPM is at least `assembly_min_tpm`
(1) in at least `assembly_min_samples` (2) samples. `novel_class_codes` then decides which
gffcompare classes reach the callers: for example `u` (intergenic), `i` (intronic) and `x`
(antisense). Unstranded single-exon models have no strand and are dropped.

## Conventions

- Coordinates are BED: 0-based, half-open. BED12 is used for spliced ORFs.
- References are Ensembl 110, GRCm39, with `chr`-prefixed chromosome names.

**Other species.** Only `refs/fetch_references.sh` and `ucsc_genome` (used for links) are
specific to mouse. The Annotator needs an Ensembl-format GTF with chromosomes named `chr`
followed by digits, `X`, `Y` or `M` (names such as `chr2L` or `chrMT` will not parse), and
ShortStop's bundled model is trained on human data.
