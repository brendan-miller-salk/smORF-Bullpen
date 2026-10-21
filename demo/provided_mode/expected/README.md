# Provided ORFs: results_demo

| File | What it is |
|---|---|
| `orfs.csv` | one row per ORF in the provided GTF |
| `summary.csv` | counts per ORF type |
| reviewed_smorfs .gtf / .nuc / .pep | smORFs whose protein is exactly a Swiss-Prot sequence (blastoff sp_identical): ORF structure, CDS nt incl. stop, protein |
| unreviewed_smorfs .gtf / .nuc / .pep | the same for everything else (novel, TrEMBL incl. exact, near-identical and other non-identical Swiss-Prot) |
| subfolders | the full per-ORF table (orfs_annotated.tsv) and every step's files |

## Headline

7,264 ORFs, 7,212 of them smORFs (<= 150 aa).

Provided mode: the ORFs come from a GTF, not from Ribo-seq. They were not called, not scored for translation (no PIF, Uniformity or Drop-off) and not collapsed by protein similarity, so they are neither Primary nor Comprehensive.

## Columns of orfs.csv

| Column | Meaning |
|---|---|
| `orf_id` | ORF identifier |
| `source_id` | the ORF's id in the provided GTF (its orf_id attribute, else transcript_id); identical structures listed together, comma-separated |
| `ucsc` | UCSC genome-browser link, opens on the ORF (blank for scaffolds) |
| `chrom` | chromosome |
| `start` | start, 0-based (BED) |
| `end` | end, half-open (BED) |
| `strand` | strand |
| `host_gene_name` | gene the ORF sits in |
| `gene_id` | Ensembl gene id |
| `smorf_type` | ORF type: the Annotator's call, relabeled from the UniProt match (annotated_CDS; Trunc-Iso: shorter, no residues of its own; N-Iso / C-Iso / NC-Iso: own residues at that end; Ins-Iso / Del-Iso: same ends, an extra / skipped internal segment; TrEMBL: identical to a TrEMBL protein). annotated_CDS / TrEMBL only when the protein is exactly that sequence; near-identical and other hits keep the Annotator's type |
| `general_smorf_type` | smorf_type grouped by transcript architecture: Upstream (uORF, uoORF, uaORF, uaoORF, udORF), Internal (iORF, riORF), Isoform (*-Iso), Downstream (dORF, doORF, daORF, daoORF), Non-coding (lncRNA, intergenic, psORF, aORF, eORF, UA, not_annotated_*), Reviewed (annotated_CDS), TrEMBL |
| `is_smorf` | True if 150 aa or shorter |
| `aa_length` | protein length, aa |
| `blastoff_class` | UniProt match: sp_/tr_ + identical (the protein is exactly that sequence), near_identical (>= blastoff_min_pident end to end, not exact), N_truncated, C_truncated, internal, N_/C_/NC_extended, insertion, deletion, partial, homolog; or novel |
| `blastoff_acc` | matched UniProt accession |
| `blastoff_gene` | matched UniProt gene |
| `shortstop_class` | ShortStop class: sam_secreted, sam_intracellular, prisms; or not_scored_length / not_scored_non_AUG / not_scored_sequence. Model is human-trained |
| `shortstop_probability` | ShortStop's probability for its class (blank if not scored) |
| `shortstop_sam_consensus` | highest SAM-Consensus level reached: SAM-Consensus-100 (all 100 ShortStop consensus models vote SAM), -90 (>= 90%), -50 (> 50%), or none. Nested; -90 or -100 is the most robust set (blank if not scored) |
| `shortstop_sam_consensus_percentage` | share of the 100 consensus models voting SAM (blank if not scored) |
| `diamond_n_species` | species whose Ensembl transcripts encode the protein: best forward-frame DIAMOND blastx hit covering >= diamond_min_cov of it at tblastn E <= diamond_evalue, against the whole transcriptome or the host gene's orthologs |
| `diamond_species` | those species, in diamond_species order; per-species call, gene, identity, coverage and E-values are in the full table |
| `diamond_n_untestable` | species where even a perfect match could not reach diamond_evalue (too short for the search space): absence there is not evidence |
| `aa_seq` | protein sequence |
