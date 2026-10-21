# Provided-mode demo: GENCODE Ribo-seq ORFs

This demo runs the 7,264 GENCODE Ribo-seq ORFs (GRCh38; Mudge et al. 2022,
[doi:10.1038/s41587-022-01369-0](https://doi.org/10.1038/s41587-022-01369-0)) through
provided mode: the Annotator, ShortStop, blastoff and DIAMOND, with no Ribo-seq. On
CARC's `epyc-64` partition the run took 19 minutes.

30 of these ORFs end in a stop codon split by an intron, written as two `stop_codon` rows.
The catalogue reads the stop across the junction, so all 7,264 ORFs are kept.

## Files

| File | Contents |
|---|---|
| `smORFsGencode.hg38.gtf.gz` | The ORFs: `transcript`, `exon`, `CDS`, `start_codon` and `stop_codon` rows |
| `config.yaml` | Run config. Replace `/path/to/refs` with your references |
| `expected/` | Outputs of the run below, for comparison |

## References

| Config key | File |
|---|---|
| `genome` | GENCODE v50 `GRCh38.primary_assembly.genome.fa` |
| `gtf` | Ensembl 116 GRCh38 GTF, chromosomes renamed to `chr` (the `rename_gtf` step of `refs/fetch_references.sh`) |
| `blastoff_fasta` | UniProt human proteome UP000005640 (release 2026_03), canonical and isoforms |
| `diamond_species`, `diamond_orthologs` | `SOURCE=human:hsapiens RELEASE=116 bash refs/fetch_conservation.sh /path/to/refs` |

`refs/fetch_references.sh` fetches mouse only, so the human genome, GTF and proteome have to be
downloaded by hand.

## Run

```bash
mkdir -p ~/runs/bullpen_demo/config ~/runs/bullpen_demo/logs && cd ~/runs/bullpen_demo
cp $BULLPEN_REPO/demo/provided_mode/config.yaml config/demo.yaml      # then edit the paths
gunzip -c $BULLPEN_REPO/demo/provided_mode/smORFsGencode.hg38.gtf.gz > smORFsGencode.hg38.gtf
sbatch -A <account> -p <partition> $BULLPEN_REPO/profiles/slurm/run_head.sbatch config/demo.yaml -n
```

The dry run should list 29 jobs and none of `trim`, `star_align`, `price`, `ribocode` or
`score`. Drop `-n` to run it.

## Expected output

`catalogue/catalogue_report.tsv`: 7,264 transcripts with CDS rows. All 7,264 end in a stop,
no structures are duplicated, and 7,264 ORFs are written. `catalogue_rejects.tsv` holds only
its header.

`orfs.csv` has 7,264 rows. 7,212 are smORFs (at most 150 aa).

| `general_smorf_type` (smORFs) | n |
|---|---|
| Upstream | 4,071 |
| Non-coding | 1,868 |
| Downstream | 811 |
| Internal | 378 |
| Isoform | 37 |
| TrEMBL | 28 |
| Reviewed | 19 |

The 19 Reviewed smORFs are exactly a Swiss-Prot protein and are written to
`reviewed_smorfs.{gtf,nuc,pep}`. The other 7,193 are in `unreviewed_smorfs.*`. Counts for
every `smorf_type` are in `expected/summary.csv`.

| ShortStop (all ORFs) | n |
|---|---|
| `prisms` | 6,598 |
| `sam_intracellular` | 578 |
| `sam_secreted` | 36 |
| `not_scored_length` | 52 |

| DIAMOND: species with a conserved hit (`diamond_n_species`) | n |
|---|---|
| 0 | 4,511 |
| 1 | 869 |
| 2 | 542 |
| 3 | 916 |
| 4 | 274 |
| 5 (mouse, rat, cow, chicken, zebrafish) | 152 |

Counts can shift a little with other Ensembl, UniProt or ShortStop releases.
