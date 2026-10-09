#!/usr/bin/env python3
"""orf_candidates.bed12 -> ShortStop predict GTF, plus the ORFs it would skip.

ShortStop's input contract (converters/gtf_to_seq.py):
  - 9 columns, no header or '#' lines
  - per ORF one `transcript` row + `CDS` rows; the id is gene_id, unquoted
  - CDS rows ascending by coordinate, each with a distinct attribute string
  - CDS rows EXCLUDE the stop codon (the 3' flank starts after it)
  - the transcript row spans the CDS +/- the 25 nt flanks, else they pad with X
ShortStop silently drops proteins that do not start with M, are outside 9-150
aa, or contain non-standard residues. Those ORFs are withheld here and given an
explicit reason instead, so none vanish without one.

Also writes the protein bullpen expects for each submitted ORF, so the
ShortStop step can confirm it translated the same sequence.
"""
import argparse
import sys

from pyfaidx import Fasta

COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
STOPS = {"TAA", "TAG", "TGA"}
FLANK = 25


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed12", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--min-aa", type=int, default=9)
    ap.add_argument("--max-aa", type=int, default=150)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-skipped", required=True, help="orf_id, reason")
    ap.add_argument("--out-expected", required=True, help="orf_id, aa_seq, strand, n_exons")
    a = ap.parse_args()

    genome = Fasta(a.genome, as_raw=True, sequence_always_upper=True)
    n_sub = 0
    reasons = {}
    with open(a.bed12) as bed, open(a.out_gtf, "w") as gtf, \
            open(a.out_skipped, "w") as skip, open(a.out_expected, "w") as exp:
        skip.write("orf_id\treason\n")
        exp.write("orf_id\taa_seq\tstrand\tn_exons\n")
        for line in bed:
            f = line.rstrip("\n").split("\t")
            chrom, start, strand, oid = f[0], int(f[1]), f[5], f[3]
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            starts = [int(x) for x in f[11].rstrip(",").split(",")]
            blocks = [(start + s, start + s + n) for s, n in zip(starts, sizes)]

            seq = "".join(genome[chrom][s:e] for s, e in blocks)
            if strand == "-":
                seq = seq.translate(COMP)[::-1]
            # Trim the stop so CDS rows exclude it, as ShortStop expects.
            if seq[-3:] in STOPS:
                if strand == "+":
                    blocks[-1] = (blocks[-1][0], blocks[-1][1] - 3)
                else:
                    blocks[0] = (blocks[0][0] + 3, blocks[0][1])
                blocks = [(s, e) for s, e in blocks if e > s]
                seq = seq[:-3]
            aa = translate(seq)

            if seq[:3] != "ATG":
                reason = "not_scored_non_AUG"
            elif not a.min_aa <= len(aa) <= a.max_aa:
                reason = "not_scored_length"
            elif set(aa) - set("ACDEFGHIKLMNPQRSTVWY"):
                reason = "not_scored_sequence"
            else:
                reason = None
            if reason:
                skip.write(f"{oid}\t{reason}\n")
                reasons[reason] = reasons.get(reason, 0) + 1
                continue

            clen = len(genome[chrom])
            t0 = max(0, blocks[0][0] - FLANK - 3)
            t1 = min(clen, blocks[-1][1] + FLANK + 3)
            attr = f"gene_id {oid}; transcript_id {oid};"
            gtf.write("\t".join([chrom, "bullpen", "transcript", str(t0 + 1), str(t1), ".",
                                 strand, ".", attr]) + "\n")
            # CDS rows ascending by coordinate, each with its own exon_number, as a
            # typical GTF has them: ShortStop keeps that order and assembles minus-strand
            # exons correctly from it. Identical attributes on every row trigger its
            # per-group re-sort instead, and spliced minus-strand ORFs come out 3' first.
            for i, (s, e) in enumerate(blocks, 1):
                gtf.write("\t".join([chrom, "bullpen", "CDS", str(s + 1), str(e), ".",
                                     strand, "0", f"{attr} exon_number {i};"]) + "\n")
            exp.write(f"{oid}\t{aa}\t{strand}\t{len(blocks)}\n")
            n_sub += 1

    print(f"shortstop_input: {n_sub:,} ORFs submitted; withheld {reasons}", file=sys.stderr)
    if n_sub == 0:
        sys.exit("shortstop_input: no ORF is eligible for ShortStop")


if __name__ == "__main__":
    main()
