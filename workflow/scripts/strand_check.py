#!/usr/bin/env python3
"""Confirm every library is sense-stranded before anything uses its reads.

Bullpen assumes sense reads everywhere:
- P-site scoring reads the fw track for + ORFs and the rc track for - ORFs.
- PRICE and RiboCode are given no strand option.
A reverse-stranded or unstranded library would therefore score as noise
without any error.

featureCounts counts each library's uniquely mapped reads (multimappers are
excluded by default) on annotated CDS twice: -s 1 (read on the feature's strand)
and -s 2 (opposite). sense_fraction = s1 / (s1 + s2). A sample is
  forward     sense_fraction >= --min-sense
  reverse     sense_fraction <= 1 - --min-sense
  ambiguous   anything between (unstranded, or mixed)
  none        no CDS reads at all
and PASSes only as forward with at least --min-assigned CDS reads. With
--strict, any FAIL exits 1, naming the samples.
"""
import argparse
import os
import sys

import pandas as pd

SUFFIX = ".genome.sorted.bam"


def assigned(summary):
    """featureCounts .summary -> {sample: Assigned count}, samples named from the BAM file."""
    df = pd.read_csv(summary, sep="\t", index_col=0)
    row = df.loc["Assigned"]
    return {os.path.basename(c).removesuffix(SUFFIX): int(v) for c, v in row.items()}


def call(sense, anti, min_sense, min_assigned):
    """(sense_fraction, call, status, reason) for one sample."""
    total = sense + anti
    frac = sense / total if total else 0.0
    if not total:
        c = "none"
    elif frac >= min_sense:
        c = "forward"
    elif frac <= 1 - min_sense:
        c = "reverse"
    else:
        c = "ambiguous"
    why = []
    if c != "forward":
        why.append(f"{c}: sense_fraction {frac:.3f} (forward needs >= {min_sense})")
    if total < min_assigned:
        why.append(f"{total} CDS reads < {min_assigned}")
    return frac, c, "FAIL" if why else "PASS", "; ".join(why)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sense", required=True, help="featureCounts -s 1 .summary")
    ap.add_argument("--antisense", required=True, help="featureCounts -s 2 .summary")
    ap.add_argument("--min-sense", type=float, default=0.9)
    ap.add_argument("--min-assigned", type=int, default=100000)
    ap.add_argument("--strict", action="store_true", help="exit 1 if any sample FAILs")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    s1, s2 = assigned(a.sense), assigned(a.antisense)
    rows = []
    for sample in sorted(s1):
        frac, c, status, why = call(s1[sample], s2[sample], a.min_sense, a.min_assigned)
        rows.append({"sample": sample, "status": status, "call": c,
                     "sense_reads": s1[sample], "antisense_reads": s2[sample],
                     "sense_fraction": round(frac, 4), "reason": why})
    df = pd.DataFrame(rows)
    df.to_csv(a.out, sep="\t", index=False)
    print(df.to_string(index=False), file=sys.stderr)

    bad = df.loc[df["status"] == "FAIL", "sample"].tolist()
    if bad and a.strict:
        sys.exit(f"strand_check: not sense-stranded or too few CDS reads: {' '.join(bad)} "
                 f"(reasons in {a.out}). Bullpen assumes sense reads throughout; drop these "
                 f"from samples.tsv, or set strand_check: warn to continue anyway.")


if __name__ == "__main__":
    main()
