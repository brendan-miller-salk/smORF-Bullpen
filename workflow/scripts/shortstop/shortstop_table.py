#!/usr/bin/env python3
"""ShortStop predictions -> annotations/shortstop.tsv (one row per ORF, keyed by orf_id).

Every ORF gets a shortstop_class: ShortStop's call (sam_intracellular,
sam_secreted, prisms) or the reason it was not scored. Also checks that the
protein ShortStop built for each ORF equals bullpen's own translation, and
writes the comparison, split by strand and exon count, to --report.
"""
import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orfs", required=True, help="orf_candidates.bed12")
    ap.add_argument("--predictions", required=True, help="predictions/shortstop_classifications.csv")
    ap.add_argument("--sequences", required=True, help="sequences/unknown_sequences.csv")
    ap.add_argument("--expected", required=True)
    ap.add_argument("--skipped", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()

    orfs = pd.read_csv(a.orfs, sep="\t", header=None, usecols=[3], names=["orf_id"]).drop_duplicates()
    pred = pd.read_csv(a.predictions).rename(columns={
        "classification": "shortstop_class", "probability": "shortstop_probability"})
    skipped = pd.read_csv(a.skipped, sep="\t")
    expected = pd.read_csv(a.expected, sep="\t")

    out = orfs.merge(pred[["orf_id", "shortstop_class", "shortstop_probability"]], on="orf_id", how="left")
    reason = out["orf_id"].map(skipped.set_index("orf_id")["reason"])
    out["shortstop_class"] = out["shortstop_class"].fillna(reason)
    lost = out["shortstop_class"].isna()
    out.loc[lost, "shortstop_class"] = "not_scored_by_shortstop"
    out.to_csv(a.out, sep="\t", index=False)

    # Did ShortStop translate what we submitted?
    seqs = pd.read_csv(a.sequences, usecols=["orf_id", "aa_seq"]).rename(columns={"aa_seq": "shortstop_aa"})
    cmp = expected.merge(seqs, on="orf_id", how="left")
    cmp["exons"] = cmp["n_exons"].map(lambda n: "multi" if n > 1 else "single")
    cmp["status"] = "match"
    cmp.loc[cmp["shortstop_aa"].isna(), "status"] = "no_sequence"
    cmp.loc[cmp["shortstop_aa"].notna() & (cmp["shortstop_aa"].astype(str).str.rstrip("*")
                                           != cmp["aa_seq"]), "status"] = "mismatch"
    table = cmp.groupby(["strand", "exons", "status"]).size().rename("n").reset_index()
    table.to_csv(a.report, sep="\t", index=False)

    print("shortstop classes:", file=sys.stderr)
    for k, v in out["shortstop_class"].value_counts().items():
        print(f"  {k:28s} {v:>9,}", file=sys.stderr)
    print("\nprotein check (bullpen vs ShortStop translation):", file=sys.stderr)
    print(table.to_string(index=False), file=sys.stderr)

    # Gate: every submitted ORF must come back, translated to the protein we
    # expect. Anything else means ShortStop scored a different sequence (or
    # silently dropped one), so its calls cannot be trusted; fail the rule.
    bad = cmp[cmp["status"] != "match"]
    if len(bad) or lost.any():
        print(f"\nERROR: protein check failed for {len(bad):,} of {len(cmp):,} submitted ORFs "
              f"({int(lost.sum()):,} without a call). First 20:", file=sys.stderr)
        cols = ["orf_id", "strand", "n_exons", "status", "aa_seq", "shortstop_aa"]
        print(bad[cols].head(20).to_string(index=False), file=sys.stderr)
        print("A split by strand/exons points at CDS row order in the input GTF; see "
              "shortstop_input.py.", file=sys.stderr)
        sys.exit(1)
    print(f"\nprotein check passed: all {len(cmp):,} submitted ORFs match.", file=sys.stderr)


if __name__ == "__main__":
    main()
