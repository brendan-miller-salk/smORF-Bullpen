#!/usr/bin/env python3
"""featureCounts matrix -> counts.tsv: one row per canonical gene and per unreviewed smORF.

Columns: feature_id, feature_type (canonical | unreviewed_smorf), length_nt (featureCounts'
union of the feature's CDS rows; for a smORF, its unique_nt), then one count
column per sample, then rpkm_<sample> for each.

RPKM = count / (length_nt / 1e3) / (library total / 1e6). length_nt is what was
counted: a gene's CDS union, or a smORF's unique part. The library total is
every read assigned in this table (canonical + smORF) for that sample, so RPKM is
comparable across the table's features and libraries, not with other runs' ppm.

smORFs whose CDS was fully subtracted have no GTF row, so featureCounts never
saw them; they are added back with length 0, zero counts and a blank RPKM.
"""
import argparse
import os
import sys

import pandas as pd

SUFFIX = ".genome.sorted.bam"
SMORF = "unreviewed_smorf"


def rpkm(counts, length_nt):
    """counts (features x samples) -> RPKM; rows with length 0 are NaN (written blank)."""
    lib = counts.sum(axis=0) / 1e6
    kb = length_nt.where(length_nt > 0) / 1e3
    return counts.div(kb, axis=0).div(lib.where(lib > 0), axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="featureCounts -o output")
    ap.add_argument("--orfs", required=True, help="build_count_gtf.py --out-orfs table")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fc = pd.read_csv(a.counts, sep="\t", comment="#")
    samples = [c for c in fc.columns[6:]]
    fc = fc.rename(columns={c: os.path.basename(c).removesuffix(SUFFIX) for c in samples})
    samples = [os.path.basename(c).removesuffix(SUFFIX) for c in samples]
    orfs = pd.read_csv(a.orfs, sep="\t")

    counted = set(orfs["orf_id"])
    fc["feature_type"] = fc["Geneid"].map(lambda g: SMORF if g in counted else "canonical")
    out = fc.rename(columns={"Geneid": "feature_id", "Length": "length_nt"})[
        ["feature_id", "feature_type", "length_nt"] + samples]
    missing = orfs.loc[~orfs["orf_id"].isin(out["feature_id"]), "orf_id"]
    zero = pd.DataFrame({"feature_id": missing, "feature_type": SMORF, "length_nt": 0,
                         **{s: 0 for s in samples}})
    out = pd.concat([out, zero], ignore_index=True)
    r = rpkm(out[samples], out["length_nt"]).round(4)
    out = pd.concat([out, r.add_prefix("rpkm_")], axis=1)
    out.to_csv(a.out, sep="\t", index=False)

    by = out.groupby("feature_type")[samples].sum()
    detected = (out.loc[out["feature_type"] == SMORF, samples] > 0).sum()
    print(f"counts_table: {(out['feature_type'] == 'canonical').sum()} canonical genes, "
          f"{(out['feature_type'] == SMORF).sum()} unreviewed smORFs ({len(missing)} fully subtracted)",
          file=sys.stderr)
    for s in samples:
        tot = by[s].sum()
        n = by.loc[SMORF, s] if SMORF in by.index else 0
        print(f"  {s}: {tot:,} reads; unreviewed smORFs {n:,} ({n / tot if tot else 0:.2%}); "
              f"smORFs with >= 1 read: {detected[s]}", file=sys.stderr)


if __name__ == "__main__":
    main()
