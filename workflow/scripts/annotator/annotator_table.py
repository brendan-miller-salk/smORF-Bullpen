#!/usr/bin/env python3
"""Annotator CSV -> annotations/annotator.tsv (one row per ORF, keyed by orf_id).

The Annotator keys its output on the renamed gene_id, so the id map from
bed12_to_annotator_gtf.py is used to get back to orf_id. ORFs the Annotator
could not see (scaffold seqnames it filters out) are kept with
smorf_type = not_annotated_scaffold, so no ORF silently drops out.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orfs", required=True, help="orf_candidates.bed12")
    p.add_argument("--annotation", required=True, help="Annotator CSV output")
    p.add_argument("--idmap", required=True)
    p.add_argument("--skipped", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)

    orfs = pd.read_csv(args.orfs, sep="\t", header=None, usecols=[3], names=["orf_id"])
    orfs = orfs.drop_duplicates()
    idmap = pd.read_csv(args.idmap, sep="\t")
    ann = pd.read_csv(args.annotation).rename(columns={
        "gene_id": "annotator_gene_id",
        "annotation": "smorf_type",
        "gene_name": "host_gene_name",
        "ensembl_gene_id": "host_gene_id",
    })

    out = (orfs
           .merge(idmap, on="orf_id", how="left")
           .merge(ann[["annotator_gene_id", "smorf_type", "host_gene_name", "host_gene_id"]],
                  on="annotator_gene_id", how="left"))

    skipped = set(pd.read_csv(args.skipped, sep="\t")["orf_id"])
    on_scaffold = out["orf_id"].isin(skipped)
    out.loc[on_scaffold & out["smorf_type"].isna(), "smorf_type"] = "not_annotated_scaffold"

    # An ORF the Annotator saw but never emitted a row for. Should be empty;
    # if it is not, something upstream is dropping ORFs and we want to know.
    still_missing = out["smorf_type"].isna()
    if still_missing.any():
        out.loc[still_missing, "smorf_type"] = "not_annotated_unknown"
        print(f"  WARNING: {int(still_missing.sum()):,} ORFs got no annotation and were not "
              f"on a scaffold. Marked 'not_annotated_unknown' - worth investigating.",
              file=sys.stderr)

    out = out.drop(columns=["annotator_gene_id"])
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(out):,} rows to {args.out}\n\nsmORF types:", file=sys.stderr)
    for t, n in out["smorf_type"].value_counts().items():
        print(f"  {t:26s} {n:>9,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
