#!/usr/bin/env python3
"""Fold smORFs that MMseqs2 clusters together into one catalogue row.

Rules:
  - only is_smorf rows are grouped, by MMseqs2 cluster; groups of >= 2 fold
  - survivor = best class, then longest aa_length, then orf_id; it keeps its
    own orf_id and coordinates (ids are not renumbered)
  - called_by_* OR-ed, pval_* best (min), samples and orf_type_native unioned
  - dropped members vanish from bed12/tsv; consensus is re-applied afterwards
In addition, `collapsed_orf_ids` lists the members folded into each survivor,
so the collapse stays traceable; and only members at
one locus (same chrom and strand, overlapping spans) and with the same splice
junctions where they overlap fold together. Folding across loci put a
processed pseudogene (Dnajc19-ps, chr7) into its parent gene's row
(Dnajc19, chr3); and folding on clusters alone re-folded splice isoforms the merge keeps
apart (797 of heart's 1,406 folds, incl. the spliced Dnajc19-ps N-isoform).
"""
import argparse
import collections
import sys

import pandas as pd

from build_catalogue import same_junctions

CLASS_RANK = {c: i for i, c in enumerate(
    ["canonical_cds", "uoORF", "uORF", "doORF", "dORF", "intORF", "novel_u", "other"])}


def strip_strand(x):
    return x[:-3] if x.endswith(("(+)", "(-)")) else x


def union(values):
    return ",".join(sorted({v for s in values for v in str(s).split(",") if v and v != "nan"}))


def bed_blocks(line):
    f = line.rstrip("\n").split("\t")
    s = int(f[1])
    sizes = [int(x) for x in f[10].rstrip(",").split(",")]
    starts = [int(x) for x in f[11].rstrip(",").split(",")]
    return tuple(sorted((s + a, s + a + b) for a, b in zip(starts, sizes)))


def same_structure_groups(rows, blocks, order):
    """Greedy complete linkage: a row joins a group only if its splice junctions match
    every member's where they overlap. `order` puts would-be survivors first."""
    groups = []
    for i in sorted(rows, key=order):
        for g in groups:
            if all(same_junctions(blocks[i], blocks[j]) for j in g):
                g.append(i)
                break
        else:
            groups.append([i])
    return groups


def same_locus_groups(df, idx):
    """Split one MMseqs2 cluster into groups of rows whose spans overlap on one strand."""
    groups = []
    for i in sorted(idx, key=lambda i: (df.at[i, "chrom"], df.at[i, "strand"], int(df.at[i, "start"]))):
        c, st, lo, hi = df.at[i, "chrom"], df.at[i, "strand"], int(df.at[i, "start"]), int(df.at[i, "end"])
        g = groups[-1] if groups else None
        if g and g["chrom"] == c and g["strand"] == st and lo < g["end"]:
            g["rows"].append(i)
            g["end"] = max(g["end"], hi)
        else:
            groups.append(dict(chrom=c, strand=st, end=hi, rows=[i]))
    return [g["rows"] for g in groups]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True, help="merged catalogue TSV")
    ap.add_argument("--bed12", required=True, help="merged catalogue BED12")
    ap.add_argument("--clusters", required=True, help="MMseqs2 <prefix>_cluster.tsv (rep, member)")
    ap.add_argument("--smorf-max-aa", type=int, required=True)
    ap.add_argument("--min-callers", type=int, default=1)
    ap.add_argument("--min-samples", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    df = pd.read_csv(a.tsv, sep="\t", dtype=str, keep_default_na=False)
    expected = ((df["aa_length"].astype(int) > 0) & (df["aa_length"].astype(int) <= a.smorf_max_aa)).astype(int)
    if not (df["is_smorf"].astype(int) == expected).all():
        sys.exit("collapse_catalogue: is_smorf disagrees with aa_length and --smorf-max-aa")
    callers = [c.removeprefix("called_by_") for c in df.columns if c.startswith("called_by_")]

    cluster_of = {}
    for line in open(a.clusters):
        rep, member = line.rstrip("\n").split("\t")[:2]
        cluster_of[strip_strand(member)] = strip_strand(rep)

    groups = collections.defaultdict(list)
    for i, r in df[df["is_smorf"] == "1"].iterrows():
        groups[cluster_of.get(r["orf_id"], r["orf_id"])].append(i)

    bed = {l.split("\t")[3]: l for l in open(a.bed12)}
    blocks = {i: bed_blocks(bed[df.at[i, "orf_id"]]) for i in df.index[df["is_smorf"] == "1"]}

    def survivor_rank(i):
        return (CLASS_RANK.get(df.at[i, "orf_class"], len(CLASS_RANK)), -int(df.at[i, "aa_length"]), df.at[i, "orf_id"])

    drop, n_folded = set(), 0
    df["collapsed_orf_ids"] = ""
    for idx in (s for c in groups.values() for g in same_locus_groups(df, c)
                for s in same_structure_groups(g, blocks, survivor_rank)):
        if len(idx) < 2:
            continue
        sub = df.loc[idx]
        order = sorted(idx, key=lambda i: (CLASS_RANK.get(df.at[i, "orf_class"], len(CLASS_RANK)),
                                           -int(df.at[i, "aa_length"]), df.at[i, "orf_id"]))
        keep, rest = order[0], order[1:]
        for c in callers:
            df.at[keep, f"called_by_{c}"] = "1" if (sub[f"called_by_{c}"] == "1").any() else "0"
            p = pd.to_numeric(sub[f"pval_{c}"], errors="coerce").dropna()
            df.at[keep, f"pval_{c}"] = f"{p.min():.6g}" if len(p) else ""
        df.at[keep, "n_callers"] = str(sum(df.at[keep, f"called_by_{c}"] == "1" for c in callers))
        df.at[keep, "samples"] = union(sub["samples"])
        df.at[keep, "n_samples"] = str(len([s for s in df.at[keep, "samples"].split(",") if s]))
        df.at[keep, "orf_type_native"] = union(sub["orf_type_native"])
        df.at[keep, "collapsed_orf_ids"] = ",".join(sorted(df.at[i, "orf_id"] for i in rest))
        drop.update(rest)
        n_folded += 1

    out = df.drop(index=sorted(drop))
    keep_ids = set(out["orf_id"])
    out.to_csv(a.out_prefix + ".tsv", sep="\t", index=False)

    cons = out[(out["n_callers"].astype(int) >= a.min_callers)
               & ((out["n_samples"].astype(int) >= a.min_samples) | (out.get("called_by_price", "0") == "1"))]
    with open(a.out_prefix + ".bed12", "w") as fh:
        fh.writelines(bed[o] for o in out["orf_id"])
    with open(a.out_prefix + ".consensus.bed12", "w") as fh:
        fh.writelines(bed[o] for o in cons["orf_id"])

    print(f"collapse_catalogue: {len(df):,} -> {len(out):,} ORFs; {n_folded:,} clusters folded "
          f"{len(drop):,} smORFs; {len(cons):,} consensus", file=sys.stderr)
    assert keep_ids <= set(bed), "survivor missing from BED12"


if __name__ == "__main__":
    main()
