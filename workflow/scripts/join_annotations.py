#!/usr/bin/env python3
"""Left-join every annotator's per-ORF table onto the scored ORF table.

Annotator contract: a TSV whose first column is orf_id, one row per ORF in
orf_candidates.bed12. Each tool's columns are kept as the tool named them; a
column name already present is an error rather than a silent overwrite. The
output always has exactly the rows of the scores table.

general_smorf_type groups the final smorf_type by transcript architecture (GENERAL_TYPE):
Upstream, Internal, Isoform, Downstream, Non-coding, Reviewed, TrEMBL. It is added next to
smorf_type whenever smorf_type exists; a label with no group stops the join.

One cross-tool step: with both annotator and blastoff joined, the Annotator's
isoform-blind iORF and overlap calls are relabeled annotated_CDS / Trunc-Iso / N-, C-,
NC-Iso / Ins-Iso / Del-Iso / TrEMBL from the
UniProt class (uniprot/blastoff.py:relabel_smorf_type). The original
call is kept in smorf_type_annotator.
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "blastoff"))
from blastoff import relabel_smorf_type  # noqa: E402

# Why a type is excluded. Trunc-Iso has no sequence of its own: every P-site is the host
# CDS's, so it passes the scores on borrowed reads (99% of heart's Primary Trunc-Isos).
TYPE_REASON = {"Trunc-Iso": "Trunc-Iso: no own sequence; every P-site is shared with the host CDS"}


# smorf_type -> architecture group. Annotator labels (smorf_annotator.py):
# a = upstream/downstream of the CDS in only some isoforms, o = also overlaps the CDS,
# riORF on a retained-intron transcript, aORF on an NMD / non-stop / CDS-not-defined
# transcript, eORF on an exon of another gene type, UA unassigned. Relabels from UniProt
# (blastoff.relabel_smorf_type): annotated_CDS (Swiss-Prot identical), TrEMBL (TrEMBL
# identical), the *-Iso isoforms. annotator_table.py adds the not_annotated_* labels.
GENERAL_TYPE = {
    "Upstream": ["uORF", "uoORF", "uaORF", "uaoORF", "udORF"],
    "Internal": ["iORF", "riORF"],
    "Isoform": ["N-Iso", "C-Iso", "NC-Iso", "Ins-Iso", "Del-Iso", "Trunc-Iso"],
    "Downstream": ["dORF", "doORF", "daORF", "daoORF"],
    "Non-coding": ["lncRNA", "intergenic", "psORF", "aORF", "eORF", "UA",
                   "not_annotated_scaffold", "not_annotated_unknown"],
    "Reviewed": ["annotated_CDS"],
    "TrEMBL": ["TrEMBL"],
}
GENERAL_OF = {t: g for g, types in GENERAL_TYPE.items() for t in types}


def general_smorf_type(types):
    """smorf_type values -> general_smorf_type list (blank for a missing type). A label
    with no group is an error, so a new Annotator label can't silently lose its row."""
    out, unknown = [], set()
    for t in types:
        if not isinstance(t, str) or not t:
            out.append("")
        elif t in GENERAL_OF:
            out.append(GENERAL_OF[t])
        else:
            unknown.add(t)
            out.append("")
    if unknown:
        raise SystemExit(f"join_annotations: smorf_type with no general_smorf_type group: "
                         f"{sorted(unknown)}; add it to GENERAL_TYPE")
    return out


def exclude_types(df, types):
    """Move ORFs of the given smorf_types out of Comprehensive and Primary: exclude_reason
    set (an earlier reason is kept) and pass_all False. Returns how many were moved."""
    if not types or "smorf_type" not in df or "pass_all" not in df:
        return 0
    hit = df["smorf_type"].isin(types)
    if "exclude_reason" not in df:
        df["exclude_reason"] = ""
    blank = df["exclude_reason"].fillna("").astype(str) == ""
    df.loc[hit & blank, "exclude_reason"] = df.loc[hit & blank, "smorf_type"].map(
        lambda t: TYPE_REASON.get(t, f"{t}: excluded type (exclude_types)"))
    df.loc[hit, "pass_all"] = False
    return int(hit.sum())


# Reviewed (blastoff sp_identical: the protein is exactly a Swiss-Prot sequence) is the known
# protein itself, so it has no sequence of its own by construction and is exempt from the
# own-sequence and RPKM rules. Exact TrEMBL matches are Unreviewed and go through both
# (~97% are annotated Ensembl CDS with no sequence of their own) unless
# config trembl_identical_exempt adds tr_identical here (--exempt-classes).
EXEMPT_CLASSES = {"sp_identical"}


def exclude_unique(df, types, min_inframe):
    """ORFs of the given types need >= min_inframe in-frame P-sites from uniquely mapped
    reads in their unshared sequence; the rest move to excluded. Meant for psORF: a
    pseudogene nearly identical to its parent otherwise counts the parent's reads. The
    type is tested on the Annotator's call too (smorf_type_annotator), so a pseudogene
    copy the UniProt relabel turned into annotated_CDS / an isoform is still tested."""
    need = {"smorf_type", "pass_all", "unshared_inframe_psites_unique"}
    if not types or not need <= set(df.columns):
        return 0
    is_type = df["smorf_type"].isin(types)
    if "smorf_type_annotator" in df:
        is_type |= df["smorf_type_annotator"].isin(types)
    hit = is_type & (df["unshared_inframe_psites_unique"].fillna(0) < min_inframe)
    if "exclude_reason" not in df:
        df["exclude_reason"] = ""
    blank = df["exclude_reason"].fillna("").astype(str) == ""
    orig = df["smorf_type_annotator"] if "smorf_type_annotator" in df else df["smorf_type"]
    df.loc[hit & blank, "exclude_reason"] = orig[hit & blank].map(
        lambda t: f"{t}: fewer than {min_inframe} in-frame P-sites from uniquely mapped reads "
                  "outside annotated CDS")
    df.loc[hit, "pass_all"] = False
    return int(hit.sum())


def exclude_no_unshared(df, exempt=EXEMPT_CLASSES):
    """ORFs with no sequence outside annotated CDS (unshared_nt 0) move to excluded unless
    their blastoff_class is exempt (default EXEMPT_CLASSES: Reviewed, the known Swiss-Prot protein). Their every P-site is
    shared with a CDS, so they have no evidence of their own; score_orfs' unshared rule
    only judges ORFs that have unshared sequence. Returns how many were moved."""
    need = {"unshared_nt", "blastoff_class", "pass_all"}
    if not need <= set(df.columns):
        print(f"  WARNING: --require-unshared needs {sorted(need - set(df.columns))}; skipped",
              file=sys.stderr)
        return 0
    hit = (df["unshared_nt"].fillna(0) == 0) & ~df["blastoff_class"].isin(set(exempt))
    if "exclude_reason" not in df:
        df["exclude_reason"] = ""
    blank = df["exclude_reason"].fillna("").astype(str) == ""
    label = df["smorf_type"] if "smorf_type" in df else pd.Series("ORF", index=df.index)
    df.loc[hit & blank, "exclude_reason"] = label[hit & blank].map(
        lambda t: f"{t}: no sequence outside annotated CDS")
    df.loc[hit, "pass_all"] = False
    return int(hit.sum())


def exclude_chroms(df, chroms):
    """ORFs on the given chromosomes move to excluded. chrM: the callers call and
    translate with the standard genetic code, so a chrM ORF's structure and protein
    are not the mitochondrial ones. Returns how many were moved."""
    if not chroms or "chrom" not in df or "pass_all" not in df:
        return 0
    hit = df["chrom"].isin(chroms)
    if "exclude_reason" not in df:
        df["exclude_reason"] = ""
    blank = df["exclude_reason"].fillna("").astype(str) == ""
    df.loc[hit & blank, "exclude_reason"] = df.loc[hit & blank, "chrom"].map(
        lambda c: f"{c}: called with the standard genetic code")
    df.loc[hit, "pass_all"] = False
    return int(hit.sum())


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scores", required=True, help="orf_scores.tsv")
    p.add_argument("--annotations", nargs="*", default=[], help="annotations/<tool>.tsv")
    p.add_argument("--out", required=True)
    p.add_argument("--unique-types", nargs="*", default=[],
                   help="smorf_types that need unique-read evidence (config unique_evidence_types)")
    p.add_argument("--unique-min-inframe", type=int, default=5)
    p.add_argument("--exclude-types", nargs="*", default=[],
                   help="smorf_types to take out of Comprehensive and Primary (config exclude_types)")
    p.add_argument("--exclude-chroms", nargs="*", default=[],
                   help="chromosomes whose ORFs move to excluded (config exclude_chroms, e.g. chrM)")
    p.add_argument("--require-unshared", action="store_true",
                   help="exclude ORFs with no unshared sequence unless sp_identical (Reviewed) "
                        "(config require_unshared_sequence)")
    p.add_argument("--exempt-classes", nargs="*", default=sorted(EXEMPT_CLASSES),
                   help="blastoff classes exempt from --require-unshared (config trembl_identical_exempt "
                        "adds tr_identical)")
    args = p.parse_args(argv)

    merged = pd.read_csv(args.scores, sep="\t")
    n = len(merged)
    for path in args.annotations:
        tool = os.path.basename(path).removesuffix(".tsv")
        ann = pd.read_csv(path, sep="\t")
        if ann.columns[0] != "orf_id":
            sys.exit(f"join_annotations: {path}: first column must be orf_id, got {ann.columns[0]!r}")
        if ann["orf_id"].duplicated().any():
            sys.exit(f"join_annotations: {path}: duplicate orf_id rows")
        clash = set(ann.columns[1:]) & set(merged.columns)
        if clash:
            sys.exit(f"join_annotations: {tool} columns already present: {sorted(clash)}")
        missing = ~merged["orf_id"].isin(ann["orf_id"])
        if missing.any():
            print(f"  WARNING: {tool}: {int(missing.sum()):,} scored ORFs have no row",
                  file=sys.stderr)
        merged = merged.merge(ann, on="orf_id", how="left")
        print(f"  joined {tool}: +{len(ann.columns) - 1} columns", file=sys.stderr)
    assert len(merged) == n, "row count changed during join"

    if {"smorf_type", "blastoff_class"} <= set(merged.columns):
        merged.insert(merged.columns.get_loc("smorf_type") + 1, "smorf_type_annotator",
                      merged["smorf_type"])
        div = merged["blastoff_diverges"] if "blastoff_diverges" in merged else [""] * n
        merged["smorf_type"] = [relabel_smorf_type(t, c, d) for t, c, d in
                                zip(merged["smorf_type_annotator"], merged["blastoff_class"], div)]
        changed = merged[merged["smorf_type"] != merged["smorf_type_annotator"]]
        moved = changed.groupby(["smorf_type_annotator", "smorf_type"]).size()
        print(f"\n  relabeled from UniProt ({len(changed):,} ORFs):", file=sys.stderr)
        for (a, b), k in moved.sort_values(ascending=False).items():
            print(f"    {a:<8} -> {b:<14} {k:>8,}", file=sys.stderr)

    if "smorf_type" in merged.columns:
        merged.insert(merged.columns.get_loc("smorf_type") + 1, "general_smorf_type",
                      general_smorf_type(merged["smorf_type"]))
        g = merged["general_smorf_type"].value_counts()
        print("\n  general_smorf_type: " + ", ".join(f"{k} {v:,}" for k, v in g.items()),
              file=sys.stderr)

    u = exclude_unique(merged, args.unique_types, args.unique_min_inframe)
    if u:
        print(f"\n  excluded for unique-read evidence ({', '.join(args.unique_types)}): {u:,} ORFs",
              file=sys.stderr)
    k = exclude_types(merged, args.exclude_types)
    if k:
        print(f"\n  excluded by type ({', '.join(args.exclude_types)}): {k:,} ORFs", file=sys.stderr)
    if args.require_unshared:
        o = exclude_no_unshared(merged, args.exempt_classes)
        print(f"\n  excluded for no sequence outside annotated CDS (not {'/'.join(args.exempt_classes)}): "
              f"{o:,} ORFs", file=sys.stderr)
    c = exclude_chroms(merged, args.exclude_chroms)
    if c:
        print(f"\n  excluded on {', '.join(args.exclude_chroms)}: {c:,} ORFs", file=sys.stderr)

    merged.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {n:,} rows to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
