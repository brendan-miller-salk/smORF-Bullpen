#!/usr/bin/env python3
"""Check pass_all against what Chothani et al. report. Warns; --strict makes a miss fatal.

The thresholds are fitted per dataset and are never copied from the paper. This
only asks whether the result lands where the paper's did:

  threshold_<score>     fitted threshold vs the range of the paper's two
                        datasets (Human body map, RiboCrypt)
  calibration_<score>   share of calibration CDS that clear their own threshold
                        vs the configured percentile (the paper's calibration rule)
  primary_fraction      Primary / Comprehensive among non-CDS ORFs vs the paper's
                        10,127 / 28,359. The paper's set holds non-CDS ORFs only
  length_bias           PIF passes more often for ORFs < 100 nt than > 100 nt
                        (paper: 63.48% vs 37.22%)
  overlap_bias          overlapping types pass_all less often than non-overlapping
                        ones (paper Limitation 2)

Comprehensive = orf_set in {primary, comprehensive} when that column exists (the
master table); otherwise every row, because Bullpen's catalogue is already
caller-filtered. Checks whose columns are missing report SKIP.
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

# Chothani et al.; see riboseq/METHODS_SPEC.md.
PAPER = {
    "pif": (0.5112, 0.7589),          # RiboCrypt, Human body map
    "uniformity": (0.7130, 0.8779),   # Human body map, RiboCrypt
    "dropoff": (0.88, 0.92),
    "primary_fraction": 10127 / 28359,
    "pif_pass_short": 0.6348,
    "pif_pass_long": 0.3722,
}
OVERLAPPING = {"uoORF", "doORF", "iORF", "uaoORF", "daoORF"}
NON_OVERLAPPING = {"uORF", "dORF", "lncRNA", "intergenic", "uaORF", "daORF"}
SCORES = ("pif", "uniformity", "dropoff")


def is_true(s):
    return s.astype(str) == "True"


def run(scores, thresholds, cds=None, catalogue=None, fraction_tol=0.15, calib_tol=0.03):
    """Return a list of dicts: check, observed, expected, status, note."""
    out = []

    def add(check, observed, expected, ok, note=""):
        out.append(dict(check=check, observed=observed, expected=expected,
                        status="SKIP" if ok is None else "OK" if ok else "WARN", note=note))

    thr = thresholds.iloc[0]
    for s in SCORES:
        lo, hi = PAPER[s]
        v = float(thr[s])
        add(f"threshold_{s}", f"{v:.3f}", f"{lo:.3f}-{hi:.3f}", lo <= v <= hi,
            "" if lo <= v <= hi else ("below" if v < lo else "above") + " the Chothani / GENCODE v2 2025 range")

    if cds is not None and "min_psites" in thr:
        used = cds[cds["psites_total"] >= float(thr["min_psites"])]
        want = float(thr["percentile"]) if thr["tail"] == "lower" else 1 - float(thr["percentile"])
        for s in ("pif", "uniformity"):
            got = (used[s] >= float(thr[s])).mean()
            add(f"calibration_{s}", f"{got:.3f}", f"{want:.3f} +/- {calib_tol}",
                abs(got - want) <= calib_tol, f"{len(used):,} CDS")
    else:
        add("calibration", "", "", None, "no cds_calibration table")

    df = scores.copy()
    if catalogue is not None and "orf_class" not in df and "orf_class" in catalogue:
        df = df.merge(catalogue[["orf_id", "orf_class"]], on="orf_id", how="left")
    comp = df[df["orf_set"].isin(["primary", "comprehensive"])] if "orf_set" in df else df
    if "orf_class" in comp:
        nc = comp[comp["orf_class"] != "canonical_cds"]
        prim = is_true(nc["pass_all"]) if "orf_set" not in nc else nc["orf_set"] == "primary"
        f = prim.mean()
        want = PAPER["primary_fraction"]
        add("primary_fraction", f"{f:.3f} ({int(prim.sum()):,}/{len(nc):,})",
            f"{want:.3f} +/- {fraction_tol}", abs(f - want) <= fraction_tol, "non-CDS ORFs")
    else:
        nc = comp
        add("primary_fraction", "", f"{PAPER['primary_fraction']:.3f}", None, "no orf_class")

    if "length_nt" in nc and "pass_pif" in nc:
        short = is_true(nc.loc[nc["length_nt"] < 100, "pass_pif"]).mean()
        long_ = is_true(nc.loc[nc["length_nt"] > 100, "pass_pif"]).mean()
        add("length_bias", f"<100 nt {short:.3f}, >100 nt {long_:.3f}",
            f"short > long ({PAPER['pif_pass_short']:.3f} vs {PAPER['pif_pass_long']:.3f})",
            short > long_, "" if short > long_ else "reversed: long ORFs pass PIF more often")

    tcol = "smorf_type_annotator" if "smorf_type_annotator" in nc else "smorf_type"
    if tcol in nc:
        ov = is_true(nc.loc[nc[tcol].isin(OVERLAPPING), "pass_all"])
        no = is_true(nc.loc[nc[tcol].isin(NON_OVERLAPPING), "pass_all"])
        if len(ov) and len(no):
            add("overlap_bias", f"overlapping {ov.mean():.3f}, non-overlapping {no.mean():.3f}",
                "overlapping < non-overlapping", ov.mean() < no.mean(),
                f"n={len(ov):,} / {len(no):,}")
    else:
        add("overlap_bias", "", "", None, "no smorf_type")
    return out


def report(rows, path, stream=sys.stderr):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    warns = [r for r in rows if r["status"] == "WARN"]
    print(f"pass_all QC vs Chothani / GENCODE v2 2025: {len(rows) - len(warns)} OK/SKIP, {len(warns)} WARN "
          f"-> {path}", file=stream)
    for r in warns:
        print(f"  WARN {r['check']}: {r['observed']} (Chothani / GENCODE v2 2025: {r['expected']}) {r['note']}", file=stream)
    return len(warns)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True, help="orf_scores_annotated.tsv or orf_master.tsv")
    ap.add_argument("--thresholds", required=True, help="thresholds.tsv")
    ap.add_argument("--cds", help="cds_calibration.tsv")
    ap.add_argument("--catalogue", help="catalogue TSV with orf_class, if --scores lacks it")
    ap.add_argument("--fraction-tol", type=float, default=0.15)
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", action="store_true", help="exit 1 on any WARN")
    a = ap.parse_args(argv)
    rows = run(pd.read_csv(a.scores, sep="\t", low_memory=False),
               pd.read_csv(a.thresholds, sep="\t"),
               pd.read_csv(a.cds, sep="\t") if a.cds else None,
               pd.read_csv(a.catalogue, sep="\t", low_memory=False) if a.catalogue else None,
               a.fraction_tol)
    n_warn = report(rows, a.out)
    return 1 if a.strict and n_warn else 0


if __name__ == "__main__":
    raise SystemExit(main())
