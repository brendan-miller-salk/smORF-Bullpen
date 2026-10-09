#!/usr/bin/env python3
"""RNA-seq library strandedness for StringTie, from STAR's ReadsPerGene.out.tab.

STAR --quantMode GeneCounts counts each library's uniquely mapped reads on
annotated genes three ways: column 2 unstranded, column 3 read 1 on the gene's
strand (StringTie --fr), column 4 read 1 opposite (--rf, dUTP / TruSeq stranded).
rf_fraction = col4 / (col3 + col4). A sample is
  rf          rf_fraction >= 0.8
  fr          rf_fraction <= 0.2
  unstranded  0.4 <= rf_fraction <= 0.6
  ambiguous   anything else, or no gene reads
With --mode auto every sample must get the same call, and it must not be
ambiguous, or the script exits 1. An explicit --mode is used as given; a
disagreeing call is only a warning. --out-flag gets the one value used.
--require-stranded (antisense 'x' models requested): every sample must call rf or fr,
and the value used must be rf or fr, whatever --mode says; else exit 1.
"""
import argparse
import os
import sys

STRANDED, UNSTRANDED = 0.8, (0.4, 0.6)


def counts(path):
    """(fr, rf) gene read totals, the N_* summary rows skipped."""
    fr = rf = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("N_"):
                continue
            fr += int(f[2])
            rf += int(f[3])
    return fr, rf


def call(fr, rf):
    total = fr + rf
    frac = rf / total if total else 0.0
    if not total:
        return frac, "ambiguous"
    if frac >= STRANDED:
        return frac, "rf"
    if frac <= 1 - STRANDED:
        return frac, "fr"
    if UNSTRANDED[0] <= frac <= UNSTRANDED[1]:
        return frac, "unstranded"
    return frac, "ambiguous"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", nargs="+", required=True,
                    help="<sample>/ReadsPerGene.out.tab; the sample is the parent directory")
    ap.add_argument("--mode", default="auto", choices=["auto", "unstranded", "rf", "fr"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-flag", required=True)
    ap.add_argument("--require-stranded", action="store_true")
    a = ap.parse_args()

    rows = []
    for path in a.counts:
        sample = os.path.basename(os.path.dirname(os.path.abspath(path)))
        fr, rf = counts(path)
        frac, c = call(fr, rf)
        rows.append((sample, fr, rf, frac, c))
    calls = sorted({r[4] for r in rows})
    used = a.mode if a.mode != "auto" else (calls[0] if len(calls) == 1 else None)

    with open(a.out, "w") as fh:
        fh.write("sample\tfr_reads\trf_reads\trf_fraction\tcall\tused\n")
        for s, fr, rf, frac, c in rows:
            fh.write(f"{s}\t{fr}\t{rf}\t{frac:.4f}\t{c}\t{used or 'NA'}\n")
    print(open(a.out).read(), file=sys.stderr, end="")

    if used is None or used == "ambiguous":
        sys.exit(f"rnaseq_strand: cannot infer one strandedness (calls: {', '.join(calls)}; "
                 f"{a.out}). Set rnaseq_strandedness to unstranded, rf or fr, or drop the odd samples.")
    if a.require_stranded and (used not in ("rf", "fr") or set(calls) - {"rf", "fr"}):
        sys.exit(f"rnaseq_strand: novel_class_codes has 'x' (antisense), which needs strand-specific "
                 f"RNA-seq, but the reads say {', '.join(calls)} (using {used}; {a.out}). "
                 f"Drop 'x' or the unstranded samples.")
    if a.mode != "auto" and calls != [a.mode]:
        print(f"WARNING: rnaseq_strandedness is {a.mode} but the reads say {', '.join(calls)}; "
              f"using {a.mode}.", file=sys.stderr)
    with open(a.out_flag, "w") as fh:
        fh.write(used + "\n")


if __name__ == "__main__":
    main()
