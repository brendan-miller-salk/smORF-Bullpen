#!/usr/bin/env python3
"""Check a sample's plastid P-site fit and keep only read lengths it supports.

plastid's psite never fails: when a read length's start-codon metagene is empty
it silently writes --default for that length. With its default median profile
that happened to every length on every library tested, so the P-site
tracks were built from a guess. psite now runs with --aggregate, and its
metagene file then holds raw 5'-end counts per length at positions -50..+49
around annotated start codons. From those counts, for each length:

  share   its fraction of the sample's start-codon reads
  offset  plastid's fitted offset (must lie within --offset-range)
  frame0  fraction of P-sites (5' end + offset) in frame 0, over codons 3-15
          of the CDS, past the initiation peak. Chance is 1/3.

A length is kept if share >= --min-share, frame0 >= --min-frame0 (on at least
--min-body-reads reads) and its offset is in range. Only kept lengths go into
the offsets file make_wiggle uses. That file has no 'default' row, so reads of
every other length are left out of the tracks.

The sample FAILs if it has fewer than --min-start-reads start-codon reads, if its
most abundant length is not kept, or if the kept lengths hold under
--min-kept-share of its reads. This script never exits non-zero on a FAIL; the
psite_report rule decides, so every sample is reported at once.

start_peak_ratio guards against a TIS library (harringtonine / lactimidomycin) left
as `elongation` in the samplesheet: P-sites per codon on codons -1..+1 of the start
codon over P-sites per codon on codons 11-19. Measured on HMC3: elongation 5.9-7.4,
TIS 17.4-25.8. At or above --tis-peak-max the sample is flagged tis_like, which
psite_report always stops on.
"""
import argparse
import sys

import pandas as pd

BODY = (6, 45)   # P-site positions (nt from the A of the AUG) counted for frame
PEAK = (-3, 6)   # start_peak_ratio: codons -1..+1 ...
PEAK_BODY = (33, 60)   # ... over codons 11-19, past the initiation peak


def read_offsets(path):
    """plastid offsets table -> {length: offset}; the 'default' row is dropped."""
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    df = df[df["length"] != "default"]
    return dict(zip(df["length"].astype(int), df["p_offset"].astype(int)))


def length_table(profiles, offsets):
    """One row per read length: reads, share, offset, frame0, body_reads."""
    rows = []
    for length, off in sorted(offsets.items()):
        col = f"{length}-mers"
        counts = profiles[col].fillna(0) if col in profiles else pd.Series(0.0, index=profiles.index)
        psite = profiles["x"] + off
        body = (psite >= BODY[0]) & (psite < BODY[1])
        in_body = counts[body]
        body_reads = in_body.sum()
        frame0 = in_body[(psite[body] % 3) == 0].sum() / body_reads if body_reads else 0.0
        rows.append({"length": length, "reads": counts.sum(), "offset": off,
                     "frame0": frame0, "body_reads": body_reads})
    t = pd.DataFrame(rows)
    total = t["reads"].sum()
    t["share"] = t["reads"] / total if total else 0.0
    return t


def start_peak_ratio(profiles, offsets, offset_range):
    """P-sites per codon at the start codon over per codon in the early CDS body,
    pooled over read lengths with an offset in range."""
    lo, hi = offset_range
    peak = body = 0.0
    for length, off in offsets.items():
        col = f"{length}-mers"
        if col not in profiles or not lo <= off <= hi:
            continue
        counts, psite = profiles[col].fillna(0), profiles["x"] + off
        peak += counts[(psite >= PEAK[0]) & (psite < PEAK[1])].sum()
        body += counts[(psite >= PEAK_BODY[0]) & (psite < PEAK_BODY[1])].sum()
    per_codon = lambda n, w: n / ((w[1] - w[0]) / 3)
    return per_codon(peak, PEAK) / per_codon(body, PEAK_BODY) if body else float("nan")


def judge(t, min_share, min_frame0, min_body_reads, offset_range):
    """Add keep/reason columns: why each length is kept or dropped."""
    lo, hi = offset_range
    reasons = []
    for r in t.itertuples():
        why = []
        if r.share < min_share:
            why.append(f"share<{min_share}")
        if not lo <= r.offset <= hi:
            why.append(f"offset {r.offset} outside {lo}-{hi}")
        if r.body_reads < min_body_reads:
            why.append(f"body_reads<{min_body_reads}")
        elif r.frame0 < min_frame0:
            why.append(f"frame0<{min_frame0}")
        reasons.append("; ".join(why))
    t = t.assign(reason=reasons)
    t["keep"] = t["reason"] == ""
    return t


def summarize(sample, t, min_start_reads, min_kept_share, peak_ratio, tis_peak_max):
    """One report row for the sample, with status PASS/FAIL and why."""
    total = int(t["reads"].sum())
    kept = t[t["keep"]]
    kept_share = kept["share"].sum() if total else 0.0
    dom = t.loc[t["reads"].idxmax()] if total else None
    why = []
    if total < min_start_reads:
        why.append(f"{total} start-codon reads < {min_start_reads}")
    if kept.empty:
        why.append("no read length passed")
    elif dom is not None and not dom["keep"]:
        why.append(f"dominant length {int(dom['length'])} not kept ({dom['reason']})")
    if not kept.empty and kept_share < min_kept_share:
        why.append(f"kept lengths hold {kept_share:.2f} of reads < {min_kept_share}")
    tis_like = bool(peak_ratio >= tis_peak_max)
    if tis_like:
        why.append(f"start-codon peak {peak_ratio:.1f}x the CDS body >= {tis_peak_max}: looks like a TIS library")
    return {
        "sample": sample,
        "status": "FAIL" if why else "PASS",
        "start_codon_reads": total,
        "dominant_length": int(dom["length"]) if dom is not None else "",
        "kept_lengths": ",".join(str(int(x)) for x in kept["length"]),
        "kept_offsets": ",".join(str(int(x)) for x in kept["offset"]),
        "kept_share": round(kept_share, 4),
        "kept_frame0": round((kept["frame0"] * kept["body_reads"]).sum()
                             / kept["body_reads"].sum(), 4) if not kept.empty else "",
        "start_peak_ratio": round(peak_ratio, 2),
        "tis_like": tis_like,
        "reason": "; ".join(why),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--offsets", required=True, help="plastid psite *_p_offsets.txt")
    ap.add_argument("--profiles", required=True, help="plastid psite --aggregate *_metagene_profiles.txt")
    ap.add_argument("--min-share", type=float, default=0.015)
    ap.add_argument("--min-frame0", type=float, default=0.45)
    ap.add_argument("--min-body-reads", type=int, default=100)
    ap.add_argument("--offset-range", type=int, nargs=2, default=(10, 16))
    ap.add_argument("--min-start-reads", type=int, default=20000)
    ap.add_argument("--min-kept-share", type=float, default=0.5)
    ap.add_argument("--tis-peak-max", type=float, default=12)
    ap.add_argument("--out-offsets", required=True, help="kept lengths only, plastid format")
    ap.add_argument("--out-lengths", required=True, help="per-length table")
    ap.add_argument("--out", required=True, help="one-row sample summary")
    a = ap.parse_args()

    profiles = pd.read_csv(a.profiles, sep="\t", comment="#")
    offsets = read_offsets(a.offsets)
    t = judge(length_table(profiles, offsets), a.min_share, a.min_frame0, a.min_body_reads, a.offset_range)
    s = summarize(a.sample, t, a.min_start_reads, a.min_kept_share,
                  start_peak_ratio(profiles, offsets, a.offset_range), a.tis_peak_max)

    with open(a.out_offsets, "w") as fh:
        fh.write("length\tp_offset\n")
        for r in t[t["keep"]].itertuples():
            fh.write(f"{r.length}\t{r.offset}\n")
    t.round(4).to_csv(a.out_lengths, sep="\t", index=False,
                      columns=["length", "reads", "share", "offset", "frame0", "body_reads", "keep", "reason"])
    pd.DataFrame([s]).to_csv(a.out, sep="\t", index=False)
    print(f"psite_check: {a.sample}: {s['status']} kept {s['kept_lengths'] or 'none'} "
          f"({s['kept_share']:.0%} of {s['start_codon_reads']:,} start-codon reads), "
          f"start peak {s['start_peak_ratio']}x {s['reason']}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
