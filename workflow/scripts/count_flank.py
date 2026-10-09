#!/usr/bin/env python3
"""Fit how far canonical-CDS reads spill past the CDS ends: the counting flank.

featureCounts assigns a read to any feature it overlaps by one base. Reads from
ribosomes on an annotated start or stop codon overhang the CDS, so a novel smORF
segment just outside a canonical CDS would collect them. build_count_gtf.py
therefore removes, from every smORF, the canonical CDS widened by a 5' and a 3'
flank. This script measures those two flanks from the data.

  windows  From one canonical transcript per gene, keep coding transcripts with
           >= --utr nt of 5' and 3' UTR, >= 2 * --cds nt of CDS, off chrM, and
           no other gene's CDS on the same strand within the windows. Write each
           transcript's two windows as transcript-ordered genomic positions:
             start side: d = 1..utr nt upstream of the first CDS nt, d <= 0 inside
             stop side:  d = 1..utr nt downstream of the last CDS nt, d <= 0 inside
           (GTF CDS excludes the stop codon, so the stop is d = 1..3: distances are
           measured from the CDS features featureCounts sees.) Also a BED of the
           windows per strand, for samtools depth.

  fit      Whole-read sense coverage of uniquely mapped reads at those positions
           (samtools depth, one file per sample and strand), summed over samples.
           Each transcript is divided by its mean CDS coverage (d in -cds..-10),
           dropping those with mean < --min-cds-cov reads, and a --trim trimmed
           mean across transcripts gives a metagene m(d) per side. Then
             cds_level = mean m(d) over the CDS part of the window
             baseline  = median m(d) over d in --baseline: the UTR just past a
                         CDS footprint's reach. Not the distal UTR, whose level
                         can differ (5' UTRs slope), which would put the band out of reach
             flank     = the first d >= 0 where a --smooth nt running median of m(d)
                         is within --excess-frac of the way from baseline to cds_level,
                         and never less than --min-flank-5p / --min-flank-3p (how far a
                         footprint physically overhangs a start or stop codon).
"""
import argparse
import re
import sys
from bisect import bisect_left
from collections import defaultdict

import numpy as np
import pandas as pd

ATTR = re.compile(r'(\S+) "([^"]*)"')


def attrs(field):
    return dict(ATTR.findall(field))


def read_gtf(path, features):
    """Rows of the given feature types -> DataFrame(chrom, feature, start, end, strand, gene_id, transcript_id).
    Coordinates stay 1-based inclusive, as in the GTF."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in features:
                continue
            a = attrs(f[8])
            rows.append((f[0], f[2], int(f[3]), int(f[4]), f[6], a.get("gene_id", ""), a.get("transcript_id", "")))
    return pd.DataFrame(rows, columns=["chrom", "feature", "start", "end", "strand", "gene_id", "transcript_id"])


def tx_positions(exons, strand):
    """Exon intervals (1-based inclusive) -> genomic positions in transcript 5'->3' order."""
    pos = np.concatenate([np.arange(s, e + 1) for s, e in sorted(exons)])
    return pos[::-1] if strand == "-" else pos


def transcript_windows(exons, cds, strand, utr, cds_nt):
    """(start_positions, stop_positions) for one transcript, each indexed d = utr..-(cds_nt-1),
    i.e. array[i] is at d = utr - i for the start side and d = i - (cds_nt - 1) for the stop side.
    None if the transcript is too short on either side."""
    pos = tx_positions(exons, strand)
    idx = np.flatnonzero(np.isin(pos, tx_positions(cds, strand)))
    if len(idx) < 2 * cds_nt:
        return None
    first, last = idx[0], idx[-1]
    if first < utr or len(pos) - 1 - last < utr:
        return None
    start = pos[first - utr:first + cds_nt]          # d = utr .. -(cds_nt-1)
    stop = pos[last - cds_nt + 1:last + utr + 1]     # d = -(cds_nt-1) .. utr
    return start, stop


def runs(positions):
    """Genomic positions -> contiguous (first, last) runs, 1-based inclusive."""
    p = np.sort(positions)
    breaks = np.flatnonzero(np.diff(p) > 1)
    starts = np.concatenate([[p[0]], p[breaks + 1]])
    ends = np.concatenate([p[breaks], [p[-1]]])
    return list(zip(starts.tolist(), ends.tolist()))


def cds_index(cds_rows):
    """(chrom, strand) -> (sorted starts, ends, gene ids) of every CDS row, for neighbour checks."""
    idx = {}
    for (chrom, strand), g in cds_rows.groupby(["chrom", "strand"]):
        g = g.sort_values("start")
        idx[(chrom, strand)] = (g["start"].to_numpy(), g["end"].to_numpy(), g["gene_id"].to_numpy())
    return idx


def has_neighbour(idx, chrom, strand, lo, hi, gene, maxlen):
    """True if another gene's CDS overlaps [lo, hi] on this strand. maxlen = longest CDS row."""
    if (chrom, strand) not in idx:
        return False
    s, e, g = idx[(chrom, strand)]
    i = bisect_left(s, lo - maxlen)
    while i < len(s) and s[i] <= hi:
        if e[i] >= lo and g[i] != gene:
            return True
        i += 1
    return False


def cmd_windows(a):
    canon = read_gtf(a.canonical_gtf, {"exon", "CDS"})
    allcds = read_gtf(a.gtf, {"CDS"})
    idx, maxlen = cds_index(allcds), int((allcds["end"] - allcds["start"]).max()) + 1
    rows, beds = [], defaultdict(list)
    kept = 0
    for tid, g in canon.groupby("transcript_id", sort=False):
        chrom, strand, gene = g["chrom"].iat[0], g["strand"].iat[0], g["gene_id"].iat[0]
        if chrom == "chrM":
            continue
        cds = g[g["feature"] == "CDS"]
        if cds.empty:
            continue
        exons = list(zip(g.loc[g["feature"] == "exon", "start"], g.loc[g["feature"] == "exon", "end"]))
        w = transcript_windows(exons, list(zip(cds["start"], cds["end"])), strand, a.utr, a.cds)
        if w is None:
            continue
        start, stop = w
        allpos = np.concatenate([start, stop])
        if has_neighbour(idx, chrom, strand, int(allpos.min()), int(allpos.max()), gene, maxlen):
            continue
        kept += 1
        d_start = a.utr - np.arange(len(start))
        d_stop = np.arange(len(stop)) - (a.cds - 1)
        for side, p, d in (("start", start, d_start), ("stop", stop, d_stop)):
            rows.append(pd.DataFrame({"transcript_id": tid, "side": side, "d": d, "chrom": chrom,
                                      "strand": strand, "pos": p}))
            beds[strand].extend((chrom, lo - 1, hi) for lo, hi in runs(p))
    if not rows:
        sys.exit("count_flank windows: no transcript passed the filters")
    pd.concat(rows).to_csv(a.out_map, sep="\t", index=False)
    for strand, path in (("+", a.out_bed_plus), ("-", a.out_bed_minus)):
        b = pd.DataFrame(beds[strand], columns=["chrom", "start", "end"]).sort_values(["chrom", "start"])
        b.to_csv(path, sep="\t", header=False, index=False)
    print(f"count_flank windows: {kept} transcripts", file=sys.stderr)


def read_depth(paths):
    """samtools depth outputs (chrom, pos, depth) -> read depth summed over libraries, indexed (chrom, pos)."""
    total = None
    for p in paths:
        df = pd.read_csv(p, sep="\t", header=None, dtype={0: str})
        s = pd.Series(df.iloc[:, 2:].sum(axis=1).to_numpy(), index=pd.MultiIndex.from_arrays([df[0], df[1]]))
        s = s.groupby(level=[0, 1]).sum()
        total = s if total is None else total.add(s, fill_value=0)
    return total


def metagene(m, cds_nt, min_cds_cov, trim=0.05):
    """Per-side profile: each transcript over its own CDS mean, then a trimmed mean across
    transcripts at each position (drops the top and bottom `trim`, so one transcript's
    read pile cannot make a block; a median would be 0 wherever most transcripts are).
    Returns (profiles, n)."""
    from scipy.stats import trim_mean
    prof, n = {}, {}
    for side, g in m.groupby("side"):
        wide = g.pivot(index="transcript_id", columns="d", values="cov").fillna(0)
        cds_cols = [d for d in wide.columns if -(cds_nt - 1) <= d <= -10]
        norm = wide[cds_cols].mean(axis=1)
        keep = norm >= min_cds_cov
        normed = wide[keep].div(norm[keep], axis=0)
        prof[side] = pd.Series(trim_mean(normed.to_numpy(), trim, axis=0), index=normed.columns)
        n[side] = int(keep.sum())
    return prof, n


def fit_side(prof, cds_nt, baseline, excess_frac, max_flank, min_flank=0, smooth=9):
    """(flank, cds_level, baseline_level, reached) from one metagene side.

    flank = the first d >= 0 at which a `smooth`-nt running median of the profile is
    within excess_frac of the way from baseline to CDS level, and at least min_flank
    (the footprint's physical overhang). `First` rather than `stays below out to
    max_flank`: a UTR whose own coverage sits near the band (5' UTRs with uORFs) would
    otherwise be cut arbitrarily far. reached is False if it never got there."""
    prof = prof.sort_index()
    cds_level = prof[[d for d in prof.index if -(cds_nt - 1) <= d <= -10]].mean()
    base = prof[[d for d in prof.index if baseline[0] <= d <= baseline[1]]].median()
    limit = excess_frac * (cds_level - base)
    sm = prof.rolling(smooth, center=True, min_periods=1).median()
    out = sm[[d for d in sm.index if 0 <= d <= max_flank]]
    below = out.index[(out - base) < limit]
    reached = len(below) > 0
    flank = int(below.min()) if reached else max_flank
    return max(flank, min_flank), float(cds_level), float(base), reached


def cmd_fit(a):
    m = pd.read_csv(a.map, sep="\t", dtype={"chrom": str})
    cov = []
    for strand, paths in (("+", a.depth_plus), ("-", a.depth_minus)):
        s = read_depth(paths)
        sub = m[m["strand"] == strand]
        key = pd.MultiIndex.from_arrays([sub["chrom"], sub["pos"]])
        cov.append(pd.Series(s.reindex(key).fillna(0).to_numpy(), index=sub.index))
    m["cov"] = pd.concat(cov)
    prof, n = metagene(m, a.cds, a.min_cds_cov, a.trim)
    res = {}
    for side, key in (("start", "5p"), ("stop", "3p")):
        floor = a.min_flank_5p if key == "5p" else a.min_flank_3p
        flank, lvl, base, reached = fit_side(prof[side], a.cds, a.baseline, a.excess_frac,
                                             a.max_flank, floor, a.smooth)
        if not reached:
            print(f"count_flank fit: WARN {key} excess never fell below the limit by {a.max_flank} nt",
                  file=sys.stderr)
        res.update({f"flank_{key}": flank, f"cds_level_{key}": round(lvl, 4),
                    f"baseline_{key}": round(base, 4), f"n_transcripts_{key}": n[side]})
    pd.DataFrame(list(res.items()), columns=["metric", "value"]).to_csv(a.out, sep="\t", index=False)
    pd.DataFrame({"d": prof["start"].index, "start": prof["start"].values}).merge(
        pd.DataFrame({"d": prof["stop"].index, "stop": prof["stop"].values}), on="d", how="outer") \
        .sort_values("d").to_csv(a.out_profile, sep="\t", index=False)
    if a.plot:
        plot(prof, res, a)
    print(f"count_flank fit: 5' flank {res['flank_5p']} nt ({res['n_transcripts_5p']} transcripts), "
          f"3' flank {res['flank_3p']} nt ({res['n_transcripts_3p']})", file=sys.stderr)


def plot(prof, res, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, side, key, title in ((axes[0], "start", "5p", "start codon (5' UTR to the left)"),
                                 (axes[1], "stop", "3p", "stop (3' UTR to the right)")):
        p = prof[side].sort_index()
        x = -p.index if side == "start" else p.index          # 5' UTR drawn to the left
        ax.plot(x, p.values, lw=1, color="#3060a0")
        ax.axhline(res[f"baseline_{key}"], color="#888", lw=0.8, ls=":")
        f = res[f"flank_{key}"]
        ax.axvline(-f if side == "start" else f, color="#c04040", lw=1, label=f"fitted flank {f} nt")
        ax.axvline(-35 if side == "start" else 35, color="#888", lw=1, ls="--", label="35 nt")
        ax.axvline(0, color="#000", lw=0.5)
        ax.set_yscale("log")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("nt from CDS edge (CDS at 0 and inward)")
        ax.legend(fontsize=8, frameon=False)
    axes[0].set_ylabel("coverage / transcript CDS mean")
    fig.tight_layout()
    fig.savefig(a.plot, dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("windows")
    w.add_argument("--canonical-gtf", required=True)
    w.add_argument("--gtf", required=True, help="full annotation, for neighbouring CDS")
    w.add_argument("--utr", type=int, default=300)
    w.add_argument("--cds", type=int, default=60)
    w.add_argument("--out-map", required=True)
    w.add_argument("--out-bed-plus", required=True)
    w.add_argument("--out-bed-minus", required=True)
    f = sub.add_parser("fit")
    f.add_argument("--map", required=True)
    f.add_argument("--depth-plus", nargs="+", required=True)
    f.add_argument("--depth-minus", nargs="+", required=True)
    f.add_argument("--cds", type=int, default=60)
    f.add_argument("--baseline", type=int, nargs=2, default=(40, 100),
                   help="UTR level, nt from the CDS edge: just past a CDS footprint's reach (~25 nt)")
    f.add_argument("--excess-frac", type=float, default=0.05)
    f.add_argument("--max-flank", type=int, default=100)
    f.add_argument("--min-cds-cov", type=float, default=5.0)
    f.add_argument("--trim", type=float, default=0.05, help="trimmed-mean fraction per tail")
    f.add_argument("--smooth", type=int, default=9, help="running-median width, nt")
    f.add_argument("--min-flank-5p", type=int, default=12, help="P-site offset: overhang upstream of a start")
    f.add_argument("--min-flank-3p", type=int, default=15, help="overhang past the last CDS nt at a stop")
    f.add_argument("--out", required=True)
    f.add_argument("--out-profile", required=True)
    f.add_argument("--plot")
    a = ap.parse_args()
    cmd_windows(a) if a.cmd == "windows" else cmd_fit(a)


if __name__ == "__main__":
    main()
