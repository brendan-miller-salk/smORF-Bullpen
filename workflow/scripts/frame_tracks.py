#!/usr/bin/env python
"""ORF-relative P-site frame tracks for one run, as UCSC-named bedGraphs.

Pools the run's per-library plastid P-site bedGraphs, scales them to P-sites per
million in-ORF P-sites, and splits each nucleotide by its frame relative to the
ORF it lies in:

  tier 1  annotated CDS of the canonical transcripts (the set the thresholds are
          fitted on). Frame honours the first CDS phase.
  tier 2  the ORF candidates (scoring/orf_candidates.bed12, stop included), only
          where no tier-1 CDS covers the nucleotide.

Frame is counted from the start codon along the spliced ORF, so frame 0 is in-frame.
A nucleotide two ORFs of one tier claim in different frames is ambiguous and dropped,
as are P-sites outside every ORF; both are reported in the summary.

The denominator is P-sites in ORFs (any frame, both strands), not all P-sites:
mitochondrial rRNA survives rRNA depletion (47% of a heart library's + strand
P-sites), so an all-P-site denominator would track each tissue's mito content.

Writes <bedgraph_dir>/<label>_{plus,minus}_f{0,1,2}.bedgraph with UCSC chromosome
names (contigs with no chromAlias entry are dropped) and --summary. Hub files are
not written here: whoever hosts the bigWigs writes the hub (colors: frame 0 blue
0,0,255, 1 red 255,0,0, 2 gold 230,170,0; one stacked multiWig row per strand).
"""
import argparse
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

UNSET, AMBIG = -1, 3
STRANDS = {"+": "plus", "-": "minus"}
CLASSES = ("outside_orf", "ambiguous", "f0", "f1", "f2")


def read_sizes(fai):
    s = pd.read_csv(fai, sep="\t", header=None, usecols=[0, 1], names=["chrom", "size"], dtype={"chrom": str})
    return dict(zip(s.chrom, s["size"]))


def read_alias(path):
    """Any name in a UCSC chromAlias row -> its UCSC name (first column)."""
    alias = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            for name in f:
                if name:
                    alias[name] = f[0]
    return alias


def cds_orfs(gtf):
    """(chrom, strand, [(start0, end)], phase of the first codon) per transcript's CDS rows."""
    tx, info = defaultdict(list), {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "CDS":
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tx[tid].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in ("0", "1", "2") else 0))
            info[tid] = (f[0], f[6])
    for tid, blocks in tx.items():
        chrom, strand = info[tid]
        blocks.sort()
        first = blocks[0] if strand == "+" else blocks[-1]
        yield chrom, strand, [(s, e) for s, e, _ in blocks], first[2]


def bed12_orfs(bed):
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            starts = [int(x) for x in f[11].rstrip(",").split(",")]
            yield f[0], f[5], [(start + s, start + s + z) for s, z in zip(starts, sizes)], 0


def paint(fmap, blocks, strand, phase):
    """Write each ORF nucleotide's frame into fmap; a disagreement marks it ambiguous."""
    pos = np.concatenate([np.arange(s, e) for s, e in blocks])
    if strand == "-":
        pos = pos[::-1]
    fr = ((np.arange(pos.size) - phase) % 3).astype(np.int8)
    cur = fmap[pos]
    fmap[pos] = np.where((cur == UNSET) | (cur == fr), fr, AMBIG)


def frame_map(size, strand, tier1, tier2):
    """Per-nucleotide frame on one chromosome strand: tier 1 wins, tier 2 fills the rest."""
    m1 = np.full(size, UNSET, np.int8)
    m2 = np.full(size, UNSET, np.int8)
    for _, st, blocks, phase in tier1:
        if st == strand:
            paint(m1, blocks, st, phase)
    for _, st, blocks, phase in tier2:
        if st == strand:
            paint(m2, blocks, st, phase)
    return np.where(m1 != UNSET, m1, m2)


def read_bedgraph(path):
    """plastid bedGraph -> (chrom, 0-based position, value), one row per nucleotide.
    make_wiggle opens with '##' comments and a 'track' line, and merges adjacent
    equal counts into one interval."""
    with open(path) as fh:
        skip = 0
        for line in fh:
            if not (line.startswith("#") or line.startswith("track")):
                break
            skip += 1
    df = pd.read_csv(path, sep="\t", header=None, skiprows=skip, names=["chrom", "start", "end", "val"],
                     dtype={"chrom": str, "start": np.int64, "end": np.int64, "val": np.float64})
    n = (df.end - df.start).to_numpy()
    df = df.loc[df.index.repeat(n)]
    df["start"] = df.start.to_numpy() + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
    return df[["chrom", "start", "val"]]


def pool(paths):
    """{chrom: Series of summed values indexed by position}."""
    df = pd.concat([read_bedgraph(p) for p in paths]).groupby(["chrom", "start"], sort=False)["val"].sum()
    return {c: s.droplevel("chrom") for c, s in df.groupby(level="chrom", sort=False)}


def assign(pooled, sizes, orfs):
    """Split pooled P-sites by frame.
    pooled: {strand: {chrom: Series}}; orfs: {chrom: (tier1, tier2)}.
    Returns ({(strand, frame): [(chrom, pos, raw value)]}, {(strand, class): P-sites})."""
    rows, tally = defaultdict(list), defaultdict(float)
    for strand in STRANDS:
        for chrom, ps in pooled[strand].items():
            tier1, tier2 = orfs.get(chrom, ([], []))
            fr = frame_map(sizes[chrom], strand, tier1, tier2)[ps.index.to_numpy()]
            pos, val = ps.index.to_numpy(), ps.to_numpy()
            for cls, mask in zip(CLASSES, (fr == UNSET, fr == AMBIG, fr == 0, fr == 1, fr == 2)):
                tally[(strand, cls)] += val[mask].sum()
            for f in range(3):
                keep = fr == f
                if keep.any():
                    rows[(strand, f)].append((chrom, pos[keep], val[keep]))
    return rows, tally


def summary_table(label, tally):
    out = []
    for strand, word in STRANDS.items():
        tot = sum(tally[(strand, c)] for c in CLASSES)
        inorf = sum(tally[(strand, f"f{f}")] for f in range(3))
        out.append({"label": label, "strand": word, "psites": round(tot),
                    "share_outside_orf": tally[(strand, "outside_orf")] / tot if tot else float("nan"),
                    "share_ambiguous": tally[(strand, "ambiguous")] / tot if tot else float("nan"),
                    **{f"f{f}_of_in_orf": tally[(strand, f"f{f}")] / inorf if inorf else float("nan")
                       for f in range(3)}})
    return pd.DataFrame(out)


def write_bedgraphs(rows, alias, scale, label, outdir):
    """UCSC-named, LC_ALL=C-sorted bedGraphs, as bedGraphToBigWig wants. Returns dropped contigs."""
    os.makedirs(outdir, exist_ok=True)
    dropped = set()
    for strand, word in STRANDS.items():
        for f in range(3):
            named = []
            for chrom, pos, val in rows[(strand, f)]:
                if chrom in alias:
                    named.append((alias[chrom], pos, val))
                else:
                    dropped.add(chrom)
            with open(os.path.join(outdir, f"{label}_{word}_f{f}.bedgraph"), "w") as out:
                for ucsc, pos, val in sorted(named, key=lambda r: r[0]):
                    o = np.argsort(pos)
                    pd.DataFrame({"c": ucsc, "s": pos[o], "e": pos[o] + 1, "v": np.round(val[o] * scale, 4)}
                                 ).to_csv(out, sep="\t", header=False, index=False)
    return dropped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fw", nargs="+", required=True, help="per-library + strand P-site bedGraphs")
    ap.add_argument("--rc", nargs="+", required=True, help="per-library - strand P-site bedGraphs")
    ap.add_argument("--canonical-gtf", required=True)
    ap.add_argument("--orfs", required=True, help="scoring/orf_candidates.bed12")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--alias", required=True, help="UCSC <genome>.chromAlias.txt")
    ap.add_argument("--label", required=True)
    ap.add_argument("--bedgraph-dir", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()

    orfs = defaultdict(lambda: ([], []))
    for o in cds_orfs(a.canonical_gtf):
        orfs[o[0]][0].append(o)
    for o in bed12_orfs(a.orfs):
        orfs[o[0]][1].append(o)

    pooled = {"+": pool(a.fw), "-": pool(a.rc)}
    rows, tally = assign(pooled, read_sizes(a.fai), orfs)
    in_orf = sum(tally[(st, f"f{f}")] for st in STRANDS for f in range(3))
    if not in_orf:
        raise SystemExit("no P-sites fall in any ORF; check the bedGraphs and ORF files share chromosome names")
    total = sum(tally.values())
    print(f"{a.label}: {total:,.0f} pooled P-sites, {in_orf:,.0f} in ORFs (the denominator)")

    dropped = write_bedgraphs(rows, read_alias(a.alias), 1e6 / in_orf, a.label, a.bedgraph_dir)
    summ = summary_table(a.label, tally)
    summ.to_csv(a.summary, sep="\t", index=False, float_format="%.4f")
    print(summ.to_string(index=False, float_format="%.3f"))
    if dropped:
        print(f"no UCSC alias, left out of the tracks: {' '.join(sorted(dropped))}")


if __name__ == "__main__":
    main()
