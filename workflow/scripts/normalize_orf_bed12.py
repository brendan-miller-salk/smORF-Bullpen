#!/usr/bin/env python3
"""Normalise an ORF catalogue BED12 before scoring.

Input is the cross-caller ORF catalogue (BED12), holding the union of every
caller's calls. Two things have to be settled before the translation-signature
metrics mean anything:

1. Stop-codon convention. The Drop-off score is defined *across* the stop
   codon: in-frame P-sites in the 15 nt before it versus the 15 nt after. If
   the reported ORFs stop at the last sense codon, that window sits one codon
   5' of where it belongs and straddles the stop, and the score is noise.

   PRICE and Rp-Bp are separate tools with separate conventions, and the
   catalogue merges them. So the convention is *measured* against the genome,
   per caller where the catalogue records one, rather than assumed. A
   catalogue that turns out to be internally inconsistent is reported loudly,
   because scoring it as-is would silently compare unlike things.

2. Duplicate block structures, which callers emit across samples and which
   would otherwise be scored repeatedly.

The catalogue itself is never modified; this writes a derived BED12.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import sys

import numpy as np

STOP_CODONS = {"TAA", "TAG", "TGA"}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def parse_bed12(path: str):
    """Yield (fields, genomic positions in transcript order)."""
    with _open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 12:
                raise ValueError(f"{path}:{lineno}: expected 12+ BED columns, got {len(f)}")
            cs = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offs = [int(x) for x in f[11].rstrip(",").split(",")]
            pos = np.concatenate([np.arange(cs + o, cs + o + n, dtype=np.int64)
                                  for o, n in zip(offs, sizes)])
            pos.sort()
            if f[5] == "-":
                pos = pos[::-1]
            yield f, pos


def positions_to_bed12(chrom, name, strand, pos, score="0", extra=None) -> str:
    asc = np.sort(pos)
    brk = np.flatnonzero(np.diff(asc) != 1)
    starts = np.concatenate(([asc[0]], asc[brk + 1]))
    ends = np.concatenate((asc[brk] + 1, [asc[-1] + 1]))
    cs, ce = int(starts[0]), int(ends[-1])
    sizes = (ends - starts).astype(np.int64)
    offs = (starts - cs).astype(np.int64)
    row = [chrom, str(cs), str(ce), name, score, strand, str(cs), str(ce), "0",
           str(len(sizes)), ",".join(map(str, sizes)) + ",", ",".join(map(str, offs)) + ","]
    if extra:
        row += extra
    return "\t".join(row)


def spliced_seq(genome, chrom, strand, pos) -> str | None:
    try:
        contig = genome[chrom]
    except KeyError:
        return None
    asc = np.sort(pos)
    brk = np.flatnonzero(np.diff(asc) != 1)
    starts = np.concatenate(([asc[0]], asc[brk + 1]))
    ends = np.concatenate((asc[brk] + 1, [asc[-1] + 1]))
    seq = "".join(str(contig[int(s):int(e)]) for s, e in zip(starts, ends)).upper()
    return revcomp(seq) if strand == "-" else seq


def downstream_codon(genome, chrom, strand, pos) -> str | None:
    try:
        contig = genome[chrom]
    except KeyError:
        return None
    last = int(pos[-1])
    if strand == "+":
        seq = str(contig[last + 1:last + 4]).upper()
    else:
        if last - 3 < 0:
            return None
        seq = revcomp(str(contig[last - 3:last]).upper())
    return seq if len(seq) == 3 else None


def extend_3prime(pos: np.ndarray, n: int = 3) -> np.ndarray:
    last = int(pos[-1])
    step = 1 if pos.size < 2 or pos[-1] > pos[0] else -1
    ext = np.arange(last + step, last + step * (n + 1), step, dtype=np.int64)
    return np.concatenate([pos, ext[ext >= 0]])


def harmonise_seqnames(records, genome, stats):
    """Rewrite catalogue seqnames to the genome's own contig names.

    MEASURED on a real cross-caller catalogue: the two
    callers disagree on convention. PRICE emits bare Ensembl names ("1", "X")
    while RiboCode emits UCSC-style ("chr1", "chrX"). Of 35,557 catalogue rows,
    9,433 PRICE rows were bare and 26,080 RiboCode rows were chr-prefixed.

    This fails SILENTLY and expensively. pyfaidx raises KeyError, the codon
    comes back None, and sniff_convention simply sees fewer samples -- so the
    stop-convention check "passes" having never examined one caller. Worse, the
    emitted BED12 keeps the bare names, so the P-site bedGraphs in
    score_orfs.py (keyed on the genome's contig names) never match either, and
    every affected ORF scores zero.

    Unresolvable names are KEPT, not dropped: scaffold ORFs are preserved by
    design elsewhere in this pipeline, so dropping them here would silently
    change the ORF set.
    """
    contigs = set(genome.keys())
    cache: dict[str, str | None] = {}

    def resolve(c: str):
        if c in cache:
            return cache[c]
        cands = [c, f"chr{c}"]
        if c.startswith("chr"):
            cands.append(c[3:])
        hit = next((x for x in cands if x in contigs), None)
        cache[c] = hit
        return hit

    out = []
    for f, pos in records:
        tgt = resolve(f[0])
        if tgt is None:
            stats["seqname_not_in_genome"] += 1
        elif tgt != f[0]:
            stats["seqname_renamed"] += 1
            f = list(f)
            f[0] = tgt
        out.append((f, pos))
    return out


def sniff_convention(records, genome, label: str, sample: int = 3000) -> tuple[str, dict]:
    """Measure whether ORFs end in a stop codon, or are followed by one."""
    incl = excl = checked = 0
    for f, pos in records[:sample]:
        seq = spliced_seq(genome, f[0], f[5], pos)
        if seq is None or len(seq) < 3:
            continue
        checked += 1
        if seq[-3:] in STOP_CODONS:
            incl += 1
        elif downstream_codon(genome, f[0], f[5], pos) in STOP_CODONS:
            excl += 1
    if not checked:
        return "unknown", {}
    frac_incl, frac_excl = incl / checked, excl / checked
    if frac_incl >= 0.9:
        conv = "stop-inclusive"
    elif frac_excl >= 0.9:
        conv = "stop-exclusive"
    else:
        conv = "mixed/unclear"
    print(f"  [{label}] {conv}: of {checked:,} sampled, {incl:,} end in a stop "
          f"({100*frac_incl:.1f}%), {excl:,} are followed by one ({100*frac_excl:.1f}%)",
          file=sys.stderr)
    return conv, {"checked": checked, "incl": incl, "excl": excl}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed12", required=True, help="Cross-caller ORF catalogue BED12.")
    p.add_argument("--genome", required=True, help="Genome FASTA.")
    p.add_argument("--out-bed12", required=True)
    p.add_argument("--caller-column", type=int, default=None,
                   help="1-based column holding the caller name, if the catalogue has one. "
                        "Enables a per-caller convention check.")
    p.add_argument("--extend-stop", choices=("auto", "always", "never"), default="auto",
                   help="auto (default): extend 3 nt only if the catalogue is stop-exclusive.")
    p.add_argument("--no-dedup", dest="dedup", action="store_false")
    p.add_argument("--report", help="Per-ORF TSV: length, start codon, terminal codon.")
    args = p.parse_args(argv)

    try:
        from pyfaidx import Fasta
    except ImportError:
        print("ERROR: needs pyfaidx (conda install -c bioconda pyfaidx).", file=sys.stderr)
        return 1
    genome = Fasta(args.genome, sequence_always_upper=True, rebuild=False)

    print(f"Reading {args.bed12} ...", file=sys.stderr)
    records = list(parse_bed12(args.bed12))
    print(f"  {len(records):,} ORF rows", file=sys.stderr)
    if not records:
        print("ERROR: catalogue is empty.", file=sys.stderr)
        return 1

    # Harmonise seqnames BEFORE sniffing the convention: a catalogue that mixes
    # "1" and "chr1" otherwise sniffs only the half that matches the genome.
    _seq_stats = collections.Counter()
    records = harmonise_seqnames(records, genome, _seq_stats)
    if _seq_stats["seqname_renamed"]:
        print(f"  renamed {_seq_stats['seqname_renamed']:,} seqnames to genome convention",
              file=sys.stderr)
    if _seq_stats["seqname_not_in_genome"]:
        print(f"  WARNING: {_seq_stats['seqname_not_in_genome']:,} ORFs on contigs absent "
              f"from the genome FASTA (kept, but they cannot be codon-checked)",
              file=sys.stderr)

    # --- convention, overall and per caller
    overall, _ = sniff_convention(records, genome, "all")

    per_caller = {}
    if args.caller_column:
        ci = args.caller_column - 1
        groups = collections.defaultdict(list)
        for rec in records:
            if len(rec[0]) > ci:
                groups[rec[0][ci]].append(rec)
        for caller, recs in sorted(groups.items()):
            per_caller[caller], _ = sniff_convention(recs, genome, caller)
        distinct = {c for c in per_caller.values() if c != "unknown"}
        if len(distinct) > 1:
            print("\n  ERROR: callers disagree on the stop-codon convention:", file=sys.stderr)
            for c, v in per_caller.items():
                print(f"    {c}: {v}", file=sys.stderr)
            print("  Scoring a mixed catalogue would compare unlike things: Drop-off would be\n"
                  "  measured across the stop for some ORFs and one codon short for others.\n"
                  "  Normalise each caller's output separately, then merge.", file=sys.stderr)
            return 2

    if args.extend_stop == "always":
        extend = True
    elif args.extend_stop == "never":
        extend = False
    else:
        extend = overall == "stop-exclusive"
        if overall == "mixed/unclear":
            print("\n  ERROR: stop-codon convention is not consistent within the catalogue.\n"
                  "  Neither extending nor leaving it alone is correct for all ORFs.\n"
                  "  Inspect the catalogue, or force with --extend-stop always|never.",
                  file=sys.stderr)
            return 2
    print(f"  extending 3 nt: {extend}", file=sys.stderr)

    stats = collections.Counter()
    start_codons = collections.Counter()
    rows = []
    seen: dict[tuple, str] = {}

    with open(args.out_bed12, "w") as out:
        for f, pos in records:
            chrom, name, strand = f[0], f[3], f[5]
            nt = int(pos.size)
            stats["in_frame_length" if nt % 3 == 0 else "length_not_mult_of_3"] += 1

            if extend:
                pos = extend_3prime(pos)
                stats["extended"] += 1

            line = positions_to_bed12(chrom, name, strand, pos, score=f[4])
            bf = line.split("\t")
            key = (bf[0], bf[1], bf[2], bf[5], bf[10], bf[11])
            if key in seen:
                stats["duplicate_orfs"] += 1
                if args.dedup:
                    stats["collapsed"] += 1
                    continue
            else:
                seen[key] = name

            if args.report:
                seq = spliced_seq(genome, chrom, strand, pos)
                sc = seq[:3] if seq and len(seq) >= 3 else ""
                lc = seq[-3:] if seq and len(seq) >= 3 else ""
                if sc:
                    start_codons[sc] += 1
                if lc in STOP_CODONS:
                    stats["ends_in_stop"] += 1
                rows.append({"orf_id": name, "chrom": chrom, "strand": strand,
                             "length_nt": int(pos.size), "start_codon": sc, "last_codon": lc})

            out.write(line + "\n")
            stats["written"] += 1

    if args.report:
        import pandas as pd
        pd.DataFrame(rows).to_csv(args.report, sep="\t", index=False)

    print(f"\nWrote {stats['written']:,} ORFs to {args.out_bed12}", file=sys.stderr)
    for k in ("in_frame_length", "length_not_mult_of_3", "ends_in_stop",
              "duplicate_orfs", "collapsed", "extended"):
        if stats[k]:
            print(f"  {k:24s} {stats[k]:>10,}", file=sys.stderr)
    if start_codons:
        total = sum(start_codons.values())
        print("\nStart codons (top 8):", file=sys.stderr)
        for c, n in start_codons.most_common(8):
            print(f"  {c}  {n:>9,}  {100*n/total:5.1f}%", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
