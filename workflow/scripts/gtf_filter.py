#!/usr/bin/env python3
"""Keep GTF lines on seqnames in the FASTA that carry a transcript_id.

The step before AGAT: drops gene
rows (no transcript_id), '#' header lines, and features on seqnames the genome
does not have.
"""
import argparse
import re
import statistics
import sys

TX = re.compile(r'transcript_id "[^"]+"')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    names = {line.split("\t", 1)[0] for line in open(a.fai)}
    with open(a.gtf) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    tabs = statistics.median(l.count("\t") for l in lines[:10000])
    if tabs != 8:
        sys.exit(f"gtf_filter: {a.gtf} does not look like a GTF (median {tabs} tabs per line)")
    kept = [l for l in lines if l.split("\t", 1)[0] in names and TX.search(l)]
    with open(a.out, "w") as out:
        out.writelines(kept)
    print(f"gtf_filter: kept {len(kept):,} of {len(lines):,} lines", file=sys.stderr)


if __name__ == "__main__":
    main()
