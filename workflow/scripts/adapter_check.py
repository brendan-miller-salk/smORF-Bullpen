#!/usr/bin/env python3
"""Confirm the adapter Trim Galore used is actually in this library.

Trim Galore's auto-detection only knows Illumina, Nextera and small RNA. On a
library with a custom linker it still picks one of them and under-trims
without failing. So, whatever the mode, this counts how many of the first
--n raw reads contain the adapter's first 13 bases (the length of Trim
Galore's own Illumina probe), and fails the sample below --min-fraction.
A Ribo-seq read is far longer than its ~30 nt footprint, so nearly every read
should carry the adapter; a wrong one gives close to 0%.
"""
import argparse
import gzip
import json
import re
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fastq", required=True, help="raw input FASTQ")
    ap.add_argument("--json", required=True, help="Trim Galore JSON report")
    ap.add_argument("--log", required=True, help="Trim Galore stdout/stderr")
    ap.add_argument("--mode", required=True, help="'auto' or the configured sequence")
    ap.add_argument("--n", type=int, default=1_000_000)
    ap.add_argument("--min-fraction", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    adapters = json.load(open(a.json))["parameters"]["adapters"]
    if not adapters:
        sys.exit(f"adapter_check: {a.sample}: Trim Galore used no adapter")
    adapter = adapters[0] if isinstance(adapters[0], str) else adapters[0].get("sequence", "")
    log = open(a.log).read()
    m = re.search(r"Using (\S+) adapter for trimming \(count: (\d+)\)\. "
                  r"Second best hit was (\S+) \(count: (\d+)\)", log)
    label, count, second, second_count = m.groups() if m else ("user-specified", "", "", "")

    probe = adapter[:13]
    seen = hit = 0
    with gzip.open(a.fastq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seen += 1
                hit += probe in line
                if seen >= a.n:
                    break
    frac = hit / seen if seen else 0.0

    with open(a.out, "w") as out:
        out.write("sample\tmode\tadapter\tdetected_as\tdetect_count\tsecond_best\t"
                  "second_best_count\treads_checked\tfraction_with_adapter\n")
        out.write(f"{a.sample}\t{'auto' if a.mode == 'auto' else 'explicit'}\t{adapter}\t{label}\t"
                  f"{count}\t{second}\t{second_count}\t{seen}\t{frac:.4f}\n")
    print(f"adapter_check: {a.sample}: {adapter} ({label}) in {frac:.1%} of {seen:,} reads",
          file=sys.stderr)
    if frac < a.min_fraction:
        sys.exit(f"adapter_check: {a.sample}: only {frac:.1%} of reads contain {probe}, below "
                 f"adapter_min_fraction {a.min_fraction}. The library likely uses a custom "
                 f"linker Trim Galore cannot detect; set `adapter:` to its sequence.")


if __name__ == "__main__":
    main()
