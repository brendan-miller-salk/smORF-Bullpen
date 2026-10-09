#!/usr/bin/env python3
"""fastp JSON reports -> one row per lane: the adapters used (detected or given) and
what was trimmed. Warns when lanes used different adapters; never fails."""
import argparse
import json
import os
import sys

COLS = ["lane", "read1_adapter", "read2_adapter", "reads_before", "reads_after",
        "adapter_trimmed_reads", "adapter_trimmed_pct"]


def row(path):
    j = json.load(open(path))
    a, s = j.get("adapter_cutting", {}), j["summary"]
    before, cut = s["before_filtering"]["total_reads"], a.get("adapter_trimmed_reads", 0)
    return [os.path.basename(path).removesuffix(".fastp.json"),
            a.get("read1_adapter_sequence", "none"), a.get("read2_adapter_sequence", "."),
            before, s["after_filtering"]["total_reads"], cut, f"{100 * cut / before:.2f}" if before else "0"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = [row(p) for p in a.json]
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    print(open(a.out).read(), file=sys.stderr, end="")
    used = {(r[1], r[2]) for r in rows}
    if len(used) > 1:
        print(f"WARNING: lanes used {len(used)} different adapter sets (a lane with few adapter "
              f"reads can come back unspecified); see {a.out}.", file=sys.stderr)


if __name__ == "__main__":
    main()
