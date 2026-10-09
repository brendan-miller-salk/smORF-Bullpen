"""psite_check.py: keeps fitted, in-frame read lengths; fails libraries with no usable fit.
Prints PASS/FAIL, exits 1 on failure."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
from psite_check import judge, length_table, summarize  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


X = np.arange(-50, 50)   # plastid --aggregate window around the A of the AUG


def periodic(n, offset):
    """5'-end counts for n reads whose P-sites (5' + offset) sit in frame 0 of the CDS."""
    c = np.zeros(len(X))
    ps = X + offset
    c[(ps >= 0) & (ps % 3 == 0)] = 1.0
    return c * n / c.sum()


def noise(n):
    """5'-end counts spread evenly: no frame preference."""
    return np.full(len(X), n / len(X))


def run(cols, offsets, **kw):
    prof = pd.DataFrame({"x": X, **{f"{k}-mers": v for k, v in cols.items()}})
    t = judge(length_table(prof, offsets), kw.get("min_share", 0.015), kw.get("min_frame0", 0.45),
              kw.get("min_body_reads", 100), kw.get("offset_range", (10, 16)))
    return t, summarize("s", t, kw.get("min_start_reads", 20000), kw.get("min_kept_share", 0.5))


# A good library: 28-30 periodic at 12, 31 periodic at 13, 27 aperiodic, 21 junk offset.
t, s = run({28: periodic(30000, 12), 29: periodic(50000, 12), 30: periodic(15000, 12),
            31: periodic(4000, 13), 27: noise(5000), 21: periodic(3000, 4)},
           {28: 12, 29: 12, 30: 12, 31: 13, 27: 12, 21: 4})
kept = dict(zip(t["length"], t["keep"]))
check("good: status", s["status"], "PASS")
check("good: kept lengths", s["kept_lengths"], "28,29,30,31")
check("good: aperiodic 27 dropped", kept[27], False)
check("good: offset-4 21 dropped", kept[21], False)
check("good: 27 reason", t.set_index("length").loc[27, "reason"], "frame0<0.45")
check("good: in-frame fraction", round(t.set_index("length").loc[29, "frame0"], 3), 1.0)

# The male-muscle shape: few reads, dominated by short aperiodic fragments.
t, s = run({22: noise(1500), 25: periodic(700, 12), 29: noise(80)}, {22: 9, 25: 12, 29: 15})
check("sparse: status", s["status"], "FAIL")
check("sparse: names read floor", "start-codon reads < 20000" in s["reason"], True)
check("sparse: names dominant length", "dominant length 22 not kept" in s["reason"], True)

# Plenty of reads, but the dominant length is aperiodic: fail even though a minor one passes.
t, s = run({29: noise(60000), 28: periodic(20000, 12)}, {29: 12, 28: 12})
check("aperiodic dominant: status", s["status"], "FAIL")
check("aperiodic dominant: 28 still kept", s["kept_lengths"], "28")

# Nothing usable: empty kept list (psite_report stops the run even in warn mode).
t, s = run({29: noise(60000)}, {29: 12})
check("none usable: kept empty", s["kept_lengths"], "")
check("none usable: reason", "no read length passed" in s["reason"], True)

# A length absent from the metagene file counts as zero reads, not a crash.
t, s = run({29: periodic(50000, 12)}, {29: 12, 40: 12})
check("missing column: kept", s["kept_lengths"], "29")

sys.exit(1 if failed else 0)
