"""strand_check.py: forward libraries pass; reverse, unstranded and thin ones fail.
Prints PASS/FAIL, exits 1 on failure."""
import os
import subprocess
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "workflow", "scripts", "strand_check.py")
sys.path.insert(0, os.path.dirname(SCRIPT))
from strand_check import call  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


check("forward", call(950000, 50000, 0.9, 100000)[1:3], ("forward", "PASS"))
check("reverse", call(40000, 960000, 0.9, 100000)[1:3], ("reverse", "FAIL"))
check("unstranded", call(500000, 500000, 0.9, 100000)[1:3], ("ambiguous", "FAIL"))
check("thin forward", call(9500, 500, 0.9, 100000)[1:3], ("forward", "FAIL"))
check("no reads", call(0, 0, 0.9, 100000)[1:3], ("none", "FAIL"))


def summary(path, counts):
    """featureCounts .summary shape: Status row header, one column per BAM path."""
    cols = {f"/x/align/{s}.genome.sorted.bam": [a, 7, 3] for s, a in counts.items()}
    pd.DataFrame(cols, index=["Assigned", "Unassigned_NoFeatures", "Unassigned_Ambiguity"]) \
        .rename_axis("Status").to_csv(path, sep="\t")


with tempfile.TemporaryDirectory() as d:
    s1, s2, out = (os.path.join(d, n) for n in ("s1.summary", "s2.summary", "r.tsv"))
    summary(s1, {"A": 900000, "B": 30000})
    summary(s2, {"A": 20000, "B": 900000})
    r = subprocess.run([sys.executable, SCRIPT, "--sense", s1, "--antisense", s2, "--out", out],
                       capture_output=True, text=True)
    rep = pd.read_csv(out, sep="\t").set_index("sample")
    check("report: samples named from BAM", sorted(rep.index), ["A", "B"])
    check("report: A forward", rep.loc["A", "call"], "forward")
    check("report: B reverse", rep.loc["B", "call"], "reverse")
    check("warn mode exits 0", r.returncode, 0)
    r = subprocess.run([sys.executable, SCRIPT, "--sense", s1, "--antisense", s2, "--out", out,
                        "--strict"], capture_output=True, text=True)
    check("strict exits 1", r.returncode, 1)
    check("strict names B only", "B" in r.stderr.splitlines()[-1] and ": A" not in r.stderr, True)

sys.exit(1 if failed else 0)
