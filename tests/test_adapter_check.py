"""adapter_check.py: passes when the adapter is in the reads, fails when it is not.
Prints PASS/FAIL, exits 1 on failure."""
import gzip
import json
import os
import subprocess
import sys
import tempfile

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts",
                      "adapter_check.py")
failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


# Trim Galore 2.3.0 shapes, as seen on real output.
LOG = ("Auto-detecting adapter type...\n"
       "Using Illumina adapter for trimming (count: 90). Second best hit was Nextera (count: 0).\n")
JSON = {"parameters": {"adapters": [{"name": "Illumina", "sequence": "AGATCGGAAGAGC"}]}}
FOOT = "ACGTTGCAACGTTGCAACGTTGCAACGTTG"

with tempfile.TemporaryDirectory() as d:
    js, lg = os.path.join(d, "r.json"), os.path.join(d, "t.log")
    json.dump(JSON, open(js, "w"))
    open(lg, "w").write(LOG)
    for name, tail, want in (("Illumina library", "AGATCGGAAGAGCACACGTCT", 0),
                             ("custom-linker library", "CTGTAGGCACCATCAAT", 1)):
        fq = os.path.join(d, "r.fq.gz")
        with gzip.open(fq, "wt") as fh:
            for i in range(100):
                seq = FOOT + tail + "A" * 20
                fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
        out = os.path.join(d, "o.tsv")
        r = subprocess.run([sys.executable, SCRIPT, "--sample", "s", "--fastq", fq, "--json", js,
                            "--log", lg, "--mode", "auto", "--out", out],
                           capture_output=True, text=True)
        check(f"{name} exit code", r.returncode, want)
        row = open(out).read().splitlines()[1].split("\t")
        check(f"{name} detected_as parsed", row[3], "Illumina")

print("\nALL TESTS PASSED" if not failed else f"\n{failed} FAILED")
sys.exit(1 if failed else 0)
