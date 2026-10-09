"""hybrid_gtf.py: novel classes split between the RiboCode hybrid and the PRICE/catalogue
full hybrid. Prints PASS/FAIL, exits 1 on failure."""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "workflow", "scripts", "hybrid_gtf.py")
failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


def row(tid, gid, feat="exon", code=None, s=100, e=200):
    extra = f' class_code "{code}";' if code else ""
    return f'chr1\tx\t{feat}\t{s}\t{e}\t.\t+\t.\tgene_id "{gid}"; transcript_id "{tid}";{extra}\n'


def tids(path):
    out = set()
    for l in open(path):
        if 'transcript_id "' in l:
            out.add(l.split('transcript_id "')[1].split('"')[0])
    return out


with tempfile.TemporaryDirectory() as d:
    ref = os.path.join(d, "ref.gtf")
    open(ref, "w").write(row("R1", "G1", "transcript") + row("R1", "G1"))
    ann = os.path.join(d, "ann.gtf")
    open(ann, "w").write("".join(row(t, "M" + t, f, c) for t, c in
                                 (("NU", "u"), ("NK", "k"), ("NO", "o"), ("NX", "x"))
                                 for f in ("transcript", "exon")))
    out = {k: os.path.join(d, k) for k in ("hyb", "full", "rep")}

    def run(codes, ribo):
        cmd = [sys.executable, SCRIPT, "--reference", ref, "--backbone", ref, "--annotated", ann,
               "--class-codes", codes, "--out-hybrid", out["hyb"], "--out-full-hybrid", out["full"],
               "--report", out["rep"]]
        if ribo:
            cmd += ["--ribocode-class-codes", ribo]
        return subprocess.run(cmd, capture_output=True, text=True)

    r = run("u,k,o", "u,o")
    check("split: exit 0", r.returncode, 0)
    check("RiboCode hybrid: backbone + u,o", tids(out["hyb"]), {"R1", "NU", "NO"})
    check("full hybrid: reference + u,k,o", tids(out["full"]), {"R1", "NU", "NK", "NO"})
    rep = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(out["rep"]).readlines()[1:]}
    check("report k: kept full, not RiboCode", rep["k"], ["1", "1", "0"])
    check("report x: kept nowhere", rep["x"], ["1", "0", "0"])

    r = run("u,k,o", None)
    check("default: RiboCode gets the full set", tids(out["hyb"]), {"R1", "NU", "NK", "NO"})

    r = run("u,o", "u,k")
    check("RiboCode classes must be a subset", (r.returncode != 0, "not in --class-codes" in r.stderr), (True, True))

sys.exit(1 if failed else 0)
