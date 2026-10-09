"""frame_tracks.py: ORF-relative frame per P-site, tier priority, bedGraph expansion,
in-ORF normalization, UCSC renaming. Prints PASS/FAIL, exits 1 on failure."""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "workflow", "scripts"))
from frame_tracks import AMBIG, UNSET, frame_map, paint, read_bedgraph  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


# + strand, spliced: exons [10,16) and [20,26). Frame runs on through the junction.
m = np.full(40, UNSET, np.int8)
paint(m, [(10, 16), (20, 26)], "+", 0)
check("+ first exon", m[10:16].tolist(), [0, 1, 2, 0, 1, 2])
check("+ second exon continues", m[20:26].tolist(), [0, 1, 2, 0, 1, 2])
check("+ intron untouched", m[16:20].tolist(), [UNSET] * 4)

# - strand: the start codon is at the high end.
m = np.full(40, UNSET, np.int8)
paint(m, [(10, 15), (20, 24)], "-", 0)
check("- counts from the right", m[20:24].tolist(), [0, 2, 1, 0])   # 23 is frame 0
check("- second exon continues", m[10:15].tolist(), [2, 1, 0, 2, 1])

# CDS phase 1: the first codon starts one nucleotide in.
m = np.full(10, UNSET, np.int8)
paint(m, [(0, 7)], "+", 1)
check("phase 1 shifts frame 0", m[:7].tolist(), [2, 0, 1, 2, 0, 1, 2])

# Two same-tier ORFs: agreeing frames stay, disagreeing ones go ambiguous.
m = np.full(20, UNSET, np.int8)
paint(m, [(0, 9)], "+", 0)
paint(m, [(3, 12)], "+", 0)      # in frame with the first
check("agreeing overlap keeps frame", m[3:9].tolist(), [0, 1, 2, 0, 1, 2])
paint(m, [(4, 10)], "+", 0)      # shifted by one
check("conflicting overlap is ambiguous", int(m[5]), AMBIG)

# Tier 1 (CDS) wins over tier 2; tier 2 fills where no CDS is; other strand ignored.
cds = [("c", "+", [(0, 6)], 0)]
orfs = [("c", "+", [(1, 10)], 0), ("c", "-", [(12, 15)], 0)]
fm = frame_map(20, "+", cds, orfs)
check("tier 1 wins", fm[1:6].tolist(), [1, 2, 0, 1, 2])
check("tier 2 fills after the CDS", fm[6:10].tolist(), [2, 0, 1, 2])
check("- strand ORF absent from + map", fm[12:15].tolist(), [UNSET] * 3)

with tempfile.TemporaryDirectory() as d:
    def w(name, text):
        p = os.path.join(d, name)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    # plastid header, and a merged 3-nt interval
    head = "## date = 'x'\n## args = {\n##   }\ntrack type=bedGraph name=x color=0,0,0\n"
    fw1 = w("a_fw.bedgraph", head + "chrZ\t0\t1\t4\nchrZ\t3\t6\t2\nctg1.1\t0\t1\t5\n")
    bg = read_bedgraph(fw1)
    check("header skipped, merged run expanded", bg.start.tolist(), [0, 3, 4, 5, 0])
    check("expanded values repeat", bg.val.tolist(), [4.0, 2.0, 2.0, 2.0, 5.0])

    # End to end. CDS chrZ:1-6 (GTF) = [0,6) +; ORF candidates chrZ [9,15) + and ctg1.1 [0,3) +
    # (no UCSC alias); minus strand has one outside-ORF P-site. In ORFs: f0 4 (pos 0), f0 2 (3),
    # f1 2 (4), f2 2 (5), f1 3 (lib 2, pos 10), f0 5 (ctg1.1): 18.
    fw2 = w("b_fw.bedgraph", head + "chrZ\t10\t11\t3\n")
    rc1 = w("a_rc.bedgraph", head + "chrZ\t30\t31\t7\n")
    rc2 = w("b_rc.bedgraph", head)
    gtf = w("c.gtf", 'chrZ\tx\tCDS\t1\t6\t.\t+\t0\ttranscript_id "T1";\n')
    bed = w("o.bed12", "chrZ\t9\t15\torf1\t0\t+\t9\t15\t0\t1\t6,\t0,\n"
                       "ctg1.1\t0\t3\torf2\t0\t+\t0\t3\t0\t1\t3,\t0,\n")
    fai = w("g.fai", "chrZ\t40\t0\t0\t0\nctg1.1\t10\t0\t0\t0\n")
    alias = w("a.txt", "# ucsc\tensembl\nchrZ\tZ\n")
    out = os.path.join(d, "bg")
    r = subprocess.run([sys.executable, os.path.join(HERE, "..", "workflow", "scripts", "frame_tracks.py"),
                        "--fw", fw1, fw2, "--rc", rc1, rc2, "--canonical-gtf", gtf, "--orfs", bed,
                        "--fai", fai, "--alias", alias, "--label", "t", "--bedgraph-dir", out,
                        "--summary", os.path.join(d, "s.tsv")],
                       capture_output=True, text=True)
    check("script exit 0", r.returncode, 0)
    if r.returncode:
        print(r.stderr)
    else:
        def load(n):
            p = os.path.join(out, n)
            return pd.read_csv(p, sep="\t", header=None) if os.path.getsize(p) else pd.DataFrame()

        f0, f1 = load("t_plus_f0.bedgraph"), load("t_plus_f1.bedgraph")
        # denominator: all 18 in-ORF P-sites, including the dropped contig's
        check("f0 positions", f0[1].tolist(), [0, 3])
        check("f0 scaled per million in-ORF", f0[3].round(1).tolist(), [round(4e6 / 18, 1), round(2e6 / 18, 1)])
        check("f1 includes the tier-2 ORF", f1[1].tolist(), [4, 10])
        check("unaliased contig dropped", "ctg1.1" in r.stdout, True)
        check("minus strand has no in-ORF P-sites", len(load("t_minus_f0.bedgraph")), 0)
        s = pd.read_csv(os.path.join(d, "s.tsv"), sep="\t").set_index("strand")
        check("plus: every P-site in an ORF", s.loc["plus", "share_outside_orf"], 0.0)
        check("minus all outside ORFs", s.loc["minus", "share_outside_orf"], 1.0)

sys.exit(1 if failed else 0)
