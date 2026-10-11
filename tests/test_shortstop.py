"""ShortStop input and gate checks. Prints PASS/FAIL, exits 1 on failure.

Pins the GTF shape ShortStop needs to assemble spliced minus-strand ORFs
5'-first: CDS rows ascending by coordinate, each with a distinct attribute
string. Identical attributes on every row make ShortStop re-sort them and
translate the exons 3'-first (seen on 1,672 real ORFs).
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "workflow", "scripts", "shortstop")
failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


with tempfile.TemporaryDirectory() as d:
    # Minus-strand two-exon ORF: protein M K L V * read 5'->3' on the minus strand.
    # Transcript-order CDS: ATG AAA CT | G GTT TAA  (exon 1 = high coords).
    exon1_tx, exon2_tx = "ATGAAACT", "GGTTTAA"
    pad = "C" * 60
    intron = "C" * 20
    # Genome (plus strand), ascending: pad, rc(exon2), intron, rc(exon1), pad
    genome = pad + rc(exon2_tx) + intron + rc(exon1_tx) + pad
    fa = os.path.join(d, "g.fa")
    with open(fa, "w") as fh:
        fh.write(">chrT\n" + genome + "\n")
    e2s = len(pad)
    e1s = e2s + len(exon2_tx) + len(intron)
    bed = os.path.join(d, "o.bed12")
    start, end = e2s, e1s + len(exon1_tx)
    with open(bed, "w") as fh:
        fh.write("\t".join(map(str, ["chrT", start, end, "orf_test", 0, "-", start, end, 0, 2,
                                     f"{len(exon2_tx)},{len(exon1_tx)},",
                                     f"0,{e1s - start},"])) + "\n")
    out = {k: os.path.join(d, k) for k in ("gtf", "skip", "exp")}
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "shortstop_input.py"),
                        "--bed12", bed, "--genome", fa, "--min-aa", "1",
                        "--out-gtf", out["gtf"], "--out-skipped", out["skip"],
                        "--out-expected", out["exp"]], capture_output=True, text=True)
    check("shortstop_input exits 0", r.returncode, 0)
    rows = [l.rstrip("\n").split("\t") for l in open(out["gtf"])]
    cds = [x for x in rows if x[2] == "CDS"]
    check("two CDS rows", len(cds), 2)
    starts = [int(x[3]) for x in cds]
    check("CDS rows ascending by coordinate", starts, sorted(starts))
    check("CDS attribute strings distinct", len({x[8] for x in cds}), len(cds))
    check("gene_id unquoted orf_id first", cds[0][8].split(";")[0], "gene_id orf_test")
    exp = open(out["exp"]).read().splitlines()[1].split("\t")
    check("expected protein, stop trimmed", exp[1], "MKLV")

    # Gate: a mismatching translation must fail shortstop_table.py.
    pred = os.path.join(d, "pred.csv")
    seqs = os.path.join(d, "seqs.csv")
    cons = os.path.join(d, "cons.csv")
    with open(pred, "w") as fh:
        fh.write("orf_id,classification,probability\norf_test,prisms,0.9\n")
    with open(cons, "w") as fh:
        fh.write("orf_id,sam_consensus_percentage,sam_consensus_50,sam_consensus_90,sam_consensus_100\n"
                 "orf_test,92.0,True,True,False\n")
    for aa, want in (("MKLV", 0), ("MVLK", 1)):
        with open(seqs, "w") as fh:
            fh.write(f"orf_id,aa_seq\norf_test,{aa}\n")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "shortstop_table.py"),
                            "--orfs", bed, "--predictions", pred, "--sequences", seqs, "--consensus", cons,
                            "--expected", out["exp"], "--skipped", out["skip"],
                            "--out", os.path.join(d, "t.tsv"), "--report", os.path.join(d, "r.tsv")],
                           capture_output=True, text=True)
        check(f"gate exit code when ShortStop translated {aa}", r.returncode, want)
        if want == 0:
            t = open(os.path.join(d, "t.tsv")).read().splitlines()
            row = dict(zip(t[0].split("\t"), t[1].split("\t")))
            check("highest consensus level", row["shortstop_sam_consensus"], "SAM-Consensus-90")
            check("consensus percentage", row["shortstop_sam_consensus_percentage"], "92.0")

print("\nALL TESTS PASSED" if not failed else f"\n{failed} FAILED")
sys.exit(1 if failed else 0)
