"""reviewed / unreviewed smORF GTF + FASTA checks. Prints PASS/FAIL, exits 1 on failure.

Synthetic 60-nt chr1. orf_plus: one block, + strand, ATG AAA TAG at 0-based 3..12.
orf_minus: two blocks on - strand, 20..24 and 30..37; read 3'->5' on the plus strand
it is ATG GCC TGA after reverse-complementing, with the stop split across blocks.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "workflow", "scripts", "smorf_sets.py")

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# minus ORF, 5'->3' on its strand: ATGGCCTGA (9 nt) = blocks 30..37 (7 nt, rc) + 20..24 (2 nt)?
# Built directly: take the ORF, reverse-complement, and lay it on plus positions
# [20,22) + [30,37): 2 + 7 = 9 nt.
MINUS = "ATGGCCTGA"
plus_view = rc(MINUS)                       # as it reads on the plus strand
seq = ["N"] * 60
seq[3:12] = "ATGAAATAG"
for p, c in zip(list(range(20, 22)) + list(range(30, 37)), plus_view):
    seq[p] = c
GENOME = "".join(seq)

BED = [
    "chr1\t3\t12\torf_plus\t0\t+\t3\t12\t0\t1\t9,\t0,",
    "chr1\t20\t37\torf_minus\t0\t-\t20\t37\t0\t2\t2,7,\t0,10,",
    "chr1\t3\t12\torf_sp\t0\t+\t3\t12\t0\t1\t9,\t0,",
    "chr1\t3\t12\torf_spid\t0\t+\t3\t12\t0\t1\t9,\t0,",
    "chr1\t3\t12\torf_long\t0\t+\t3\t12\t0\t1\t9,\t0,",
    "chr1\t3\t12\torf_excl\t0\t+\t3\t12\t0\t1\t9,\t0,",
]
ORFS = [  # orf_id, orf_set, is_smorf, blastoff_class, aa_seq, smorf_type, host
    ("orf_plus", "primary", "True", "novel", "MK", "uORF", "Gene1"),
    ("orf_minus", "comprehensive", "True", "tr_identical", "MA", "TrEMBL", ""),
    ("orf_sp", "comprehensive", "True", "sp_partial", "MK", "C-Iso", "Gene2"),
    ("orf_spid", "primary", "True", "sp_identical", "MK", "annotated_CDS", "Gene2"),
    ("orf_long", "primary", "False", "novel", "MK", "lncRNA", ""),
    ("orf_excl", "excluded", "True", "novel", "MK", "uORF", ""),
]

with tempfile.TemporaryDirectory() as d:
    j = lambda x: os.path.join(d, x)  # noqa: E731
    open(j("g.fa"), "w").write(">chr1\n" + GENOME + "\n")
    open(j("o.bed12"), "w").write("\n".join(BED) + "\n")
    with open(j("orfs.csv"), "w") as fh:
        fh.write("orf_id,ucsc,orf_set,is_smorf,blastoff_class,aa_seq,smorf_type,host_gene_name\n")
        for o, st, sm, bc, aa, ty, h in ORFS:
            fh.write(f'{o},"=HYPERLINK(""x"",""y"")",{st},{sm},{bc},{aa},{ty},{h}\n')
    r = subprocess.run([sys.executable, SCRIPT, "--orfs", j("orfs.csv"), "--bed12", j("o.bed12"),
                        "--genome", j("g.fa"), "--outdir", d], capture_output=True, text=True)
    check("exit", r.returncode, 0)
    if r.returncode:
        print(r.stderr)

    def fasta(p):
        out, k = {}, None
        for l in open(p):
            if l.startswith(">"):
                k = l[1:].split()[0]
                out[k] = ""
            else:
                out[k] += l.strip()
        return out

    un_nuc, un_pep = fasta(j("unreviewed_smorfs.nuc")), fasta(j("unreviewed_smorfs.pep"))
    rv_nuc = fasta(j("reviewed_smorfs.nuc"))
    check("unreviewed = all but sp_identical, Comprehensive smORFs only", list(un_nuc),
          ["orf_plus", "orf_minus", "orf_sp"])
    check("reviewed = sp_identical only", list(rv_nuc), ["orf_spid"])
    check("plus nuc incl. stop", un_nuc["orf_plus"], "ATGAAATAG")
    check("minus nuc reverse-complemented", un_nuc["orf_minus"], MINUS)
    check("pep", (un_pep["orf_plus"], un_pep["orf_minus"]), ("MK", "MA"))
    hdr = open(j("unreviewed_smorfs.pep")).readline()
    check("header", hdr.strip(), ">orf_plus chr1:4-12(+) smorf_type=uORF host_gene=Gene1 "
                                 "blastoff_class=novel orf_set=primary aa=2")

    rows = [l.rstrip("\n").split("\t") for l in open(j("unreviewed_smorfs.gtf"))]
    got = {(r[8].split('"')[1], r[2]): [] for r in rows}
    for r in rows:
        got[(r[8].split('"')[1], r[2])].append((int(r[3]), int(r[4]), r[7]))
    check("plus transcript", got[("orf_plus", "transcript")], [(4, 12, ".")])
    check("plus CDS excludes stop", got[("orf_plus", "CDS")], [(7 - 3, 9, "0")])
    check("plus start/stop", (got[("orf_plus", "start_codon")], got[("orf_plus", "stop_codon")]),
          ([(4, 6, "0")], [(10, 12, "0")]))
    # minus: 5'->3' positions are 36..30 then 21,20 (0-based); CDS = first 6 = 36..31;
    # stop = 30, 21, 20 -> split across both exons.
    check("minus exons", got[("orf_minus", "exon")], [(21, 22, "."), (31, 37, ".")])
    check("minus CDS", got[("orf_minus", "CDS")], [(32, 37, "0")])
    # transcription order: 31 is the stop's 1st base (frame 0), 21-22 its last two (frame 2)
    check("minus stop split", got[("orf_minus", "stop_codon")], [(31, 31, "0"), (21, 22, "2")])
    check("minus start", got[("orf_minus", "start_codon")], [(35, 37, "0")])
    check("GTF grouped per ORF, sorted", [r[8].split('"')[1] for r in rows][0], "orf_plus")
    check("attributes", 'smorf_type "uORF"; blastoff_class "novel"' in rows[0][8], True)

    # a disagreeing aa_seq must fail the run
    bad = open(j("orfs.csv")).read().replace(",MK,uORF,Gene1", ",MQ,uORF,Gene1")
    open(j("bad.csv"), "w").write(bad)
    r = subprocess.run([sys.executable, SCRIPT, "--orfs", j("bad.csv"), "--bed12", j("o.bed12"),
                        "--genome", j("g.fa"), "--outdir", d], capture_output=True, text=True)
    check("aa_seq mismatch fails", r.returncode != 0 and "MISMATCH orf_plus" in r.stderr, True)

print(f"\n{'ALL PASS' if not failed else f'{failed} FAILED'}")
sys.exit(1 if failed else 0)
