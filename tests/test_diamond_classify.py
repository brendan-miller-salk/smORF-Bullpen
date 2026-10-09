"""DIAMOND conservation classifier checks. Prints PASS/FAIL, exits 1 on failure.

Synthetic blastx rows (transcripts vs ORF proteins). Human has a 3 Mb filler gene,
so the whole transcriptome is ~1e6 aa, and an ortholog table mapping ENSMUSG1 to
ORTHO (300 nt = 100 aa). Rat has the same transcripts but no ortholog table.
Poly-A ORFs: a perfect 15-aa match is ~29.2 bits, E ~0.03 against the transcriptome
(untestable) but ~3e-6 against ORTHO (testable).
"""
import gzip
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "workflow", "scripts", "diamond_classify.py")
sys.path.insert(0, os.path.dirname(SCRIPT))
from diamond_classify import perfect_bits  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


def title(tid, gid, sym, biotype="protein_coding"):
    return (f">{tid} cdna chromosome:GRCh38:1:100:200:1 gene:{gid}.3 gene_biotype:{biotype} "
            f"transcript_biotype:{biotype}" + (f" gene_symbol:{sym}" if sym else "")
            + " description:x [Source:HGNC Symbol;Acc:HGNC:1]")


TRANSCRIPTS = [(title("ENST_O", "ENSG_O", "ORTHO"), 300),
               (title("ENST_BIG", "ENSG_BIG", "BIG"), 3_000_000),
               (title("ENST_NC", "ENSG_NC", "", "lncRNA"), 300)]
LETTERS = sum(n for _, n in TRANSCRIPTS) / 3


def hit(tid, orf, ss, se, sl, bit, frame=1, pident="90.0"):
    return "\t".join(map(str, [tid, orf, pident, se - ss + 1, 1, 90, 300, frame,
                                ss, se, sl, "1e-30", bit]))


PROTEINS = {"orf_long": "M" + "A" * 59, "orf_part": "M" + "A" * 59,
            "orf_short_orth": "A" * 15, "orf_short_noorth": "A" * 15,
            "orf_short_orth_fail": "A" * 15}
ORFS = [*PROTEINS, "orf_nostop"]
HUMAN = [
    hit("ENST_BIG", "orf_long", 1, 60, 60, 80),               # outside orthologs, passes
    hit("ENST_NC", "orf_long", 1, 60, 60, 90, frame=-2),      # reverse frame: ignored
    hit("ENST_BIG", "orf_part", 1, 10, 60, 70),               # cov 0.17: fails
    hit("ENST_NC", "orf_part", 1, 40, 60, 60),                # cov 0.67: passes, lower bit
    hit("ENST_O", "orf_short_orth", 1, 15, 15, 27),           # passes only as ortholog
    hit("ENST_BIG", "orf_short_noorth", 1, 15, 15, 27),       # fails, cannot pass: untestable
]
CATALOGUE = [("orf_long", ""), ("orf_part", ""), ("orf_short_orth", "MSTRG.1"),
             ("orf_short_noorth", "MSTRG.2"), ("orf_short_orth_fail", "ENSMUSG1.5"),
             ("orf_nostop", "ENSMUSG1")]
ANNOTATOR = [("orf_short_orth", "ENSMUSG1.4"), ("orf_short_orth_fail", "")]

with tempfile.TemporaryDirectory() as d:
    j = lambda x: os.path.join(d, x)  # noqa: E731
    with gzip.open(j("human.tsv.gz"), "wt") as fh:   # DIAMOND --compress 1
        fh.write("\n".join(HUMAN) + "\n")
    open(j("rat.tsv"), "w").close()
    with gzip.open(j("t.fa.gz"), "wt") as fh:
        for t, n in TRANSCRIPTS:
            fh.write(t + "\n" + ("ACGT" * (n // 4)) + "\n")
    open(j("orth.tsv"), "w").write("mouse_gene_id\thuman_gene_id\torthology_type\n"
                                   "ENSMUSG1\tENSG_O\tortholog_one2one\nENSMUSG9\tENSG_BIG\tx\n")
    open(j("orfs.faa"), "w").write("".join(f">{k}\n{v}\n" for k, v in PROTEINS.items()))
    open(j("cat.tsv"), "w").write("orf_id\tchrom\tgene_id\n"
                                  + "".join(f"{o}\tchr1\t{g}\n" for o, g in CATALOGUE))
    open(j("ann.tsv"), "w").write("orf_id\tsmorf_type\thost_gene_id\n"
                                  + "".join(f"{o}\tx\t{g}\n" for o, g in ANNOTATOR))
    with open(j("o.bed12"), "w") as fh:
        for oid in ORFS:
            fh.write(f"chr1\t0\t30\t{oid}\t0\t+\t0\t30\t0\t1\t30,\t0,\n")
    r = subprocess.run([sys.executable, SCRIPT, "--hits", f"human={j('human.tsv.gz')}",
                        f"rat={j('rat.tsv')}", "--transcripts", f"human={j('t.fa.gz')}",
                        f"rat={j('t.fa.gz')}", "--orthologs", f"human={j('orth.tsv')}",
                        "--proteins", j("orfs.faa"), "--catalogue", j("cat.tsv"),
                        "--annotator", j("ann.tsv"), "--bed12", j("o.bed12"),
                        "--evalue", "1e-3", "--min-cov", "0.5", "--out", j("o.tsv")],
                       capture_output=True, text=True)
    check("classify exit", r.returncode, 0)
    if r.returncode:
        print(r.stderr)
    rows = [l.rstrip("\n").split("\t") for l in open(j("o.tsv"))]
    hdr, body = rows[0], {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}

check("orf_id first", hdr[0], "orf_id")
check("one row per ORF, bed order", [r[0] for r in rows[1:]], ORFS)
check("per-species columns", hdr[4:12], [f"diamond_human_{c}" for c in
      ("call", "ortholog", "gene", "biotype", "pident", "cov", "evalue", "ortholog_evalue")])
check("perfect 15 A bits (+1.5 slack)", round(perfect_bits("A" * 15), 1), 29.2)

c = body["orf_long"]
check("long conserved outside orthologs", (c["diamond_human_call"], c["diamond_species"]),
      ("conserved", "human"))
check("long reverse frame ignored", c["diamond_human_gene"], "BIG")
check("long tblastn E", c["diamond_human_evalue"], f"{60 * LETTERS * 2.0 ** -80:.2e}")
check("long no ortholog E", c["diamond_human_ortholog_evalue"], "")
check("long rat: no hit but testable", c["diamond_rat_call"], "not_conserved")

p = body["orf_part"]
check("passing hit beats higher bitscore", (p["diamond_human_gene"], p["diamond_human_cov"]),
      ("ENSG_NC", "0.667"))
check("no symbol -> gene id; biotype", (p["diamond_human_gene"], p["diamond_human_biotype"]),
      ("ENSG_NC", "lncRNA"))

s = body["orf_short_orth"]
check("short conserved via ortholog", s["diamond_human_call"], "conserved")
check("annotator host gene -> ortholog", s["diamond_human_ortholog"], "ORTHO")
check("short transcriptome E fails", float(s["diamond_human_evalue"]) > 1e-3, True)
check("short ortholog E", s["diamond_human_ortholog_evalue"], f"{15 * 100 * 2.0 ** -27:.2e}")
check("short rat: no orthologs -> untestable", s["diamond_rat_call"], "untestable")

u = body["orf_short_noorth"]
check("short without orthologs untestable", u["diamond_human_call"], "untestable")
check("untestable still reports its hit", u["diamond_human_gene"], "BIG")

f = body["orf_short_orth_fail"]
check("catalogue gene_id fallback -> ortholog", f["diamond_human_ortholog"], "ORTHO")
check("short testable in ortholog, no hit", f["diamond_human_call"], "not_conserved")

n = body["orf_nostop"]
check("no protein untestable", (n["diamond_n_species"], n["diamond_n_untestable"],
                                n["diamond_human_call"]), ("0", "2", "untestable"))

print(f"\n{'ALL PASS' if not failed else f'{failed} FAILED'}")
sys.exit(1 if failed else 0)
