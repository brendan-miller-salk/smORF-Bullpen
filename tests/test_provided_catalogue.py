"""Synthetic checks for provided_catalogue.py (provided mode). Prints PASS/FAIL, exits 1 on failure."""
import contextlib
import io
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
from provided_catalogue import main  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


# chr1, filler C. A (+): ATGAAAA | intron | AA + stop_codon row TAA, 0-based [10,17) and [40,45).
# B (-): sense ATGCCC then TAG, no stop_codon row -> extended; genome [100,109) = CTA GGGCAT,
# its CDS 1-based 104-109.
g = ["C"] * 200
g[10:17] = "ATGAAAA"
g[40:45] = "AATAA"
g[100:109] = "CTAGGGCAT"
GENOME = ">chr1\n" + "".join(g) + "\n"


def gtf_row(chrom, feat, s1, e1, strand, attrs):
    return "\t".join([chrom, "t", feat, str(s1), str(e1), ".", strand, ".", attrs]) + "\n"


BASE = [
    gtf_row("1", "transcript", 11, 45, "+", 'gene_id "gA"; transcript_id "tA";'),
    gtf_row("1", "CDS", 11, 17, "+", 'gene_id "gA"; transcript_id "tA";'),
    gtf_row("1", "CDS", 41, 42, "+", 'gene_id "gA"; transcript_id "tA";'),
    gtf_row("1", "stop_codon", 43, 45, "+", 'gene_id "gA"; transcript_id "tA";'),
    # same structure again, own orf_id attribute: folded into A's row
    gtf_row("1", "transcript", 11, 45, "+", 'gene_id "gA"; transcript_id "tA2"; orf_id "myORF_7";'),
    gtf_row("1", "CDS", 11, 17, "+", 'gene_id "gA"; transcript_id "tA2";'),
    gtf_row("1", "CDS", 41, 45, "+", 'gene_id "gA"; transcript_id "tA2";'),
    gtf_row("1", "transcript", 104, 109, "-", 'gene_id "gB"; transcript_id "tB";'),
    gtf_row("1", "CDS", 104, 109, "-", 'gene_id "gB"; transcript_id "tB";'),
    gtf_row("1", "transcript", 150, 190, "+", 'gene_id "gN"; transcript_id "tNC";'),
    gtf_row("1", "exon", 150, 190, "+", 'gene_id "gN"; transcript_id "tNC";'),
]


def run(d, rows):
    open(os.path.join(d, "g.fa"), "w").write(GENOME)
    open(os.path.join(d, "o.gtf"), "w").writelines(rows)
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err):
            code = main(["--gtf", os.path.join(d, "o.gtf"), "--genome", os.path.join(d, "g.fa"),
                         "--out-prefix", os.path.join(d, "cat"), "--report", os.path.join(d, "rep.tsv")])
    except SystemExit as e:
        return str(e.code)
    return code


with tempfile.TemporaryDirectory() as d:
    check("runs", run(d, BASE), 0)
    bed = [l.split("\t") for l in open(os.path.join(d, "cat.bed12")).read().splitlines()]
    tsv = [l.split("\t") for l in open(os.path.join(d, "cat.tsv")).read().splitlines()]
    rep = dict(l.split("\t") for l in open(os.path.join(d, "rep.tsv")).read().splitlines()[1:])
    check("duplicate folded, non-coding ignored: 2 ORFs", len(bed), 2)
    check("bare seqname -> genome's", bed[0][0], "chr1")
    check("stop_codon row merged into last block", bed[0][1:3] + bed[0][9:12],
          ["10", "45", "2", "7,5,", "0,30,"])
    check("minus strand extended 3 nt to its stop", bed[1][1:3], ["100", "109"])
    check("consensus == catalogue", open(os.path.join(d, "cat.consensus.bed12")).read(),
          open(os.path.join(d, "cat.bed12")).read())
    h = tsv[0]
    row = [dict(zip(h, r)) for r in tsv[1:]]
    check("orf ids", [r["orf_id"] for r in row], ["orf_00000001", "orf_00000002"])
    check("source_id: orf_id attribute, else transcript_id, joined", row[0]["source_id"], "myORF_7,tA")
    check("transcript_ids follow source order", row[0]["transcript_id"], "tA2,tA")
    check("aa_length excludes the stop", [r["aa_length"] for r in row], ["3", "2"])
    check("orf_class", {r["orf_class"] for r in row}, {"provided"})
    check("report", (rep["orfs"], rep["duplicate_structures_folded"], rep["extended_3nt_to_stop"]),
          ("2", "1", "1"))

with tempfile.TemporaryDirectory() as d:
    rows = BASE[:4] + [gtf_row("1", "CDS", 60, 61, "+", 'gene_id "gX"; transcript_id "tX";')]
    check("CDS without transcript row aborts", "no transcript row" in str(run(d, rows)), True)
with tempfile.TemporaryDirectory() as d:
    rows = BASE[:3] + [gtf_row("1", "stop_codon", 43, 44, "+", 'gene_id "gA"; transcript_id "tA";')]
    check("not a multiple of 3 aborts", "multiple of 3" in str(run(d, rows)), True)
with tempfile.TemporaryDirectory() as d:
    rows = [gtf_row("1", "transcript", 151, 159, "+", 'gene_id "gC"; transcript_id "tC";'),
            gtf_row("1", "CDS", 151, 159, "+", 'gene_id "gC"; transcript_id "tC";')]
    check("no stop aborts", "stop codon" in str(run(d, rows)), True)
with tempfile.TemporaryDirectory() as d:
    rows = [r.replace("1\tt", "9\tt", 1) for r in BASE[:4]]
    check("unknown seqname aborts", "not in the genome" in str(run(d, rows)), True)

print(f"\n{'ALL PASS' if not failed else f'{failed} FAILED'}")
sys.exit(1 if failed else 0)
