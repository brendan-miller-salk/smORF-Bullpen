"""Synthetic checks for build_catalogue.py coordinate handling. Prints PASS/FAIL, exits 1 on failure."""
import collections
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
from build_catalogue import (SeqnameMap, bh_p_cutoff, cluster_by_overlap, duplicate_inputs, merge,  # noqa: E402
                             overlap, price_gene_from_transcript, project, read_ribotish, same_junctions)

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


# Exons 100-110 and 200-220 (0-based half-open); transcript length 30.
ex = [(100, 110), (200, 220)]
check("plus, within exon 1", project(ex, "+", 2, 8), ((102, 108),))
check("plus, spans junction", project(ex, "+", 5, 15), ((105, 110), (200, 205)))
# Minus strand: transcript 5' end is genomic 220.
check("minus, within last exon", project(ex, "-", 0, 6), ((214, 220),))
check("minus, spans junction", project(ex, "-", 15, 25), ((105, 110), (200, 205)))
check("past transcript end -> None", project(ex, "+", 25, 35), None)

check("overlap across blocks", overlap(((0, 10), (20, 30)), ((5, 25),)), 10)

m = SeqnameMap(["chr1", "chrM", "GL456233.2"])
check("bare -> chr", m("1"), "chr1")
check("MT -> chrM", m("MT"), "chrM")
check("M -> chrM", m("M"), "chrM")
check("scaffold kept", m("GL456233.2"), "GL456233.2")
check("unknown -> None", m("99"), None)

P = lambda k, b, cls="novel_u", aa=10: dict(key=k, blocks=b, orf_class=cls, aa_length=aa)
cl = cluster_by_overlap([P("a", ((0, 100),)), P("b", ((5, 100),)), P("c", ((50, 150),))], 0.8)
check("0.8 reciprocal overlap clusters", sorted(len(c) for c in cl), [1, 2])

# Benjamini-Hochberg: p_(i) <= i/m * fdr, keep up to the largest such i.
check("bh keeps step-up", bh_p_cutoff([0.001, 0.008, 0.039, 0.041, 0.9], 0.1), 0.041)
check("bh none pass", bh_p_cutoff([0.5, 0.6, 0.9], 0.1), -1.0)
check("bh all pass", bh_p_cutoff([0.0, 0.01, 0.02], 0.1), 0.02)

with tempfile.TemporaryDirectory() as d:
    def w(name, text):
        p = os.path.join(d, name)
        open(p, "w").write(text)
        return p
    a1, a2, b, e1, e2 = (w("a1.txt", "x\n"), w("a2.txt", "x\n"), w("b.txt", "y\n"),
                         w("e1.txt", ""), w("e2.txt", ""))
    check("identical samples flagged", duplicate_inputs([a1, a2, b]), [[a1, a2]])
    check("distinct samples pass", duplicate_inputs([a1, b]), [])
    check("empty samples ignored", duplicate_inputs([e1, e2, b]), [])


# Splice isoforms stay apart; start-codon variants in the same exon structure merge.
unspliced = ((1860, 2211),)                       # Dnajc19-ps 116 aa
spliced = ((1860, 2156), (2208, 2245))            # MSTRG.55403.1 110 aa: 52-nt intron inside
later_start = ((1860, 2150),)                     # same exon, downstream start
check("junction inside the shared span differs", same_junctions(unspliced, spliced), False)
check("start variant in the same exon", same_junctions(unspliced, later_start), True)
check("ends inside the other's intron", same_junctions(((100, 190),), ((100, 150), (200, 290))), False)
check("start in an upstream exon vs downstream start", same_junctions(((0, 50), (100, 200)), ((120, 200),)), True)
cl2 = cluster_by_overlap([P("u", unspliced), P("s", spliced), P("l", later_start)], 0.8)
check("spliced and unspliced never share a cluster",
      any({"u", "s"} <= {p["key"] for p in c} for c in cl2), False)
ins_a = ((0, 300), (400, 700))
ins_b = ((0, 300), (350, 392), (400, 700))        # extra 42-nt exon: 93% overlap
check("extra internal exon not merged", len(cluster_by_overlap([P("a", ins_a), P("b", ins_b)], 0.8)), 2)
check("non-overlapping ORFs are not the same", same_junctions(((0, 50),), ((60, 120),)), False)
check("abutting ORFs are not the same", same_junctions(((0, 60),), ((60, 120),)), False)


# --merge exact: one row per distinct structure; both callers credited on an exact match.
def call(caller, blocks, tx="T1", cls="uORF", nid="x"):
    return dict(chrom="chr1", strand="+", blocks=blocks, orf_class=cls, caller=caller,
                native_id=nid, transcript_id=tx, gene_id="G1", bed_score=0)


calls = [call("price", ((1000, 1063),)), call("ribocode", ((1000, 1063),), tx="T2"),
         call("ribocode", ((1012, 1063),)),                 # downstream start, same stop
         call("price", ((2000, 2600),), cls="canonical_cds"),
         call("ribocode", ((2003, 2600),), cls="canonical_cds")]   # 1 codon shorter
ex = merge(calls, 0.8, "exact")
check("exact: one cluster per structure", len(ex), 4)
same = [c for c in ex if c[0]["blocks"] == ((1000, 1063),)][0]
check("exact: identical calls from both callers share a row",
      sorted(m["caller"] for m in same[0]["members"]), ["price", "ribocode"])
check("overlap: the canonical start variants fold", len(merge(calls, 0.8, "overlap")), 3)

# PRICE's Gene column is replaced by the gene of the transcript PRICE names.
tx = {"T1": ("chr7", "-", [(0, 10)], "G_pseudo"), "T2": ("chr7", "-", [(0, 10)], "G2")}
rows = [dict(transcript_id="T1", native_id="T1_ncRNA_0", gene_id="G_olfr_G_hbb"),
        dict(transcript_id="T2", native_id="T2_CDS_0", gene_id="G2")]
n = price_gene_from_transcript(rows, tx)
check("price gene from transcript", [r["gene_id"] for r in rows], ["G_pseudo", "G2"])
check("price genes changed", n, 1)

print("\nALL TESTS PASSED" if not failed else f"\n{failed} FAILED")


# Ribo-TISH predict --blocks: GenomePos chr:start-end:strand, Blocks 0-based half-open,
# stop codon included; ':CDSFrameOverlap' / ':Known' suffixes do not change the class.
RT_HEADER = ("Gid\tTid\tSymbol\tGeneType\tGenomePos\tStartCodon\tStart\tStop\tTisType\tTISGroup\t"
             "TISCounts\tTISPvalue\tRiboPvalue\tRiboPStatus\tFisherPvalue\tTISQvalue\tFrameQvalue\t"
             "FisherQvalue\tAALen\tAASeq\tBlocks\n")


def rt_row(gpos, ttype, p, blocks, codon="ATG", tis_q="0.01", aalen=None):
    n = sum(int(b.split("-")[1]) - int(b.split("-")[0]) for b in blocks.split(","))
    aalen = n // 3 - 1 if aalen is None else aalen
    return "\t".join(["G1", "T1", "SYM", "protein_coding", gpos, codon, "10", "40", ttype, "1", "5",
                      "1e-4", "1e-3", "T", str(p), tis_q, "0.01", "0.01", str(aalen), "MAAAAAAAA", blocks]) + "\n"


with tempfile.TemporaryDirectory() as d:
    path = os.path.join(d, "HMC3_LTM_rep1.txt")
    open(path, "w").write(RT_HEADER + rt_row("1:101-160:+", "5'UTR:CDSFrameOverlap", 1e-6, "100-110,200-220")
                          + rt_row("chr1:301-330:-", "Novel", 0.02, "300-330", "CTG"))
    rows = read_ribotish(path, "HMC3_LTM_rep1", SeqnameMap(["chr1"]), 0.05, collections.Counter())
    check("ribotish: rows", len(rows), 2)
    check("ribotish: chrom harmonised, strand", (rows[0]["chrom"], rows[0]["strand"]), ("chr1", "+"))
    check("ribotish: blocks as given", rows[0]["blocks"], ((100, 110), (200, 220)))
    check("ribotish: class without suffix", (rows[0]["orf_class"], rows[1]["orf_class"]), ("uORF", "novel_u"))
    check("ribotish: p and sample", (rows[0]["pval"], rows[0]["sample"]), (1e-6, "HMC3_LTM_rep1"))
    check("ribotish: native id keeps start codon", rows[1]["native_id"], "T1:10-40:CTG")
    open(path, "w").write(RT_HEADER + rt_row("chr1:1-30:+", "Annotated", 1e-9, "0-30", tis_q="0.2")
                          + rt_row("chr1:1-30:+", "Annotated", 1e-9, "0-30", aalen=10)
                          + rt_row("chr1:1-30:+", "Annotated", 1e-9, "0-30", tis_q="None"))
    st = collections.Counter()
    check("ribotish: TIS q above cut-off, no stop, no TIS q -> dropped",
          len(read_ribotish(path, "S", SeqnameMap(["chr1"]), 0.05, st)), 0)
    check("ribotish: drops counted", (st["ribotish_dropped_tis_q:S"], st["ribotish_dropped_no_stop:S"]), (2, 1))
    check("ribotish: no TIS filter keeps the stopped calls",
          len(read_ribotish(path, "S", SeqnameMap(["chr1"]), None, collections.Counter())), 2)

sys.exit(1 if failed else 0)
