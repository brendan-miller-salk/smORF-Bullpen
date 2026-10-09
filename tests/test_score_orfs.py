"""Synthetic tests for score_orfs.py: known inputs, hand-computed expectations."""
import sys, os, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
from score_orfs import PsiteIndex, Orf, score_orf, parse_bed12, cds_from_gtf, calibrate, add_ppm, score_many, CdsMask
import pandas as pd

tmp = tempfile.mkdtemp()
fails = []
def check(name, got, want, tol=1e-9):
    ok = (want is None and got is None) or (isinstance(want, float) and np.isnan(want) and np.isnan(got)) or \
         (isinstance(got, (int, float, np.floating)) and abs(got - want) < tol) or (got == want)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")
    if not ok: fails.append(name)

def make_bg(path, rows):
    with open(path, "w") as fh:
        for chrom, start, end, val in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{val}\n")

# ---------------------------------------------------------------- test 1
# Perfect + strand ORF, 30 nt (10 codons incl stop). All P-sites in frame 0.
fwd = os.path.join(tmp, "f.bedgraph"); rev = os.path.join(tmp, "r.bedgraph")
rows = [("chr1", 100 + 3*i, 101 + 3*i, 10) for i in range(10)]   # frame 0 only
make_bg(fwd, rows); make_bg(rev, [])
idx = PsiteIndex.from_bedgraphs([fwd], [rev])
orf = Orf("perfect", "chr1", "+", np.arange(100, 130, dtype=np.int64), {})
r = score_orf(orf, idx)
check("t1 PIF (all in frame)", r["pif"], 1.0)
check("t1 Uniformity (all codons covered in frame)", r["uniformity"], 1.0)
check("t1 psites_total", r["psites_total"], 100.0)
check("t1 n_codons", r["n_codons"], 10)

# ---------------------------------------------------------------- test 2
# Pure noise: equal counts at every position -> PIF = 1/3, no codon exceeds 33.3%
rows = [("chr1", 100 + i, 101 + i, 9) for i in range(30)]
make_bg(fwd, rows)
idx = PsiteIndex.from_bedgraphs([fwd], [rev])
r = score_orf(orf, idx)
check("t2 PIF (uniform noise)", r["pif"], 1/3)
check("t2 Uniformity (no codon >33.33%)", r["uniformity"], 0.0)

# ---------------------------------------------------------------- test 3
# Half the codons in frame, half zero-coverage -> zero-codon policy matters.
rows = [("chr1", 100 + 3*i, 101 + 3*i, 10) for i in range(5)]   # first 5 codons only
make_bg(fwd, rows)
idx = PsiteIndex.from_bedgraphs([fwd], [rev])
check("t3 Uniformity zero=fail", score_orf(orf, idx, "fail")["uniformity"], 0.5)
check("t3 Uniformity zero=ignore", score_orf(orf, idx, "ignore")["uniformity"], 1.0)

# ---------------------------------------------------------------- test 4
# Minus strand mirror of test 1. Transcript order descends in genomic coords,
# so frame-0 positions are 129, 126, ... and plastid writes negative values.
rows = [("chr1", 129 - 3*i, 130 - 3*i, -10) for i in range(10)]
make_bg(fwd, []); make_bg(rev, rows)
idx = PsiteIndex.from_bedgraphs([fwd], [rev])
orf_m = Orf("minus", "chr1", "-", np.arange(100, 130, dtype=np.int64)[::-1], {})
r = score_orf(orf_m, idx)
check("t4 minus-strand PIF", r["pif"], 1.0)
check("t4 minus-strand Uniformity", r["uniformity"], 1.0)
check("t4 minus-strand abs() applied", r["psites_total"], 100.0)

# ---------------------------------------------------------------- test 5
# Drop-off. Stop codon = last 3 nt (127,128,129). Before-window = 15 nt at
# 112..126, after-window = 130..144. Put 10 in-frame before, 0 after -> 1.0
before = [("chr1", 112 + 3*i, 113 + 3*i, 10) for i in range(5)]  # 112,115,118,121,124 = frame 0
make_bg(fwd, before); make_bg(rev, [])
idx = PsiteIndex.from_bedgraphs([fwd], [rev])
check("t5 dropoff all before", score_orf(orf, idx)["dropoff"], 1.0)

after = [("chr1", 130 + 3*i, 131 + 3*i, 10) for i in range(5)]   # frame 0 continues downstream
make_bg(fwd, before + after)
idx = PsiteIndex.from_bedgraphs([fwd], [rev])
check("t5 dropoff half after", score_orf(orf, idx)["dropoff"], 0.5)

make_bg(fwd, after)
idx = PsiteIndex.from_bedgraphs([fwd], [rev])
check("t5 dropoff none before", score_orf(orf, idx)["dropoff"], 0.0)

# ---------------------------------------------------------------- test 6
# Empty ORF -> NaNs, no crash.
make_bg(fwd, []); make_bg(rev, [])
idx = PsiteIndex.from_bedgraphs([fwd], [rev])
r = score_orf(orf, idx)
check("t6 no coverage -> NaN pif", np.isnan(r["pif"]), True)

# ---------------------------------------------------------------- test 7
# BED12 parsing, spliced 2-block ORF on minus strand.
bed = os.path.join(tmp, "orfs.bed12")
with open(bed, "w") as fh:
    # chromStart=100, blocks of 6 and 9 at offsets 0 and 20 -> 100-105, 120-128
    fh.write("chr1\t100\t129\tsplicedORF\t0\t-\t100\t129\t0\t2\t6,9,\t0,20,\n")
orfs = list(parse_bed12(bed))
check("t7 one ORF parsed", len(orfs), 1)
check("t7 spliced length", orfs[0].positions.size, 15)
check("t7 minus-strand first pos is highest", int(orfs[0].positions[0]), 128)
check("t7 block gap excluded", int(orfs[0].positions[-1]), 100)
check("t7 106-119 not present", bool(np.any((orfs[0].positions > 105) & (orfs[0].positions < 120))), False)

# ---------------------------------------------------------------- test 8
# GTF CDS extraction, one-per-gene picks the longest; the stop codon is appended, and
# incomplete CDS (NF tags, no stop, non-zero first phase, length not % 3) are skipped.
gtf = os.path.join(tmp, "t.gtf")
with open(gtf, "w") as fh:
    fh.write('chr1\ts\tCDS\t101\t130\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n')
    fh.write('chr1\ts\tstop_codon\t131\t133\t.\t+\t0\tgene_id "G1"; transcript_id "T1";\n')
    fh.write('chr1\ts\tCDS\t101\t160\t.\t+\t0\tgene_id "G1"; transcript_id "T2";\n')
    fh.write('chr1\ts\tstop_codon\t161\t163\t.\t+\t0\tgene_id "G1"; transcript_id "T2";\n')
    fh.write('chr1\ts\tCDS\t101\t190\t.\t+\t0\tgene_id "G1"; transcript_id "T2nf"; tag "cds_end_NF";\n')
    fh.write('chr2\ts\tCDS\t204\t230\t.\t-\t0\tgene_id "G2"; transcript_id "T3";\n')
    fh.write('chr2\ts\tstop_codon\t201\t203\t.\t-\t0\tgene_id "G2"; transcript_id "T3";\n')
    fh.write('chr3\ts\tCDS\t301\t330\t.\t+\t0\tgene_id "G3"; transcript_id "T4";\n')   # no stop
    fh.write('chr3\ts\tCDS\t401\t430\t.\t+\t1\tgene_id "G4"; transcript_id "T5";\n')   # phase 1
    fh.write('chr3\ts\tstop_codon\t431\t433\t.\t+\t0\tgene_id "G4"; transcript_id "T5";\n')
    fh.write('chr3\ts\tCDS\t501\t531\t.\t+\t0\tgene_id "G5"; transcript_id "T6";\n')   # 31 nt
    fh.write('chr3\ts\tstop_codon\t532\t534\t.\t+\t0\tgene_id "G5"; transcript_id "T6";\n')
c = sorted(cds_from_gtf(gtf), key=lambda o: o.name)
check("t8 complete CDS only, one per gene", [o.name for o in c], ["T2", "T3"])
check("t8 GTF 1-based -> 0-based", int(c[0].positions[0]), 100)
check("t8 length incl. stop", c[0].positions.size, 63)
check("t8 minus: starts at the CDS start", int(c[1].positions[0]), 229)
check("t8 minus: ends with the stop", int(c[1].positions[-1]), 200)

# ---------------------------------------------------------------- test 9
# Threshold calibration direction.
rng = np.random.default_rng(0)
cds_df = pd.DataFrame({
    "pif": np.clip(rng.normal(0.85, 0.05, 2000), 0, 1),
    "uniformity": np.clip(rng.normal(0.80, 0.06, 2000), 0, 1),
    "dropoff": np.clip(rng.normal(0.90, 0.04, 2000), 0, 1),
    "psites_total": np.full(2000, 100.0),
})
lo = calibrate(cds_df, tail="lower")
hi = calibrate(cds_df, tail="upper")
check("t9 lower tail below mean", lo["pif"] < lo["pif_mu"], True)
check("t9 upper tail above mean", hi["pif"] > hi["pif_mu"], True)
check("t9 lower ~ mu-1.645sigma", abs(lo["pif"] - (lo["pif_mu"] - 1.6449*lo["pif_sigma"])) < 0.002, True)
check("t9 dropoff uses mean", abs(lo["dropoff"] - lo["dropoff_mu"]) < 1e-9, True)
frac_cds_passing = float((cds_df["pif"] >= lo["pif"]).mean())
check("t9 lower tail passes ~95% of CDS", 0.92 < frac_cds_passing < 0.98, True)

# ---------------------------------------------------------------- test 10
# PPM sums to 1e6 across scored features.
df = pd.DataFrame({"psites_total": [100.0, 200.0, 0.0], "length_nt": [30, 60, 30]})
df = add_ppm(df)
check("t10 PPM sums to 1e6", abs(df["ppm"].sum() - 1e6) < 1e-6, True)

print()

# ---------------------------------------------------------------- unshared evidence
# Annotated CDS chr2:1-12 (+) and its stop 13-15 (GTF). ORF A = [0, 21): an isoform whose
# last 6 nt lie outside the CDS. ORF B = the CDS exactly. ORF C on the - strand is not
# covered by the + strand CDS.
gtf2 = os.path.join(tmp, "mask.gtf")
with open(gtf2, "w") as fh:
    fh.write('chr2\tx\tCDS\t1\t12\t.\t+\t0\ttranscript_id "T"; gene_id "G";\n')
    fh.write('chr2\tx\tstop_codon\t13\t15\t.\t+\t0\ttranscript_id "T"; gene_id "G";\n')
mask = CdsMask(gtf2)
f2 = os.path.join(tmp, "f2.bedgraph"); r2 = os.path.join(tmp, "r2.bedgraph")
make_bg(f2, [("chr2", 0, 1, 50), ("chr2", 15, 16, 4), ("chr2", 16, 17, 7), ("chr2", 18, 19, 3)])
make_bg(r2, [])
idx2 = PsiteIndex.from_bedgraphs([f2], [r2])
ra = score_orf(Orf("A", "chr2", "+", np.arange(0, 21, dtype=np.int64), {}), idx2, cds_mask=mask)
check("unshared nt: outside CDS+stop", ra["unshared_nt"], 6)
check("unshared psites", ra["unshared_psites"], 14.0)
check("unshared in-frame psites (positions 15, 18)", ra["unshared_inframe_psites"], 7.0)
rb = score_orf(Orf("B", "chr2", "+", np.arange(0, 15, dtype=np.int64), {}), idx2, cds_mask=mask)
check("annotated CDS itself has no unshared nt", rb["unshared_nt"], 0)
rc = score_orf(Orf("C", "chr2", "-", np.arange(0, 9, dtype=np.int64)[::-1], {}), idx2, cds_mask=mask)
check("other strand is unshared", rc["unshared_nt"], 9)
# at least N in-frame P-sites in the unshared sequence
from score_orfs import apply_unshared
tab = pd.DataFrame(dict(unshared_nt=[0, 6, 6, 6], unshared_inframe_psites=[0.0, 7.0, 4.0, 0.0],
                        pass_all=[True, True, True, False]))
check("min 5: ORFs failing", apply_unshared(tab, "drop", 5), 2)
check("min 5: pass_unshared", tab["pass_unshared"].tolist(), [True, True, False, False])
check("min 5: pass_all", tab["pass_all"].tolist(), [True, True, False, False])
check("min 5: reason", tab.at[2, "exclude_reason"], "fewer than 5 in-frame P-sites outside annotated CDS")
t2 = tab.assign(pass_all=True)
apply_unshared(t2, "report", 5)
check("report mode leaves pass_all", t2["pass_all"].tolist(), [True, True, True, True])
# unique-read track: only its P-sites count toward *_unique
fu = os.path.join(tmp, "fu.bedgraph")
make_bg(fu, [("chr2", 15, 16, 4)])          # of ORF A's unshared P-sites, only pos 15 is unique
idxu = PsiteIndex.from_bedgraphs([fu], [r2])
rau = score_orf(Orf("A", "chr2", "+", np.arange(0, 21, dtype=np.int64), {}), idx2, cds_mask=mask, psites_unique=idxu)
check("unique P-sites in ORF", rau["psites_unique"], 4.0)
check("unique unshared in-frame", rau["unshared_inframe_psites_unique"], 4.0)
check("all-read columns unchanged", rau["unshared_inframe_psites"], 7.0)

# pooling: overlapping intervals from two samples add (sample 1 [0,10) = 2, sample 2 [5,6) = 3)
p1, p2, pr = (os.path.join(tmp, x) for x in ("p1.bg", "p2.bg", "pr.bg"))
make_bg(p1, [("chr3", 0, 10, 2)]); make_bg(p2, [("chr3", 5, 6, 3)]); make_bg(pr, [])
idxp = PsiteIndex.from_bedgraphs([p1, p2], [pr, pr])
rp = score_orf(Orf("P", "chr3", "+", np.arange(0, 12, dtype=np.int64), {}), idxp)
check("pooled overlapping intervals sum", rp["psites_total"], 23.0)

# frame-aware mask: sequence a CDS covers in ANOTHER frame is the ORF's own
blind = CdsMask(gtf2, frame_aware=False)
shift = Orf("D", "chr2", "+", np.arange(1, 16, dtype=np.int64), {})       # CDS frame + 1
check("frame-aware: other frame is unshared", score_orf(shift, idx2, cds_mask=mask)["unshared_nt"], 15)
check("frame-blind: other frame is shared", score_orf(shift, idx2, cds_mask=blind)["unshared_nt"], 1)
inner = Orf("E", "chr2", "+", np.arange(3, 21, dtype=np.int64), {})       # same frame, later start
check("frame-aware: same frame is shared", score_orf(inner, idx2, cds_mask=mask)["unshared_nt"], 6)
# minus strand CDS chr4:4-15 (-) + stop 1-3; spliced + strand CDS chr5:1-10 + 21-31 (phase 2) + stop 32-34
gtf3 = os.path.join(tmp, "mask3.gtf")
with open(gtf3, "w") as fh:
    fh.write('chr4\tx\tCDS\t4\t15\t.\t-\t0\ttranscript_id "M"; gene_id "GM";\n')
    fh.write('chr4\tx\tstop_codon\t1\t3\t.\t-\t0\ttranscript_id "M"; gene_id "GM";\n')
    fh.write('chr5\tx\tCDS\t1\t10\t.\t+\t0\ttranscript_id "S"; gene_id "GS";\n')
    fh.write('chr5\tx\tCDS\t21\t31\t.\t+\t2\ttranscript_id "S"; gene_id "GS";\n')
    fh.write('chr5\tx\tstop_codon\t32\t34\t.\t+\t0\ttranscript_id "S"; gene_id "GS";\n')
m3 = CdsMask(gtf3)
def un(chrom, strand, pos):
    return score_orf(Orf("x", chrom, strand, np.asarray(pos, dtype=np.int64), {}), idx2, cds_mask=m3)["unshared_nt"]
check("minus: the CDS itself", un("chr4", "-", range(14, -1, -1)), 0)
check("minus: same frame, later start", un("chr4", "-", range(11, -1, -1)), 0)
check("minus: other frame", un("chr4", "-", range(13, -1, -1)), 14)
check("spliced: same frame in exon 2 (phase 2)", un("chr5", "+", range(22, 34)), 0)
check("spliced: other frame in exon 2", un("chr5", "+", range(21, 33)), 12)
# body-codon density (above_cds_median): ORF A = 7 codons; body = codons 1-5 (start, stop out);
# codons 0-4 lie in the CDS (frame-blind), so only codon 5 (nt 15-17) is host-free.
check("body codons", rau["body_codons"], 5)
check("body in-frame unique (nt 3..15)", rau["body_inframe_psites_unique"], 4.0)
check("free body codons", rau["free_body_codons"], 1)
check("free body in-frame unique", rau["free_body_inframe_psites_unique"], 4.0)
check("in-frame unique", rau["psites_inframe_unique"], 4.0)
rsh = score_orf(shift, idx2, cds_mask=mask, psites_unique=idxu)          # CDS other frame: still not free
check("free = frame-blind, even with the frame-aware mask", rsh["free_body_codons"], 0)
from score_orfs import cds_median_density
cdst = pd.DataFrame(dict(psites_total=[200, 150, 50], psites_inframe_unique=[100.0, 60.0, 40.0],
                         body_codons=[10, 20, 5], body_inframe_psites_unique=[80.0, 40.0, 30.0]))
d = cds_median_density(cdst, 100)          # densities 8 and 2 per codon over 200 in-frame -> 4e4, 1e4 per million
check("CDS median density (calibration CDS only)", d["cds_median_density_unique"], 25000.0)
check("CDS density denominator: all CDS", d["inframe_cds_psites_unique"], 200.0)
check("CDS density n", d["n_cds_density"], 2)
check("no mask: no columns", "unshared_nt" in score_orf(Orf("A", "chr2", "+", np.arange(0, 21, dtype=np.int64), {}), idx2), False)

print(f"{'ALL TESTS PASSED' if not fails else 'FAILURES: ' + ', '.join(fails)}")
sys.exit(1 if fails else 0)
