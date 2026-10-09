"""Counting: count GTF (stop trim, flanks, subtraction), flank fit, counts table.
Prints PASS/FAIL, exits 1 on failure."""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "workflow", "scripts")
sys.path.insert(0, SCRIPTS)
from build_count_gtf import canonical_masks, subtract, trim_3p  # noqa: E402
from count_flank import fit_side, runs, transcript_windows  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


# --- stop-codon trim (BED12 is stop-inclusive; blocks 0-based half-open) -------
check("trim +", trim_3p([(100, 130)], "+"), [(100, 127)])
check("trim -", trim_3p([(100, 130)], "-"), [(103, 130)])
check("trim + across junction", trim_3p([(100, 130), (200, 202)], "+"), [(100, 129)])
check("trim - across junction", trim_3p([(100, 101), (200, 230)], "-"), [(202, 230)])

# --- subtraction ---------------------------------------------------------------
m = ([10, 50], [[10, 20], [50, 60]])
check("subtract: split in two", subtract([(0, 100)], m), [(0, 10), (20, 50), (60, 100)])
check("subtract: fully covered", subtract([(12, 18)], m), [])
check("subtract: starts inside a mask", subtract([(15, 30)], m), [(20, 30)])
check("subtract: no mask on strand", subtract([(1, 5)], None), [(1, 5)])

# --- canonical flanks: 5' flank upstream of the first CDS, 3' downstream of the last
with tempfile.TemporaryDirectory() as d:
    gtf = os.path.join(d, "a.gtf")
    with open(gtf, "w") as fh:
        fh.write('c\tE\tCDS\t101\t200\t.\t+\t0\ttranscript_id "T1"; gene_id "G1";\n')
        fh.write('c\tE\tCDS\t301\t400\t.\t+\t0\ttranscript_id "T1"; gene_id "G1";\n')
        fh.write('c\tE\tCDS\t1001\t1100\t.\t-\t0\ttranscript_id "T2"; gene_id "G2";\n')
        fh.write('c\tE\texon\t1\t2000\t.\t+\t.\ttranscript_id "T1"; gene_id "G1";\n')
    masks, lines = canonical_masks(gtf, 12, 20)
    check("flank +: 12 up of first, 20 down of last", masks[("c", "+")][1], [[88, 200], [300, 420]])
    check("flank -: 5' is the high end", masks[("c", "-")][1], [[980, 1112]])
    check("only CDS rows copied", len(lines), 3)

# --- windows: d runs utr..-(cds-1) on the start side, transcript-ordered ----------
exons, cds = [(1, 400), (1001, 1400)], [(301, 400), (1001, 1100)]
w = transcript_windows(exons, cds, "+", utr=100, cds_nt=50)
check("window +: start side begins 100 nt upstream", (int(w[0][0]), int(w[0][100])), (201, 301))
check("window +: stop side ends 100 nt past the CDS", (int(w[1][49]), int(w[1][-1])), (1100, 1200))
w = transcript_windows([(1, 400)], [(101, 300)], "-", utr=100, cds_nt=50)
check("window -: start side runs down the genome", (int(w[0][0]), int(w[0][100])), (400, 300))
check("window: too-short UTR rejected", transcript_windows(exons, cds, "+", utr=301, cds_nt=50), None)
check("runs", runs(np.array([5, 3, 4, 10, 11])), [(3, 5), (10, 11)])

# --- flank fit: excess decays over 15 nt, background from d=100 ---------------------
d = np.arange(-59, 301)
prof = pd.Series(np.where(d <= 0, 1.0, 0.02 + 0.98 * np.exp(-d / 3.0)), index=d)
flank, lvl, base, ok = fit_side(prof, 60, (100, 300), 0.05, 100)
check("fit: first d where excess < 5% of CDS-background gap", flank, 9)
check("fit: baseline", round(base, 3), 0.02)
check("fit: floor applies", fit_side(prof, 60, (100, 300), 0.05, 100, min_flank=15)[0], 15)
# A 10-nt block at 30-39 (one transcript's read pile) must not move the flank.
blk = prof.copy(); blk[(blk.index >= 30) & (blk.index < 40)] = 0.15
check("fit: later block ignored", fit_side(blk, 60, (100, 300), 0.05, 100)[0], 9)
# A 5' UTR at 70% of CDS level: step down within 15 nt, then flat near the band.
utr = pd.Series(np.where(d <= 0, 1.0, np.where(d < 15, 0.85, 0.70)), index=d)
utr[(utr.index > 15) & (utr.index % 7 == 0)] += 0.03          # noise that crosses the band
check("fit: high flat UTR -> step, not noise", fit_side(utr, 60, (100, 300), 0.05, 100)[0], 15)
check("fit: never reached -> max, flagged", fit_side(pd.Series(np.ones(len(d)), index=d),
      60, (100, 300), 0.05, 100)[3], False)

# --- counts table: fully subtracted ORFs come back as zero rows -------------------
with tempfile.TemporaryDirectory() as d:
    fc, orfs, out = (os.path.join(d, n) for n in ("fc.txt", "orfs.tsv", "counts.tsv"))
    with open(fc, "w") as fh:
        fh.write("# Program:featureCounts v2.0.8\n")
        fh.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\t/x/align/F1.genome.sorted.bam\n")
        fh.write("G1\tc\t101\t400\t+\t200\t50\n")
        fh.write("orf_1\tc\t500\t530\t+\t30\t4\n")
    pd.DataFrame({"orf_id": ["orf_1", "orf_2"], "aa_length": [10, 12], "cds_nt": [30, 36],
                  "unique_nt": [30, 0]}).to_csv(orfs, sep="\t", index=False)
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "counts_table.py"), "--counts", fc,
                        "--orfs", orfs, "--out", out], capture_output=True, text=True)
    t = pd.read_csv(out, sep="\t").set_index("feature_id")
    check("table: count then RPKM columns", list(t.columns), ["feature_type", "length_nt", "F1", "rpkm_F1"])
    # library = 50 + 4 = 54 reads; G1 = 50 / 0.2 kb / 54e-6
    check("table: RPKM of a gene", round(t.loc["G1", "rpkm_F1"], 1), round(50 / 0.2 / 54e-6, 1))
    check("table: RPKM of a novel ORF uses its unique length", round(t.loc["orf_1", "rpkm_F1"], 1),
          round(4 / 0.03 / 54e-6, 1))
    check("table: subtracted ORF has blank RPKM", bool(pd.isna(t.loc["orf_2", "rpkm_F1"])), True)
    check("table: types", t["feature_type"].to_dict(),
          {"G1": "canonical", "orf_1": "unreviewed_smorf", "orf_2": "unreviewed_smorf"})
    check("table: subtracted ORF is a zero row", (int(t.loc["orf_2", "F1"]), int(t.loc["orf_2", "length_nt"])),
          (0, 0))

sys.exit(1 if failed else 0)
