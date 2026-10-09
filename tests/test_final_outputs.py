"""final_outputs: slim table, summary and README. Prints PASS/FAIL, exits 1 on failure."""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
from final_outputs import assign_sets, exclude_low_counts, flag_above_cds_median, slim, write  # noqa: E402
from ucsc_links import parse_session, ucsc_link  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


full = pd.DataFrame(dict(
    orf_id=["a", "b", "c"], chrom="chr1", start=[0, 10, 20], end=[30, 40, 50], strand="+",
    orf_set=["primary", "comprehensive", "not_supported"], smorf_type=["uORF", "iORF", "iORF"],
    is_smorf=[1, 1, 0], aa_length=[10, 10, 200],
    called_by_price=["True", "False", "True"], called_by_ribocode=["True", "True", "False"],
    score_price=[0.001, np.nan, 0.5], pif=[0.81234, 0.2, 0.9], uniformity=[0.9, 0.9, 0.9],
    dropoff=[0.95, 0.95, np.nan], pass_pif=["True", "False", "True"],
    pass_uniformity=["True", "True", "True"], pass_dropoff=["True", "True", "False"],
    pass_all=["True", "False", "False"], internal_col=1))
o = slim(full)
check("internal columns dropped", "internal_col" in o, False)
check("callers", list(o["callers"]), ["PRICE+RiboCode", "RiboCode", "PRICE"])
check("failed", list(o["failed"]), ["", "pif", "dropoff n/a"])
check("scores rounded", o["pif"][0], 0.812)
check("is_smorf bool", list(o["is_smorf"]), [True, True, False])

# unshared_evidence: drop -> exclude_reason -> orf_set excluded, failed says why
ex = full.drop(columns=["orf_set"]).assign(exclude_reason=["", "no in-frame P-sites outside annotated CDS", np.nan],
                                           unshared_nt=[0, 12, 0], unshared_inframe_psites=[0.0, 0.0, 0.0])
check("excluded set", list(assign_sets(ex)), ["primary", "excluded", "comprehensive"])
check("failed names exclusion", slim(ex)["failed"].tolist()[1],
      "pif, excluded: no in-frame P-sites outside annotated CDS")
check("unshared columns kept", [c for c in ("unshared_nt", "unshared_inframe_psites") if c in slim(ex)],
      ["unshared_nt", "unshared_inframe_psites"])

# Bullpen's catalogue: 0/1 flags and pval_price instead of True/False and score_price.
b = full.drop(columns=["score_price"]).assign(
    called_by_price=[1, 0, 1], called_by_ribocode=[1, 1, 0], pval_price=[0.001, np.nan, 0.5])
ob = slim(b)
check("callers from 0/1 flags", list(ob["callers"]), ["PRICE+RiboCode", "RiboCode", "PRICE"])
check("price_p from pval_price", ob["price_p"].tolist()[0], 0.001)
oi = slim(b.assign(isoform_group=["chr1:+:130", "chr1:+:130", "chr2:-:9"]))
check("isoform_group kept, after strand", list(oi.columns[list(oi.columns).index("strand") + 1:][:1]),
      ["isoform_group"])

plain = parse_session("", "mm39")
sess = parse_session("brmiller/Panda_quad", "mm39")
check("url: BED -> 1-based, padded",
      ucsc_link("chr1", 100, 200, plain, pad=50, fmt="url"),
      "https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39&position=chr1:51-250")
check("bare seqname mapped", ucsc_link("1", 100, 200, plain, pad=0, fmt="url").endswith("chr1:101-200"), True)
check("MT -> chrM", "chrM:" in ucsc_link("MT", 0, 30, plain), True)
check("session short link, no db",
      ucsc_link("chrX", 0, 30, sess, pad=0, fmt="url"),
      "https://genome.ucsc.edu/s/brmiller/Panda_quad?position=chrX:1-30")
check("short link URL accepted", parse_session("https://genome.ucsc.edu/s/brmiller/Panda_quad", "mm39"), sess)
check("hyperlink cell", ucsc_link("chr2", 9, 12, plain, pad=0),
      '=HYPERLINK("https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39&position=chr2:10-12","chr2:10-12")')
check("scaffold -> blank", ucsc_link("GL456233.2", 0, 30, plain), "")
try:
    parse_session("https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=123", "mm39")
    check("live-cart hgsid refused", False, True)
except ValueError:
    check("live-cart hgsid refused", True, True)

thr = pd.DataFrame([dict(pif=0.4, uniformity=0.5, dropoff=0.9, n_cds_used=100, tail="lower",
                         percentile=0.95)])
qc = [dict(check="threshold_pif", observed="0.400", expected="0.511-0.759", status="WARN", note="")]
with tempfile.TemporaryDirectory() as d:
    orfs, summ = write(full, thr, qc, d, "test")
    check("three files", sorted(os.listdir(d)), ["README.md", "orfs.csv", "summary.csv"])
    row = summ[summ["item"] == "all ORFs"].iloc[0]
    check("counts per set", [int(row[s]) for s in ("primary", "comprehensive", "not_supported", "all")],
          [1, 1, 1, 3])
    check("counts are integers", str(summ["primary"].dtype), "Int64")
    text = open(os.path.join(d, "README.md")).read()
    check("README lists the WARN", "WARN threshold_pif" in text, True)
    check("README documents orf_set", "| `orf_set` |" in text, True)
    back = pd.read_csv(os.path.join(d, "orfs.csv"))
    check("csv round-trips", list(back["orf_id"]), ["a", "b", "c"])

# Provided mode: no scores, no orf_set, no thresholds or QC.
prov = pd.DataFrame(dict(orf_id=["o1", "o2"], source_id=["myORF_7,tA", "tB"], chrom="chr1",
                         start=[10, 100], end=[45, 109], strand=["+", "-"], smorf_type=["uORF", "novel"],
                         is_smorf=[1, 1], aa_length=[3, 2]))
with tempfile.TemporaryDirectory() as d:
    orfs, summ = write(prov, None, [], d, "prov")
    check("provided: source_id kept, after orf_id", list(orfs.columns[:2]), ["orf_id", "source_id"])
    check("provided: no orf_set, callers or scores",
          [c for c in ("orf_set", "callers", "pif", "failed") if c in orfs], [])
    check("provided: summary has no per-set columns", list(summ.columns), ["section", "item", "all", "note"])
    check("provided: no threshold rows", "threshold" in set(summ["section"]), False)
    text = open(os.path.join(d, "README.md")).read()
    check("provided: README says unscored", "not scored for translation" in text, True)
    check("provided: README headline", "2 ORFs, 2 of them smORFs" in text, True)
    check("provided: no Primary in README", "Primary |" in text, False)

# count gate: unreviewed smORFs need pooled RPKM > 5 in their unique part
# library total = 1e6 reads over both samples, so RPKM = reads / kb
f = pd.DataFrame(dict(orf_id=["lo", "hi", "rev", "old", "zero", "edge"], pass_all=[True] * 6,
                      exclude_reason=["", "", "", "psORF: earlier", "", ""]))
cnt = pd.DataFrame(dict(feature_id=["G1", "lo", "hi", "old", "zero", "edge"],
                        feature_type=["canonical"] + ["unreviewed_smorf"] * 5,
                        length_nt=[900, 1000, 100, 30, 0, 1000],
                        S1=[499_990, 2, 3, 0, 0, 5], S2=[499_998, 2, 0, 0, 0, 0],
                        rpkm_S1=[1.0, 99.0, 1.0, 0, np.nan, 1.0], rpkm_S2=[1.0] * 6))
check("count gate: moved", exclude_low_counts(f, cnt, 5.0), 3)
check("count gate: unique_reads summed", f["unique_reads"].tolist()[:2], [4, 3])
check("count gate: pooled rpkm, not the rpkm_ columns", [round(x, 3) for x in f["unique_rpkm"][:2]],
      [4.0, 30.0])
check("count gate: pass_all", f["pass_all"].tolist(), [False, True, True, False, False, False])
check("count gate: reason", f.at[0, "exclude_reason"], "pooled RPKM 5 or less outside flanked annotated CDS")
check("count gate: RPKM exactly 5 fails (strictly above)", f.at[5, "exclude_reason"] != "", True)
check("count gate: not counted (reviewed) untouched", (pd.isna(f.at[2, "unique_reads"]), f.at[2, "exclude_reason"]),
      (True, ""))
check("count gate: earlier reason kept", f.at[3, "exclude_reason"], "psORF: earlier")
check("count gate: fully subtracted excluded", assign_sets(f).tolist(),
      ["excluded", "primary", "primary", "excluded", "excluded", "excluded"])

# above_cds_median: CDS median density 1 per codon per million, 1e6 in-frame CDS P-sites -> lambda = codons
fl = pd.DataFrame(dict(
    orf_set=["primary"] * 4 + ["comprehensive", "primary"], is_smorf=[1, 1, 1, 1, 1, 0],
    blastoff_class=["sp_identical", "novel", "novel", "novel", "novel", "novel"],
    body_codons=[10, 10, 10, 10, 10, 10], body_inframe_psites_unique=[100.0, 1.0, 1000.0, 300.0, 500.0, 500.0],
    free_body_codons=[0, 10, 0, 5, 10, 10], free_body_inframe_psites_unique=[0.0, 1.0, 0.0, 200.0, 500.0, 500.0]))
thr1 = pd.Series(dict(cds_median_density_unique=1.0, inframe_cds_psites_unique=1e6))
check("cds median: flagged", flag_above_cds_median(fl, thr1, 0.05), 2)
check("cds median: Reviewed own codons / sparse / overlap with no free codons / free codons / not tested",
      [None if pd.isna(x) else bool(x) for x in fl["above_cds_median"]], [True, False, False, True, None, None])
check("cds median: ratio from the codons used (NaN -> -1)",
      fl["cds_median_density_ratio"].fillna(-1).tolist()[:4], [10.0, 0.1, -1.0, 40.0])
check("cds median: no free codons -> no FDR", bool(np.isnan(fl.at[2, "cds_median_fdr"])), True)
check("cds median: kept in orfs.csv", {"above_cds_median", "cds_median_density_ratio"} <= set(slim(fl).columns), True)
check("cds median: no thresholds -> skipped", flag_above_cds_median(fl.copy(), None, 0.05), None)

print(f"\n{'ALL PASS' if not failed else f'{failed} FAILED'}")
sys.exit(1 if failed else 0)
