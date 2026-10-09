#!/usr/bin/env python3
"""Write the three files a reader opens: orfs.csv, summary.csv, README.md.

  orfs.csv     one row per ORF, the columns below. Filter orf_set == primary
               for the Primary set. The full table stays in details/. CSV so
               the ucsc column's =HYPERLINK(...) cells click through in Excel.
  summary.csv  counts per set and ORF type, the fitted thresholds, and the
               pass_all QC against Chothani et al.
  README.md    what the files and columns mean, with the headline numbers.

Input is the full per-ORF table: orf_master.tsv here, or Bullpen's
orf_scores_annotated.tsv joined to its catalogue. Columns it lacks are left out
of orfs.csv, not invented. Provided mode (no --thresholds/--qc): the table is
Bullpen's orfs_annotated.tsv, with no scores and no orf_set.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

from ucsc_links import parse_session, ucsc_link

# (output column, meaning). Built by slim() below; order is the file's order.
SMORF_MAX_AA = 150   # set from --smorf-max-aa


def smorf_label():
    return f"smORFs (<= {SMORF_MAX_AA} aa)"


COLUMNS = [
    ("orf_id", "ORF identifier"),
    ("source_id", "the ORF's id in the provided GTF (its orf_id attribute, else transcript_id); "
                  "identical structures listed together, comma-separated"),
    ("ucsc", "UCSC genome-browser link, opens on the ORF (blank for scaffolds)"),
    ("chrom", "chromosome"),
    ("start", "start, 0-based (BED)"),
    ("end", "end, half-open (BED)"),
    ("strand", "strand"),
    ("isoform_group", "the ORF's 3' end after the stop codon, chrom:strand:position (BED): ORFs "
                      "with the same value share the stop, i.e. are start or splice variants of "
                      "one ORF. Each distinct structure is its own row (catalogue_merge: exact)"),
    ("host_gene_name", "gene the ORF sits in"),
    ("gene_id", "Ensembl gene id"),
    ("orf_set", "primary: caller-supported and passes all three scores (Chothani Primary set). "
                "comprehensive: caller-supported, fails a score. excluded: too few in-frame P-sites "
                "in the sequence it doesn't share with an annotated CDS (unshared_min_inframe_psites), "
                "no sequence outside annotated CDS unless its protein is exactly a Swiss-Prot "
                "sequence (Reviewed; require_unshared_sequence), an unreviewed smORF whose "
                "pooled RPKM there is not above min_unique_rpkm "
                "(unique_rpkm), an excluded type such as Trunc-Iso (exclude_types), or an excluded "
                "chromosome such as chrM (exclude_chroms); reason in failed. "
                "not_supported: in neither set"),
    ("smorf_type", "ORF type: the Annotator's call, relabeled from the UniProt match "
                   "(annotated_CDS; Trunc-Iso: shorter, no residues of its own; N-Iso / C-Iso / "
                   "NC-Iso: own residues at that end; Ins-Iso / Del-Iso: same ends, an extra / "
                   "skipped internal segment; TrEMBL: identical to a TrEMBL protein). "
                   "annotated_CDS / TrEMBL only when the protein is exactly that sequence; "
                   "near-identical and other hits keep the Annotator's type"),
    ("general_smorf_type", "smorf_type grouped by transcript architecture: Upstream (uORF, uoORF, "
                           "uaORF, uaoORF, udORF), Internal (iORF, riORF), Isoform (*-Iso), "
                           "Downstream (dORF, doORF, daORF, daoORF), Non-coding (lncRNA, "
                           "intergenic, psORF, aORF, eORF, UA, not_annotated_*), Reviewed "
                           "(annotated_CDS), TrEMBL"),
    ("is_smorf", "True if {max_aa} aa or shorter"),
    ("aa_length", "protein length, aa"),
    ("callers", "which ORF callers found it: PRICE, RiboCode, or both"),
    ("price_p", "PRICE p-value (blank if PRICE did not call it)"),
    ("ppm", "P-sites per million, scaled so the column sums to 1e6"),
    ("pif", "proportion of P-sites in frame"),
    ("uniformity", "share of codons whose P-sites are mostly in frame"),
    ("dropoff", "in-frame P-sites before vs after the stop codon (blank if not measurable)"),
    ("unshared_nt", "nt of the ORF that no annotated CDS (any Ensembl transcript) reads in the "
                    "same frame; sequence a CDS covers in another frame is the ORF's own "
                    "(unshared_frame: aware)"),
    ("unshared_psites", "pooled P-sites in those nt"),
    ("unshared_inframe_psites", "pooled P-sites in those nt, in the ORF's frame: its own evidence"),
    ("psites_unique", "pooled P-sites from uniquely mapped reads (STAR MAPQ 255) in the ORF"),
    ("ppm_unique", "ppm from uniquely mapped reads only (column sums to 1e6)"),
    ("unshared_inframe_psites_unique", "unshared_inframe_psites from uniquely mapped reads only; "
                                       "a pseudogene's parent-gene reads drop out"),
    ("unique_reads", "featureCounts reads in the part no flanked annotated CDS reaches (any frame), "
                     "summed over libraries (unreviewed smORFs only; blank if not counted)"),
    ("unique_rpkm", "pooled RPKM of those reads: unique_reads per kb of that part per million reads "
                    "assigned in counts.tsv, all libraries pooled (blank if not counted or no such part)"),
    ("above_cds_median", "Primary smORFs only: True if its in-frame P-site density from uniquely mapped "
                         "reads, over the codons between start and stop, is above this dataset's median "
                         "annotated CDS (one-sided Poisson, BH FDR < cds_median_fdr). Reviewed and ORFs "
                         "no annotated CDS touches use all those codons; other ORFs only the codons no "
                         "annotated CDS touches in any frame (none: False). Blank if not a Primary smORF"),
    ("cds_median_density_ratio", "that density over the CDS median (1 = typical annotated CDS)"),
    ("cds_median_fdr", "BH FDR of the above-CDS-median test"),
    ("failed", "scores below threshold; 'dropoff n/a' = too few reads at the stop to measure; "
               "'excluded: ...' = why the ORF is out of Comprehensive and Primary"),
    ("blastoff_class", "UniProt match: sp_/tr_ + identical (the protein is exactly that sequence), "
                       "near_identical (>= blastoff_min_pident end to end, not exact), N_truncated, "
                       "C_truncated, internal, N_/C_/NC_extended, insertion, deletion, partial, "
                       "homolog; or novel"),
    ("blastoff_acc", "matched UniProt accession"),
    ("blastoff_gene", "matched UniProt gene"),
    ("shortstop_class", "ShortStop class: sam_secreted, sam_intracellular, prisms; or "
                        "not_scored_length / not_scored_non_AUG / not_scored_sequence. "
                        "Model is human-trained"),
    ("shortstop_probability", "ShortStop's probability for its class (blank if not scored)"),
    ("diamond_n_species", "species whose Ensembl transcripts encode the protein: best forward-frame "
                          "DIAMOND blastx hit covering >= diamond_min_cov of it at tblastn E <= "
                          "diamond_evalue, against the whole transcriptome or the host gene's orthologs"),
    ("diamond_species", "those species, in diamond_species order; per-species call, gene, identity, "
                        "coverage and E-values are in the full table"),
    ("diamond_n_untestable", "species where even a perfect match could not reach diamond_evalue "
                             "(too short for the search space): absence there is not evidence"),
    ("aa_seq", "protein sequence"),
]
SETS = ["primary", "comprehensive", "excluded", "not_supported"]
GENERAL_ORDER = ["Upstream", "Internal", "Isoform", "Downstream", "Non-coding", "Reviewed", "TrEMBL"]


def _true(v):
    """Flags arrive as True/False (orf_master) or 1/0 (Bullpen's catalogue)."""
    return str(v) in ("True", "1", "1.0")


def _callers(row):
    c = [n for n, col in (("PRICE", "called_by_price"), ("RiboCode", "called_by_ribocode"))
         if _true(row.get(col))]
    return "+".join(c)


def _failed(row):
    out = [s for s in ("pif", "uniformity") if str(row.get(f"pass_{s}")) == "False"]
    if pd.isna(row.get("dropoff")):
        out.append("dropoff n/a")
    elif str(row.get("pass_dropoff")) == "False":
        out.append("dropoff")
    if str(row.get("exclude_reason", "")) not in ("", "nan"):
        out.append(f"excluded: {row['exclude_reason']}")
    return ", ".join(out)


def assign_sets(full):
    """Bullpen's catalogue is already caller-filtered: every row is Comprehensive, Primary
    if it passes, and excluded if scoring gave it an exclude_reason."""
    sets = full["pass_all"].astype(str).map({"True": "primary"}).fillna("comprehensive")
    if "exclude_reason" in full:
        sets[full["exclude_reason"].fillna("").astype(str) != ""] = "excluded"
    return sets


def exclude_low_counts(full, counts, min_rpkm):
    """Unreviewed smORFs need a pooled RPKM above min_rpkm in the part no flanked canonical
    CDS reaches (counts.tsv): reads summed over libraries, per kb of that part, per million
    reads assigned in the table over all libraries. The rest, and counted smORFs with no
    such part, move to excluded; an earlier reason is kept. Adds unique_reads and
    unique_rpkm. Returns how many were moved."""
    sm = counts[counts["feature_type"] != "canonical"]
    samples = [c for c in counts.columns
               if c not in ("feature_id", "feature_type", "length_nt") and not c.startswith("rpkm_")]
    reads = sm.set_index("feature_id")[samples].sum(axis=1)
    kb = sm.set_index("feature_id")["length_nt"].where(lambda n: n > 0) / 1e3
    lib_m = counts[samples].to_numpy().sum() / 1e6
    full["unique_reads"] = full["orf_id"].map(reads)
    full["unique_rpkm"] = full["orf_id"].map(reads / kb / lib_m if lib_m else reads * float("nan"))
    counted = full["unique_reads"].notna()
    if "exclude_reason" not in full:
        full["exclude_reason"] = ""
    blank = full["exclude_reason"].fillna("").astype(str) == ""
    hit = counted & ~(full["unique_rpkm"] > min_rpkm)      # NaN (no unique part) fails
    full.loc[hit & blank, "exclude_reason"] = (
        f"pooled RPKM {min_rpkm:g} or less outside flanked annotated CDS")
    full.loc[hit, "pass_all"] = False
    return int((hit & blank).sum())


def flag_above_cds_median(full, thr, fdr):
    """Adds above_cds_median / cds_median_density_ratio / cds_median_fdr to Primary smORFs (see
    COLUMNS). The density is in-frame unique P-sites per body codon per million in-frame unique
    CDS P-sites (score_orfs.cds_median_density). Returns how many are flagged, or None if the
    table or thresholds lack the inputs (provided mode, older runs)."""
    need = ("body_codons", "body_inframe_psites_unique", "free_body_codons", "free_body_inframe_psites_unique")
    if thr is None or "cds_median_density_unique" not in thr or "orf_set" not in full \
            or any(c not in full for c in need):
        return None
    from scipy import stats
    med, tot = float(thr["cds_median_density_unique"]), float(thr["inframe_cds_psites_unique"])
    sm = full["is_smorf"].astype(str).isin(["1", "True"]) if "is_smorf" in full else True
    test = (full["orf_set"] == "primary") & sm
    rev = full["blastoff_class"].astype(str).eq("sp_identical") if "blastoff_class" in full else False
    own = rev | (full["free_body_codons"] == full["body_codons"])
    k = np.where(own, full["body_inframe_psites_unique"], full["free_body_inframe_psites_unique"]).astype(float)
    n = np.where(own, full["body_codons"], full["free_body_codons"]).astype(float)
    ok = test.to_numpy() & (n > 0)
    p = np.full(len(full), np.nan)
    p[ok] = stats.poisson.sf(k[ok] - 1, med * n[ok] * tot / 1e6)
    q = np.full(len(full), np.nan)
    if ok.any():
        q[ok] = stats.false_discovery_control(p[ok], method="bh")
    flag = pd.Series(pd.NA, index=full.index, dtype="boolean")
    flag[test] = False
    flag[ok & (q < fdr)] = True
    full["above_cds_median"] = flag
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = k / n / tot * 1e6 / med
    full["cds_median_density_ratio"] = np.where(ok, np.round(ratio, 3), np.nan)
    full["cds_median_fdr"] = q
    return int(flag.sum())


def slim(full):
    df = pd.DataFrame(index=full.index)
    for col, _ in COLUMNS:
        if col == "callers" and "called_by_price" in full:
            df[col] = full.apply(_callers, axis=1)
        elif col == "price_p" and ("score_price" in full or "pval_price" in full):
            df[col] = full["score_price" if "score_price" in full else "pval_price"]
        elif col == "failed" and "pass_pif" in full:
            df[col] = full.apply(_failed, axis=1)
        elif col == "is_smorf" and col in full:
            df[col] = full[col].astype(str).isin(["1", "True"])
        elif col in ("pif", "uniformity", "dropoff", "shortstop_probability") and col in full:
            df[col] = full[col].round(3)
        elif col == "ppm" and col in full:
            df[col] = full[col].round(4)
        elif col in full:
            df[col] = full[col]
    return df


def summary(orfs, thresholds, qc):
    """thresholds None (provided mode): counts only, and no per-set columns."""
    rows = []
    sets = SETS if "orf_set" in orfs else []

    def counts(section, item, sub):
        c = sub["orf_set"].value_counts() if sets else pd.Series(dtype=int)
        rows.append(dict(section=section, item=item, **{s: int(c.get(s, 0)) for s in sets},
                         all=len(sub), note=""))

    counts("orfs", "all ORFs", orfs)
    if "is_smorf" in orfs:
        counts("orfs", smorf_label(), orfs[orfs["is_smorf"]])
    if "above_cds_median" in orfs:
        counts("orfs", "smORFs above CDS median", orfs[orfs["above_cds_median"].astype(str) == "True"])
    if "smorf_type" in orfs:
        order = orfs.loc[orfs["orf_set"] == "primary", "smorf_type"].value_counts().index \
            if "orf_set" in orfs else orfs["smorf_type"].value_counts().index
        rest = [t for t in orfs["smorf_type"].dropna().unique() if t not in set(order)]
        for t in list(order) + sorted(rest):
            counts("type (all ORFs)", t, orfs[orfs["smorf_type"] == t])
        if "is_smorf" in orfs:
            sm = orfs[orfs["is_smorf"]]
            for t in list(order) + sorted(rest):
                if (sm["smorf_type"] == t).any():
                    counts("type (smORFs)", t, sm[sm["smorf_type"] == t])
    if "general_smorf_type" in orfs and "is_smorf" in orfs:
        sm = orfs[orfs["is_smorf"]]
        for g in GENERAL_ORDER:
            if (sm["general_smorf_type"] == g).any():
                counts("general type (smORFs)", g, sm[sm["general_smorf_type"] == g])
    thr = thresholds.iloc[0] if thresholds is not None else None
    for s in ("pif", "uniformity", "dropoff") if thr is not None else ():
        how = "mean of CDS" if s == "dropoff" else f"{thr['tail']} tail, {thr['percentile']}"
        rows.append(dict(section="threshold", item=s,
                         note=f"{float(thr[s]):.3f}: fitted on {int(thr['n_cds_used']):,} "
                              f"annotated CDS ({how})"))
    if thr is not None and "cds_median_density_unique" in thr:
        rows.append(dict(section="threshold", item="cds_median_density",
                         note=f"{float(thr['cds_median_density_unique']):.4g}: median in-frame unique P-sites "
                              f"per body codon per million in-frame unique CDS P-sites, "
                              f"{int(thr['n_cds_density']):,} annotated CDS"))
    for r in qc:
        rows.append(dict(section="qc vs Chothani / GENCODE v2 2025", item=r["check"], note=
                         f"{r['status']}: {r['observed']} (Chothani / GENCODE v2 2025: {r['expected']}) {r['note']}".strip()))
    out = pd.DataFrame(rows, columns=["section", "item", *sets, "all", "note"])
    return out.astype({c: "Int64" for c in (*sets, "all")})


def _file_rows(details):
    """README file-table rows from 'name: description' (one string or a list)."""
    return [f"| {d.split(': ', 1)[0]} | {d.split(': ', 1)[1]} |"
            for d in ([details] if isinstance(details, str) else details)]


def readme(name, orfs, summ, qc, details="`details/`: the full per-ORF table (`orfs_full.tsv`) and every intermediate file"):
    n = summ[(summ.section == "orfs")].set_index("item")
    def c(item, s):
        return f"{int(n.loc[item, s]):,}" if item in n.index else "?"
    if "orf_set" not in orfs:
        return readme_provided(name, orfs, n, details)
    warns = [r for r in qc if r["status"] == "WARN"]
    lines = [
        f"# {name}", "",
        "| File | What it is |", "|---|---|",
        "| `orfs.csv` | one row per ORF. Filter `orf_set == primary` for the Primary set |",
        "| `summary.csv` | counts per set and ORF type, fitted thresholds, QC vs Chothani / GENCODE v2 2025 |",
        *_file_rows(details), "",
        "## Headline", "",
        "| | Primary | Comprehensive | All ORFs |", "|---|---:|---:|---:|",
        f"| ORFs | {c('all ORFs', 'primary')} | "
        f"{int(n.loc['all ORFs', 'primary'] + n.loc['all ORFs', 'comprehensive']):,} | {c('all ORFs', 'all')} |",
    ]
    if smorf_label() in n.index:
        s = smorf_label()
        lines.append(f"| {s} | {c(s, 'primary')} | "
                     f"{int(n.loc[s, 'primary'] + n.loc[s, 'comprehensive']):,} | {c(s, 'all')} |")
    lines += ["", "Comprehensive = caller-supported (RiboCode, or PRICE at its default FDR 0.1). "
              "Primary = Comprehensive and passes PIF, Uniformity and Drop-off, each thresholded "
              "on this dataset's annotated CDS (Chothani et al.).", ""]
    if "smORFs above CDS median" in n.index:
        lines += [f"Above CDS median: {c('smORFs above CDS median', 'primary')} Primary smORFs whose "
                  "in-frame ribosome density (uniquely mapped reads, codons between start and stop) is "
                  "significantly above this dataset's median annotated CDS (`above_cds_median`).", ""]
    lines += ["## QC vs Chothani / GENCODE v2 2025", "",
              f"{len(qc) - len(warns)} checks OK, {len(warns)} WARN (details in `summary.csv`)."]
    lines += [f"- WARN {r['check']}: {r['observed']} (Chothani / GENCODE v2 2025: {r['expected']})" for r in warns]
    lines += ["", "## Columns of orfs.csv", "", "| Column | Meaning |", "|---|---|"]
    lines += [f"| `{col}` | {desc.replace('{max_aa}', str(SMORF_MAX_AA))} |" for col, desc in COLUMNS if col in orfs]
    return "\n".join(lines) + "\n"


def readme_provided(name, orfs, n, details):
    head = f"{int(n.loc['all ORFs', 'all']):,} ORFs"
    if smorf_label() in n.index:
        head += f", {int(n.loc[smorf_label(), 'all']):,} of them {smorf_label()}"
    lines = [
        f"# {name}", "",
        "| File | What it is |", "|---|---|",
        "| `orfs.csv` | one row per ORF in the provided GTF |",
        "| `summary.csv` | counts per ORF type |",
        *_file_rows(details), "",
        "## Headline", "",
        head + ".", "",
        "Provided mode: the ORFs come from a GTF, not from Ribo-seq. They were not called, "
        "not scored for translation (no PIF, Uniformity or Drop-off) and not collapsed by "
        "protein similarity, so they are neither Primary nor Comprehensive.", "",
        "## Columns of orfs.csv", "", "| Column | Meaning |", "|---|---|"]
    lines += [f"| `{col}` | {desc.replace('{max_aa}', str(SMORF_MAX_AA))} |" for col, desc in COLUMNS if col in orfs]
    return "\n".join(lines) + "\n"


def write(full, thresholds, qc, outdir, name, details=None):
    orfs = slim(full)
    summ = summary(orfs, thresholds, qc)
    orfs.to_csv(os.path.join(outdir, "orfs.csv"), index=False)
    summ.to_csv(os.path.join(outdir, "summary.csv"), index=False)
    with open(os.path.join(outdir, "README.md"), "w") as fh:
        fh.write(readme(name, orfs, summ, qc, *([details] if details else [])))
    return orfs, summ


def read_aa_fasta(path):
    """orf_id -> protein. Headers are orf_id(+) / orf_id(-) (bedtools -s -nameOnly)."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = re.sub(r"\([+-]\)$", "", line[1:].split()[0]), []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="full per-ORF table")
    ap.add_argument("--catalogue", help="catalogue TSV to fill coordinates / callers the table lacks")
    ap.add_argument("--thresholds", help="scoring/thresholds.tsv; omit in provided mode")
    ap.add_argument("--qc", help="pass_all_qc.tsv; omit in provided mode")
    ap.add_argument("--aa-fasta", help="catalogue protein FASTA; fills aa_seq by orf_id")
    ap.add_argument("--ucsc-genome", help="UCSC assembly (e.g. mm39); fills the ucsc column "
                    "when the table lacks one")
    ap.add_argument("--ucsc-session", default="",
                    help="saved session 'user/name' or https://genome.ucsc.edu/s/<user>/<name>; "
                         "must be on the same assembly")
    ap.add_argument("--ucsc-pad", type=int, default=50, metavar="BP")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--name", default="Ribo-seq ORFs")
    ap.add_argument("--details", nargs="+",
                    help="README rows for the other files, each 'name: description'")
    ap.add_argument("--smorf-max-aa", type=int, default=150, help="is_smorf cutoff, for the labels")
    ap.add_argument("--counts", help="counts.tsv; unreviewed smORFs at or below "
                                     "--min-unique-rpkm move to excluded")
    ap.add_argument("--min-unique-rpkm", type=float, default=5.0)
    ap.add_argument("--cds-median-fdr", type=float, default=0.05,
                    help="BH FDR for above_cds_median (Primary smORFs vs the CDS median density)")
    a = ap.parse_args(argv)
    full = pd.read_csv(a.table, sep="\t", low_memory=False)
    if a.catalogue:
        cat = pd.read_csv(a.catalogue, sep="\t", low_memory=False)
        extra = [c for c in cat.columns if c not in full.columns]
        full = full.merge(cat[["orf_id", *extra]], on="orf_id", how="left")
    if "stop_key" in full and "isoform_group" not in full:
        full["isoform_group"] = full["stop_key"]
    if a.aa_fasta:
        aa = read_aa_fasta(a.aa_fasta)
        full["aa_seq"] = full["orf_id"].map(aa)
        missing = int(full["aa_seq"].isna().sum())
        if missing:
            raise SystemExit(f"final_outputs: {missing} ORFs have no protein in {a.aa_fasta}")
    if "ucsc" not in full and (a.ucsc_genome or a.ucsc_session):
        session = parse_session(a.ucsc_session, a.ucsc_genome)
        full["ucsc"] = [ucsc_link(c, s, e, session, a.ucsc_pad)
                        for c, s, e in zip(full["chrom"], full["start"], full["end"])]
    if a.counts and "pass_all" in full:
        k = exclude_low_counts(full, pd.read_csv(a.counts, sep="\t"), a.min_unique_rpkm)
        print(f"final_outputs: {k:,} unreviewed smORFs excluded for pooled RPKM "
              f"{a.min_unique_rpkm:g} or less", file=sys.stderr)
    if "orf_set" not in full and "pass_all" in full:
        full["orf_set"] = assign_sets(full)
    global SMORF_MAX_AA
    SMORF_MAX_AA = a.smorf_max_aa
    if bool(a.thresholds) != bool(a.qc):
        raise SystemExit("final_outputs: give both --thresholds and --qc, or neither (provided mode)")
    qc = pd.read_csv(a.qc, sep="\t").fillna("").to_dict("records") if a.qc else []
    thr = pd.read_csv(a.thresholds, sep="\t") if a.thresholds else None
    k = flag_above_cds_median(full, thr.iloc[0] if thr is not None else None, a.cds_median_fdr)
    if k is not None:
        print(f"final_outputs: {k:,} Primary smORFs above the CDS median density "
              f"(FDR < {a.cds_median_fdr:g})", file=sys.stderr)
    write(full, thr, qc, a.outdir, a.name, a.details)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
