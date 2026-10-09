#!/usr/bin/env python3
"""Build the featureCounts annotation: canonical CDS plus the unreviewed smORF parts that
no canonical CDS reaches.

Unreviewed = Comprehensive (every catalogue ORF without an exclude_reason) whose
blastoff_class is not sp_identical (Reviewed; tr_identical is counted), and aa_length <=
--max-aa (0 = any length). final_outputs.py reads the counts back: an unreviewed smORF with too low a
pooled RPKM here moves to excluded (min_unique_rpkm).
Each one's CDS comes from its BED12 blocks with the stop codon trimmed, since GTF
CDS rows exclude the stop and the two feature types should be defined alike.

Every canonical transcript's CDS is widened by --flank-5p upstream of its first
CDS row and --flank-3p downstream of its last (strand-aware, from count_flank.py),
and that is subtracted from each smORF on the same strand. What remains is written
as CDS rows with gene_id = transcript_id = orf_id, one row per remaining segment.

Canonical CDS rows are copied unchanged. The subtraction already guarantees that no
smORF row overlaps a canonical CDS row on its strand, so one featureCounts pass over
this GTF counts the same as separate canonical and smORF passes concatenated.

A second table gives each counted smORF's cds_nt and unique_nt (0 when fully subtracted:
those have no GTF row, and counts_table.py adds them back with zero counts).
"""
import argparse
import re
import sys
from bisect import bisect_left
from collections import defaultdict

import pandas as pd

ATTR = re.compile(r'transcript_id "([^"]*)"')


def bed12_blocks(line):
    """BED12 line -> (orf_id, chrom, strand, [(start, end)] 0-based half-open, genomic order)."""
    f = line.rstrip("\n").split("\t")
    chrom, start, name, strand = f[0], int(f[1]), f[3], f[5]
    sizes = [int(x) for x in f[10].rstrip(",").split(",")]
    offs = [int(x) for x in f[11].rstrip(",").split(",")]
    return name, chrom, strand, [(start + o, start + o + s) for o, s in zip(offs, sizes)]


def trim_3p(blocks, strand, n=3):
    """Drop the 3'-most n nt (the stop codon), strand-aware, across block boundaries."""
    blocks = [list(b) for b in sorted(blocks)]
    order = range(len(blocks) - 1, -1, -1) if strand == "+" else range(len(blocks))
    for i in order:
        s, e = blocks[i]
        take = min(n, e - s)
        if strand == "+":
            blocks[i][1] = e - take
        else:
            blocks[i][0] = s + take
        n -= take
        if n == 0:
            break
    return [tuple(b) for b in blocks if b[1] > b[0]]


def canonical_masks(gtf, flank_5p, flank_3p):
    """REF_GTF CDS rows -> ({(chrom, strand): merged 0-based intervals}, CDS lines as text).
    Each transcript's first and last CDS rows are widened by the flanks."""
    rows, lines = defaultdict(list), []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            lines.append(line if line.endswith("\n") else line + "\n")
            m = ATTR.search(f[8])
            rows[(m.group(1) if m else "", f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    mask = defaultdict(list)
    for (tid, chrom, strand), iv in rows.items():
        iv = sorted(iv)
        up, down = (flank_5p, flank_3p) if strand == "+" else (flank_3p, flank_5p)
        iv[0] = (max(0, iv[0][0] - up), iv[0][1])
        iv[-1] = (iv[-1][0], iv[-1][1] + down)
        mask[(chrom, strand)].extend(iv)
    return {k: merge(v) for k, v in mask.items()}, lines


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    starts = [s for s, _ in out]
    return starts, out


def subtract(blocks, mask):
    """Intervals minus a merged (starts, intervals) mask."""
    if mask is None:
        return list(blocks)
    starts, iv = mask
    out = []
    for s, e in blocks:
        cur = s
        i = max(0, bisect_left(starts, s) - 1)
        while i < len(iv) and iv[i][0] < e:
            ms, me = iv[i]
            if me > cur:
                if ms > cur:
                    out.append((cur, min(ms, e)))
                cur = max(cur, me)
            if cur >= e:
                break
            i += 1
        if cur < e:
            out.append((cur, e))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True, help="orf_scores_annotated.tsv (needs blastoff_class)")
    ap.add_argument("--bed12", required=True, help="scoring/orf_candidates.bed12 (stop-inclusive)")
    ap.add_argument("--gtf", required=True, help="reference GTF; its CDS rows are the canonical features")
    ap.add_argument("--flank", help="count_flank.py fit table (flank_5p, flank_3p)")
    ap.add_argument("--flank-5p", type=int, help="override: nt upstream of each canonical CDS")
    ap.add_argument("--flank-3p", type=int, help="override: nt downstream of each canonical CDS")
    ap.add_argument("--max-aa", type=int, default=150, help="longest smORF counted, aa (0 = any)")
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-orfs", required=True, help="orf_id, aa_length, cds_nt, unique_nt")
    ap.add_argument("--exempt-classes", nargs="*", default=["sp_identical"],
                    help="blastoff classes not counted, so exempt from the RPKM gate")
    a = ap.parse_args()

    fit = {}
    if a.flank:
        fit = dict(pd.read_csv(a.flank, sep="\t").itertuples(index=False))
    f5 = a.flank_5p if a.flank_5p is not None else int(fit["flank_5p"])
    f3 = a.flank_3p if a.flank_3p is not None else int(fit["flank_3p"])

    scores = pd.read_csv(a.scores, sep="\t", dtype=str)
    # Reviewed (sp_identical) is the known protein itself: not counted, so exempt from the
    # RPKM gate (join_annotations.EXEMPT_CLASSES; --exempt-classes may add tr_identical)
    keep = ~scores["blastoff_class"].fillna("").isin(set(a.exempt_classes))
    if "exclude_reason" in scores:   # excluded = no longer Comprehensive
        keep &= scores["exclude_reason"].fillna("") == ""
    counted = set(scores.loc[keep, "orf_id"])
    masks, canon_lines = canonical_masks(a.gtf, f5, f3)

    out_rows, table = [], []
    with open(a.bed12) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            orf, chrom, strand, blocks = bed12_blocks(line)
            if orf not in counted:
                continue
            aa = (sum(e - s for s, e in blocks) - 3) // 3
            if a.max_aa and not 0 < aa <= a.max_aa:
                continue
            cds = trim_3p(blocks, strand)
            seg = subtract(cds, masks.get((chrom, strand)))
            table.append((orf, aa, sum(e - s for s, e in cds), sum(e - s for s, e in seg)))
            att = f'gene_id "{orf}"; transcript_id "{orf}";'
            out_rows.extend(f"{chrom}\tBullpen\tCDS\t{s + 1}\t{e}\t.\t{strand}\t.\t{att}\n" for s, e in seg)

    with open(a.out_gtf, "w") as fh:
        fh.writelines(canon_lines)
        fh.writelines(out_rows)
    t = pd.DataFrame(table, columns=["orf_id", "aa_length", "cds_nt", "unique_nt"])
    t.to_csv(a.out_orfs, sep="\t", index=False)
    print(f"build_count_gtf: flanks 5' {f5} / 3' {f3} nt; {len(t)} unreviewed smORFs, "
          f"{(t['unique_nt'] > 0).sum()} with a unique part, {(t['unique_nt'] == 0).sum()} fully subtracted; "
          f"{len(out_rows)} smORF rows + {len(canon_lines)} canonical CDS rows", file=sys.stderr)


if __name__ == "__main__":
    main()
