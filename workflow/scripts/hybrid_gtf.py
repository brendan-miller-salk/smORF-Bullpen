#!/usr/bin/env python3
"""Build the two hybrid GTFs.

  hybrid       = backbone + novel transcripts of --ribocode-class-codes, a
                 synthetic `gene` row before each gene's first transcript row
                 (concat.awk). Feeds RiboCode and the hybrid STAR index.
  full_hybrid  = the whole filtered reference + every line of backbone + novel
                 transcripts of --class-codes whose transcript_id (or, for gene
                 rows, gene_id) the reference lacks (build_full_hybrid_gtf.awk).
                 Feeds GEDI/PRICE and the catalogue.

The RiboCode set must be a subset of the full set (the catalogue projects
RiboCode's transcripts through full_hybrid). It can be smaller because novel
transcripts that repeat a reference intron chain (k, m, n) make transcriptome
reads multi-map, and RiboCode's metaplots drop multi-mapped reads: on the u,i,j,k,o,m,n
set it found no periodic read length in libraries that had 50-75k calls before.

The backbone is reference.canonical.gtf (one transcript per gene) by
default, or the full filtered reference for `ribocode_annotation: full`.
Novel transcripts come from gffcompare's annotated GTF (filter.awk): unstranded
lines dropped, a transcript kept if any of its lines has a kept class_code.
With no --annotated, hybrid is the backbone unchanged and full_hybrid the
reference.
"""
import argparse
import collections
import re
import sys

TID = re.compile(r'transcript_id "([^"]*)"')
GID = re.compile(r'gene_id "([^"]*)"')
CODE = re.compile(r'class_code "([^"]*)"')


def data_lines(path):
    with open(path) as fh:
        return [l if l.endswith("\n") else l + "\n" for l in fh if not l.startswith("#")]


def filter_novel(path, codes):
    stored = collections.OrderedDict()
    code_of = {}
    for line in data_lines(path):
        f = line.split("\t")
        if len(f) < 9 or f[6] == ".":
            continue
        t = TID.search(f[8])
        if not t:
            continue
        stored.setdefault(t.group(1), []).append(line)
        c = CODE.search(f[8])
        if c:
            code_of.setdefault(t.group(1), set()).add(c.group(1))
    counts, kept_counts, out = collections.Counter(), collections.Counter(), []
    for tid, lines in stored.items():
        cs = code_of.get(tid, {"?"})
        for c in cs:
            counts[c] += 1
        if cs & codes:
            out.extend(lines)
            for c in cs & codes:
                kept_counts[c] += 1
    return out, sum(1 for t in stored if code_of.get(t, set()) & codes), counts, kept_counts


def concat(lines):
    """concat.awk: synthetic gene row before each gene's first transcript row."""
    gene_seen = {GID.search(l.split("\t")[8]).group(1)
                 for l in lines if l.split("\t")[2] == "gene" and GID.search(l.split("\t")[8])}
    span = {}
    for l in lines:
        f = l.split("\t")
        if f[2] != "transcript":
            continue
        g = GID.search(f[8])
        if not g or g.group(1) in gene_seen:
            continue
        gid = g.group(1)
        s, e = int(f[3]), int(f[4])
        if gid in span:
            chrom, src, lo, hi, strand = span[gid]
            span[gid] = (chrom, src, min(lo, s), max(hi, e), strand)
        else:
            span[gid] = (f[0], f[1], s, e, f[6])
    out, done = [], set()
    for l in lines:
        f = l.split("\t")
        if f[2] == "transcript":
            g = GID.search(f[8])
            if g and g.group(1) in span and g.group(1) not in done:
                chrom, src, lo, hi, strand = span[g.group(1)]
                out.append("\t".join([chrom, src, "gene", str(lo), str(hi), ".", strand, ".",
                                      f'gene_id "{g.group(1)}";']) + "\n")
                done.add(g.group(1))
        out.append(l)
    return out, len(done)


def full_hybrid(ref_lines, hybrid_lines):
    """build_full_hybrid_gtf.awk."""
    ref_tx, ref_genes = set(), set()
    for l in ref_lines:
        a = l.split("\t")[8]
        t, g = TID.search(a), GID.search(a)
        if t:
            ref_tx.add(t.group(1))
        if g:
            ref_genes.add(g.group(1))
    out = list(ref_lines)
    for l in hybrid_lines:
        a = l.split("\t")[8]
        t = TID.search(a)
        if t:
            if t.group(1) not in ref_tx:
                out.append(l)
        else:
            g = GID.search(a)
            if g and g.group(1) not in ref_genes:
                out.append(l)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, help="filtered reference GTF (gtf_filter.py)")
    ap.add_argument("--backbone", required=True, help="reference.canonical.gtf, or the reference itself")
    ap.add_argument("--annotated", help="gffcompare <prefix>.annotated.gtf; omit for no novel models")
    ap.add_argument("--class-codes", default="u,i,j", help="novel classes in full_hybrid (PRICE, catalogue)")
    ap.add_argument("--ribocode-class-codes", help="novel classes in hybrid (RiboCode); default --class-codes")
    ap.add_argument("--out-hybrid", required=True)
    ap.add_argument("--out-full-hybrid", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()

    ref, backbone = data_lines(a.reference), data_lines(a.backbone)
    codes = set(a.class_codes.split(","))
    ribo_codes = set((a.ribocode_class_codes or a.class_codes).split(","))
    if not ribo_codes <= codes:
        sys.exit(f"hybrid_gtf: RiboCode classes {sorted(ribo_codes - codes)} are not in --class-codes")
    if a.annotated:
        novel, n_kept, counts, kept = filter_novel(a.annotated, codes)
        novel_ribo, n_ribo, _, kept_ribo = filter_novel(a.annotated, ribo_codes)
        hybrid, n_gene_rows = concat(backbone + novel_ribo)
        full = full_hybrid(ref, concat(backbone + novel)[0])
    else:
        counts = kept = kept_ribo = {}
        n_kept = n_ribo = n_gene_rows = 0
        hybrid = backbone
        full = full_hybrid(ref, hybrid)
    with open(a.report, "w") as rep:
        rep.write("class_code\tn_transcripts\tkept\tkept_ribocode\n")
        for c in sorted(counts):
            rep.write(f"{c}\t{counts[c]}\t{kept.get(c, 0)}\t{kept_ribo.get(c, 0)}\n")
    with open(a.out_hybrid, "w") as fh:
        fh.writelines(hybrid)
    with open(a.out_full_hybrid, "w") as fh:
        fh.writelines(full)
    print(f"hybrid_gtf: {n_kept:,} novel transcripts in full hybrid, {n_ribo:,} in hybrid; "
          f"hybrid {len(hybrid):,} lines "
          f"({n_gene_rows:,} synthetic gene rows); full hybrid {len(full):,} lines", file=sys.stderr)


if __name__ == "__main__":
    main()
