#!/usr/bin/env python3
"""Reviewed / unreviewed smORFs: GTF + nucleotide + peptide FASTA of the Comprehensive
smORFs, split by their UniProt match (blastoff_class).

Both sets: orf_set primary or comprehensive (provided mode, no orf_set: every ORF)
and is_smorf.
    reviewed_smorfs     blastoff_class sp_identical  (protein = a Swiss-Prot sequence)
    unreviewed_smorfs   everything else: novel, any tr_* (tr_identical included), and
                        non-identical sp_* (near_identical, isoforms, partial, homolog, ...)

Per set, in --outdir:
    <set>.gtf  per ORF: transcript, exon (stop included), CDS (stop excluded, with
                  phase), start_codon, stop_codon; gene_id = transcript_id = orf_id.
                  Rows of a split CDS / start / stop codon are in transcription order,
                  each with its frame (Ensembl: bases before the next codon starts)
    <set>.nuc  CDS nucleotides from the genome, stop codon included
    <set>.pep  orfs.csv aa_seq (literal first codon, stop trimmed); must equal the
                  translation of .nuc
FASTA records follow orfs.csv order; GTF records are grouped per ORF, ORFs sorted by
chrom and start. BED12 blocks give the structure (0-based half-open,
stop codon included); the GTF is 1-based inclusive.
"""
import argparse
import os
import sys

import pandas as pd

COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
STOPS = {"TAA", "TAG", "TGA"}


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def _true(v):
    return str(v) in ("True", "1", "1.0")


SETS = {
    "reviewed_smorfs": lambda bc: bc == "sp_identical",
    "unreviewed_smorfs": lambda bc: bc != "sp_identical",
}


def select(orfs, name):
    keep = orfs["is_smorf"].map(_true)
    if "orf_set" in orfs:
        keep &= orfs["orf_set"].isin(["primary", "comprehensive"])
    keep &= SETS[name](orfs["blastoff_class"].fillna("").astype(str))
    return orfs[keep]


def read_bed12(path, wanted):
    """orf_id -> (chrom, strand, [(start, end)] genomic 0-based half-open, ascending)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[3] not in wanted:
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            starts = [int(x) for x in f[11].rstrip(",").split(",")]
            out[f[3]] = (f[0], f[5], [(start + s, start + s + z) for s, z in zip(starts, sizes)])
    return out


def positions(blocks, strand):
    """Genomic 0-based positions in transcription order."""
    pos = [p for s, e in blocks for p in range(s, e)]
    return pos[::-1] if strand == "-" else pos


def runs(pos):
    """Contiguous runs of a position set -> [(start0, end0)] ascending."""
    pos = sorted(pos)
    out = []
    for p in pos:
        if out and p == out[-1][1]:
            out[-1][1] = p + 1
        else:
            out.append([p, p + 1])
    return [tuple(r) for r in out]


def framed(segs, strand, feat):
    """Rows for one feature's segments, in transcription order, each with its frame."""
    segs = segs[::-1] if strand == "-" else segs
    out, done = [], 0
    for s, e in segs:
        out.append((s, e, feat, str((3 - done % 3) % 3)))
        done += e - s
    return out


def gtf_rows(oid, chrom, strand, blocks, attrs):
    pos = positions(blocks, strand)
    cds, start, stop = pos[:-3], pos[:3], pos[-3:]
    att = f'gene_id "{oid}"; transcript_id "{oid}";' + "".join(
        f' {k} "{v}";' for k, v in attrs.items() if v not in ("", "nan", None))
    lo, hi = blocks[0][0], blocks[-1][1]
    rows = [(chrom, lo, hi, "transcript", ".")]
    rows += [(chrom, s, e, "exon", ".") for s, e in blocks]
    for segs, feat in ((runs(cds), "CDS"), (runs(start), "start_codon"), (runs(stop), "stop_codon")):
        rows += [(chrom, s, e, f, ph) for s, e, f, ph in framed(segs, strand, feat)]
    return ["\t".join([c, "Bullpen", feat, str(s + 1), str(e), ".", strand, ph, att])
            for c, s, e, feat, ph in rows]


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_set(sel, bed, genome, prefix):
    """-> number of ORFs whose genome translation disagrees with aa_seq."""
    gtf, n_bad = [], 0   # (chrom, start, rows)
    with open(f"{prefix}.nuc", "w") as fn, open(f"{prefix}.pep", "w") as fp:
        for r in sel.itertuples(index=False):
            chrom, strand, blocks = bed[r.orf_id]
            nt = "".join(genome[chrom][s:e] for s, e in blocks)
            if strand == "-":
                nt = nt.translate(COMP)[::-1]
            aa = r.aa_seq
            if translate(nt).rstrip("*") != aa or nt[-3:] not in STOPS:
                n_bad += 1
                print(f"  MISMATCH {r.orf_id}: translation {translate(nt)[:20]}... vs aa_seq "
                      f"{aa[:20]}..., last codon {nt[-3:]}", file=sys.stderr)
            stype, host = getattr(r, "smorf_type", ""), getattr(r, "host_gene_name", "")
            desc = (f"{r.orf_id} {chrom}:{blocks[0][0] + 1}-{blocks[-1][1]}({strand}) "
                    f"smorf_type={stype or '-'} host_gene={host or '-'} "
                    f"blastoff_class={r.blastoff_class} "
                    + (f"orf_set={r.orf_set} " if "orf_set" in sel else "")
                    + f"aa={len(aa)}")
            fn.write(f">{desc}\n{wrap(nt)}\n")
            fp.write(f">{desc}\n{wrap(aa)}\n")
            gtf.append((chrom, blocks[0][0], gtf_rows(r.orf_id, chrom, strand, blocks, {
                "gene_name": host, "smorf_type": stype,
                "blastoff_class": r.blastoff_class, "orf_set": getattr(r, "orf_set", ""),
                "isoform_group": getattr(r, "isoform_group", "")})))
    gtf.sort(key=lambda x: (x[0], x[1]))
    with open(f"{prefix}.gtf", "w") as fh:
        for _, _, rows in gtf:
            fh.write("\n".join(rows) + "\n")
    return n_bad


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orfs", required=True, help="orfs.csv")
    p.add_argument("--bed12", required=True, help="scoring/orf_candidates.bed12")
    p.add_argument("--genome", required=True)
    p.add_argument("--outdir", required=True, help="writes <set>.gtf/.nuc/.pep for both sets")
    a = p.parse_args(argv)

    from pyfaidx import Fasta
    orfs = pd.read_csv(a.orfs, dtype=str, keep_default_na=False)
    sets = {name: select(orfs, name) for name in SETS}
    wanted = set().union(*(set(s["orf_id"]) for s in sets.values()))
    bed = read_bed12(a.bed12, wanted)
    missing = sorted(wanted - set(bed))
    if missing:
        sys.exit(f"smorf_sets: {len(missing)} ORFs not in {a.bed12}, e.g. {missing[:3]}")
    genome = Fasta(a.genome, as_raw=True, sequence_always_upper=True)

    n_bad = 0
    for name, sel in sets.items():
        n_bad += write_set(sel, bed, genome, os.path.join(a.outdir, name))
        bc = sel["blastoff_class"]
        sp = bc.str.startswith("sp_")
        parts = {"Swiss-Prot identical": (bc == "sp_identical").sum(),
                 "Swiss-Prot non-identical": (sp & (bc != "sp_identical")).sum(),
                 "TrEMBL": bc.str.startswith("tr_").sum(), "novel": (bc == "novel").sum()}
        print(f"{name}: {len(sel):,} (" + ", ".join(f"{v:,} {k}" for k, v in parts.items() if v)
              + ")", file=sys.stderr)
    if n_bad:
        sys.exit(f"smorf_sets: {n_bad} ORFs whose genome translation disagrees with aa_seq")
    return 0


if __name__ == "__main__":
    sys.exit(main())
