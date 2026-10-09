#!/usr/bin/env python3
"""Provided mode: an ORF GTF -> the catalogue files build_catalogue.py writes.

One ORF per transcript_id with CDS rows. The GTF needs `transcript` and `CDS`
features at minimum; `stop_codon` rows, when present, are folded into the ORF.
The BED12 is stop-inclusive: an ORF one codon short of a stop (Ensembl's CDS
convention) is extended 3 nt, so normalize_orf_bed12.py sees one convention.

source_id is the transcript's `orf_id` attribute, else its transcript_id.
Identical block structures are written once, their source_ids joined with ','.
Nothing is merged by overlap or protein similarity: no MMseqs2 collapse.

Aborts, listing offenders, on: CDS rows with no transcript row, a transcript on
two chromosomes or strands, seqnames not in the genome, a CDS not a multiple of
3, or an ORF that neither ends in nor is followed by a stop.
"""
from __future__ import annotations

import argparse
import collections
import sys

from pyfaidx import Fasta

from build_catalogue import ATTR, COMP, STOPS, SeqnameMap, last_codon

COLS = ["orf_id", "source_id", "chrom", "start", "end", "strand", "gene_id", "transcript_id",
        "orf_class", "aa_length", "is_smorf"]


def die(what, bad):
    shown = ", ".join(map(str, bad[:10])) + (f" ... (+{len(bad) - 10})" if len(bad) > 10 else "")
    sys.exit(f"provided_catalogue: {len(bad):,} {what}: {shown}")


def read_gtf(path):
    """-> {transcript_id: dict(chrom, strand, gene_id, source_id, blocks)} and the CDS-only ids."""
    tx, cds, bad_loc = {}, collections.defaultdict(list), set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "CDS", "stop_codon"):
                continue
            a = dict(ATTR.findall(f[8]))
            tid = a.get("transcript_id")
            if not tid:
                continue
            if f[2] == "transcript":
                tx[tid] = dict(chrom=f[0], strand=f[6], gene_id=a.get("gene_id", ""),
                               source_id=a.get("orf_id", tid))
            else:
                cds[tid].append((f[0], f[6], int(f[3]) - 1, int(f[4])))
    orphans = sorted(t for t in cds if t not in tx)
    for tid, rows in cds.items():
        if tid in tx and {(c, s) for c, s, _, _ in rows} != {(tx[tid]["chrom"], tx[tid]["strand"])}:
            bad_loc.add(tid)
    for tid in list(tx):
        if tid not in cds:
            del tx[tid]   # non-coding transcript: nothing to annotate
        else:
            tx[tid]["blocks"] = merge_blocks([(s, e) for _, _, s, e in cds[tid]])
    return tx, orphans, sorted(bad_loc)


def merge_blocks(blocks):
    out = []
    for s, e in sorted(blocks):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def next_codon(genome, chrom, strand, blocks):
    if strand == "+":
        e = blocks[-1][1]
        return str(genome[chrom][e:e + 3]).upper()
    s = blocks[0][0]
    return str(genome[chrom][max(0, s - 3):s]).upper().translate(COMP)[::-1] if s >= 3 else ""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, help="provided ORF GTF (transcript + CDS rows)")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--smorf-max-aa", type=int, default=150)
    ap.add_argument("--out-prefix", required=True, help="writes <prefix>.bed12/.tsv/.consensus.bed12")
    ap.add_argument("--report", required=True)
    a = ap.parse_args(argv)

    tx, orphans, bad_loc = read_gtf(a.gtf)
    if orphans:
        die("transcript_ids have CDS rows but no transcript row", orphans)
    if bad_loc:
        die("transcripts have CDS rows off their transcript's chromosome or strand", bad_loc)
    if not tx:
        sys.exit(f"provided_catalogue: no transcript with CDS rows in {a.gtf}")

    genome = Fasta(a.genome, as_raw=True, sequence_always_upper=True)
    seqmap = SeqnameMap(genome.keys())
    for t in tx.values():
        t["chrom"] = seqmap(t["chrom"]) or t["chrom"]
    if seqmap.unknown:
        die("seqnames are not in the genome", sorted(seqmap.unknown))

    stats = collections.Counter(transcripts_with_cds=len(tx))
    stats.update({f"seqname_renamed:{k}": v for k, v in seqmap.renamed.items()})
    no_stop, bad_len = [], []
    for tid, t in tx.items():
        if last_codon(genome, t["chrom"], t["strand"], t["blocks"]) in STOPS:
            stats["ends_in_stop"] += 1
        elif next_codon(genome, t["chrom"], t["strand"], t["blocks"]) in STOPS:
            b = t["blocks"]
            t["blocks"] = (b[:-1] + [(b[-1][0], b[-1][1] + 3)] if t["strand"] == "+"
                           else [(b[0][0] - 3, b[0][1])] + b[1:])
            stats["extended_3nt_to_stop"] += 1
        else:
            no_stop.append(tid)
        if sum(e - s for s, e in t["blocks"]) % 3:
            bad_len.append(tid)
    if bad_len:
        die("ORFs are not a multiple of 3 nt (stop included)", sorted(bad_len))
    if no_stop:
        die("ORFs neither end in nor are followed by a stop codon", sorted(no_stop))

    by_struct = collections.defaultdict(list)
    for tid, t in tx.items():
        by_struct[(t["chrom"], t["strand"], tuple(t["blocks"]))].append(tid)
    stats["duplicate_structures_folded"] = len(tx) - len(by_struct)

    keys = sorted(by_struct, key=lambda k: (k[0], k[2][0][0], k[2][-1][1], k[1],
                                            min(tx[t]["source_id"] for t in by_struct[k])))
    with open(a.out_prefix + ".bed12", "w") as bed, open(a.out_prefix + ".tsv", "w") as tsv:
        tsv.write("\t".join(COLS) + "\n")
        for i, key in enumerate(keys, 1):
            chrom, strand, bl = key
            tids = sorted(by_struct[key], key=lambda t: tx[t]["source_id"])
            oid, start, end = f"orf_{i:08d}", bl[0][0], bl[-1][1]
            bed.write("\t".join(map(str, [
                chrom, start, end, oid, 0, strand, start, end, 0, len(bl),
                ",".join(str(e - s) for s, e in bl) + ",",
                ",".join(str(s - start) for s, e in bl) + ","])) + "\n")
            aa = sum(e - s for s, e in bl) // 3 - 1
            genes = ",".join(dict.fromkeys(tx[t]["gene_id"] for t in tids))
            tsv.write("\t".join(map(str, [
                oid, ",".join(tx[t]["source_id"] for t in tids), chrom, start, end, strand, genes,
                ",".join(tids), "provided", aa, int(0 < aa <= a.smorf_max_aa)])) + "\n")
    with open(a.out_prefix + ".bed12") as src, open(a.out_prefix + ".consensus.bed12", "w") as dst:
        dst.write(src.read())
    stats["orfs"] = len(keys)
    with open(a.report, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in stats.items():
            fh.write(f"{k}\t{v}\n")
            print(f"  {k}: {v:,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
