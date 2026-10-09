#!/usr/bin/env python3
"""Classify ORF proteins against UniProt: known protein, isoform, homolog, novel.

Replaces PANDA's blastoff.py, which emitted only the IDs with any bitscore >= 40
hit. It never split Swiss-Prot / TrEMBL / novel, and never looked at coverage,
so it could not see truncations.

  query     orf_candidates.bed12 + genome -> query_short.fa (<= --short-aa) and
            query_long.fa. Short queries go to `blastp -task blastp-short`.
  classify  blastp outfmt 6 (BLAST_COLS below) + the UniProt FASTA -> one row
            per ORF, orf_id first.

A hit is one ORF-protein pair. BLAST may split an alignment across an inserted or
skipped segment (long queries run with several HSPs); HSPs in the same order on
both sequences are chained into one hit first. Identity is pident, or, when the
alignment has gaps, identity over the ungapped columns: gaps are judged as
insertions / deletions by length, not as mismatches.

Classes (<db> = sp | tr).
  <db>_identical     EXACT: the ORF's protein (literal first residue, stop trimmed)
                     equals the UniProt sequence, by direct string comparison, not
                     from BLAST. Swiss-Prot canonical > Swiss-Prot isoform > TrEMBL.
A hit is near-identical when its identity >= --min-pident and it covers the query to
within --slack aa at both ends; every class below except homolog also needs an
alignment of >= --min-shared-aa aa, so a 5-9 aa query can't be called truncated or
internal on a chance 4-6 aa match (it is exact, a homolog, or novel):
  <db>_near_identical near-identical, covers the subject end to end, lengths agree
                     within --slack, but not the exact sequence (paralog copy, SNP,
                     1-2 aa different start or stop). Keeps the Annotator's type
  <db>_insertion     the same, but the ORF is longer: an extra exon or inserted
                     segment inside the protein. offset_aa = extra aa
  <db>_deletion      the same, but the ORF is shorter: a skipped exon or deleted
                     segment. offset_aa = missing aa
  <db>_N_truncated   reaches the subject's C-terminus but starts downstream of its
                     M. That is a downstream start. offset_aa = aa missing at the N-terminus
  <db>_C_truncated   starts at the subject's M, stops early. offset_aa = aa
                     missing at the C-terminus. Check against the PRICE/RiboCode stop convention
  <db>_internal      both subject ends missing
  <db>_N_extended    covers the subject end to end, with extra residues before
                     its M (an upstream start). offset_aa = extra aa
  <db>_C_extended    the same, with extra residues after its stop (read-through or
                     another last exon). offset_aa = extra aa
  <db>_NC_extended   extra at both ends. offset_aa = total extra aa.
                     None of the extended classes is a short isoform
  <db>_partial       shares an in-frame stretch of >= --min-shared-aa at
                     >= --min-pident, but the rest of the ORF differs (another
                     exon, frame or retained intron at one end). offset_aa = where
                     the shared stretch starts in the subject
  <db>_homolog       any other hit with evalue <= --homolog-evalue
  novel              none of the above
Best hit per ORF: exact identity (either db) first, then near-identical / truncated /
internal, then extended / insertion / deletion, then partial, then homolog; within a
tier, near-identical before any truncation, then Swiss-Prot before TrEMBL, then
N > C > internal, then canonical before isoform, then bitscore.
`identical_to` lists every accession whose sequence the ORF's protein equals, from
either db, Swiss-Prot canonical first (blastoff_acc is the first).
`pass_bit40` is PANDA's criterion: any hit with bitscore >= 40 and evalue <= 1e-4.

`diverges` says which end of the ORF differs from the matched protein: the
end(s) where the ORF has residues of its own (an upstream start, read-through,
another exon or frame). For a pure truncation, it is the end(s) it is missing.
N, C, NC, or blank (identical, or not in frame).

relabel_smorf_type() is the Annotator relabel. Only an exact match gets a database
label; near-identical, homologs, non-identical TrEMBL matches and novel ORFs keep their
type:
  any type                 Swiss-Prot identical -> annotated_CDS
                           TrEMBL identical -> TrEMBL
  iORF uoORF doORF         Swiss-Prot truncated / internal -> Trunc-Iso (no residues of
  uaoORF daoORF uaORF      its own: neither upstream of the start nor past the stop)
  daORF                    Swiss-Prot extended / partial -> N-Iso / C-Iso / NC-Iso, by
                           the end(s) where the ORF has its own residues (`diverges`)
                           Swiss-Prot insertion -> Ins-Iso, deletion -> Del-Iso
join_annotations.py imports it when both annotator and blastoff are enabled.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict

BLAST_COLS = ("qseqid sseqid pident length mismatch gapopen qstart qend qlen "
              "sstart send slen evalue bitscore").split()
OUT_COLS = ["orf_id", "blastoff_class", "blastoff_db", "blastoff_acc", "blastoff_entry",
            "blastoff_gene", "blastoff_is_isoform", "blastoff_fragment", "blastoff_offset_aa",
            "blastoff_diverges",
            "blastoff_pident", "blastoff_qcov", "blastoff_scov", "blastoff_evalue",
            "blastoff_bitscore", "blastoff_identical_to", "blastoff_pass_bit40"]
NEAR_ORDER = {"near_identical": 0, "N_truncated": 1, "C_truncated": 2, "internal": 3}

# Annotator types every ORF that shares a CDS start, or starts inside a CDS, as iORF. It
# does not know isoforms, so full-length proteins and truncated isoforms land
# there with true internal ORFs.
TRUNCATED = {"N_truncated", "C_truncated", "internal"}
EXTENDED = {"N_extended", "C_extended", "NC_extended"}
INDEL = {"insertion": "Ins-Iso", "deletion": "Del-Iso"}
IN_FRAME = {"identical", "near_identical", "partial"} | TRUNCATED | EXTENDED | set(INDEL)
CHAIN_OVERLAP_AA = 10   # chained HSPs may overlap this much at a splice boundary
# Overlapping and isoform-ambiguous types: an in-frame match makes them an N- or
# C-terminal isoform of the known protein.
ISO_TYPES = {"uoORF", "doORF", "uaoORF", "daoORF", "uaORF", "daORF"}
ISO_LABEL = {"N": "N-Iso", "C": "C-Iso", "NC": "NC-Iso"}


def relabel_smorf_type(smorf_type, blastoff_class, diverges=""):
    """Annotator type -> annotated_CDS / Trunc-Iso / N-Iso / C-Iso / NC-Iso / Ins-Iso /
    Del-Iso / TrEMBL for UniProt matches (rules in the module docstring); else unchanged."""
    if not isinstance(blastoff_class, str) or "_" not in blastoff_class:
        return smorf_type
    db, kind = blastoff_class.split("_", 1)
    if kind == "identical":
        return "TrEMBL" if db == "tr" else "annotated_CDS"
    if kind == "near_identical":   # a paralog copy or a 1-2 aa different end: not the protein
        return smorf_type
    if smorf_type not in ISO_TYPES | {"iORF"} or kind not in IN_FRAME or db == "tr":
        return smorf_type
    if kind in TRUNCATED:     # lies wholly inside the protein
        return "Trunc-Iso"
    if kind in INDEL:
        return INDEL[kind]
    if not isinstance(diverges, str) or diverges not in ISO_LABEL:
        diverges = kind.split("_")[0] if kind in EXTENDED else "NC"
    return ISO_LABEL[diverges]


def diverging_ends(h, slack):
    """The end(s) where the ORF carries residues the protein lacks. For a pure
    truncation (none of its own), the end(s) where it is missing the protein's."""
    n = h["qstart"] > 1 + slack
    c = h["qend"] < h["qlen"] - slack
    if not (n or c):
        n = h["sstart"] > 1 + slack
        c = h["send"] < h["slen"] - slack
    return ("N" if n else "") + ("C" if c else "")


# --- query -----------------------------------------------------------------
COMP = str.maketrans("ACGTN", "TGCAN")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def cmd_query(a):
    from pyfaidx import Fasta
    genome = Fasta(a.genome, as_raw=True, sequence_always_upper=True)
    n = Counter()
    with open(a.bed12) as bed, open(a.out_short, "w") as fs, open(a.out_long, "w") as fl:
        for line in bed:
            f = line.rstrip("\n").split("\t")
            chrom, start, strand, oid = f[0], int(f[1]), f[5], f[3]
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            starts = [int(x) for x in f[11].rstrip(",").split(",")]
            nt = "".join(genome[chrom][start + s:start + s + z] for s, z in zip(starts, sizes))
            if strand == "-":
                nt = nt.translate(COMP)[::-1]
            aa = translate(nt).rstrip("*")
            if not aa or "*" in aa:
                n["skipped"] += 1   # internal stop: not a clean ORF, classify as novel
                continue
            (fs if len(aa) <= a.short_aa else fl).write(f">{oid}\n{aa}\n")
            n["short" if len(aa) <= a.short_aa else "long"] += 1
    print(f"queries: {n['short']:,} <= {a.short_aa} aa, {n['long']:,} longer, "
          f"{n['skipped']:,} skipped (internal stop)", file=sys.stderr)
    return 0


# --- classify --------------------------------------------------------------
HDR = re.compile(r"^>(sp|tr)\|([^|]+)\|(\S+)\s*(.*)$")


def read_uniprot(path, seqs=None):
    """accession -> (db, entry, gene, is_isoform, is_fragment). With a dict `seqs`,
    also fills accession -> sequence."""
    meta, acc, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                if acc is not None and seqs is not None:
                    parts.append(line.strip())
                continue
            if acc is not None and seqs is not None:
                seqs[acc] = "".join(parts)
            acc, parts = None, []
            m = HDR.match(line.rstrip("\n"))
            if not m:
                continue
            db, acc, entry, desc = m.groups()
            g = re.search(r"\bGN=(\S+)", desc)
            iso = "-" in acc or desc.startswith("Isoform of")
            meta[acc] = (db, entry, g.group(1) if g else "", iso, "(Fragment)" in desc)
    if acc is not None and seqs is not None:
        seqs[acc] = "".join(parts)
    return meta


def read_fasta(paths):
    out, k = {}, None
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    k = line[1:].split()[0]
                    out[k] = []
                elif k is not None:
                    out[k].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def accession(sseqid):
    # With -parse_seqids BLAST may print sp|ACC|ENTRY or bare ACC; handle both.
    return sseqid.split("|")[1] if "|" in sseqid else sseqid


def chain(hsps):
    """One hit from one ORF-protein pair's HSPs: the best HSP plus every other HSP that
    lies wholly before or after all chained ones on both sequences. `pid` is the
    identity classification uses: over ungapped columns when there are gaps."""
    kept = []
    for h in sorted(hsps, key=lambda x: -x["bitscore"]):
        if all((h["qend"] <= c["qstart"] + CHAIN_OVERLAP_AA and h["send"] <= c["sstart"] + CHAIN_OVERLAP_AA)
               or (h["qstart"] >= c["qend"] - CHAIN_OVERLAP_AA and h["sstart"] >= c["send"] - CHAIN_OVERLAP_AA)
               for c in kept):
            kept.append(h)
    ident = sum(h["pident"] * h["length"] / 100 for h in kept)
    length = sum(h["length"] for h in kept)
    mismatch = sum(h["mismatch"] for h in kept)
    gaps = sum(h["gapopen"] for h in kept) + len(kept) - 1
    m = dict(kept[0])
    m.update(qstart=min(h["qstart"] for h in kept), qend=max(h["qend"] for h in kept),
             sstart=min(h["sstart"] for h in kept), send=max(h["send"] for h in kept),
             length=length, mismatch=mismatch, gapopen=gaps, pident=100 * ident / length,
             evalue=min(h["evalue"] for h in kept), bitscore=sum(h["bitscore"] for h in kept),
             hsps=len(kept))
    m["pid"] = 100 * ident / (ident + mismatch) if gaps and ident + mismatch else m["pident"]
    return m


def call(h, a):
    """(kind, offset_aa) for one (chained) hit, or None if it is not a reportable hit.
    Exact identity is decided in cmd_classify by sequence, not here."""
    q_full = h["qstart"] <= 1 + a.slack and h["qend"] >= h["qlen"] - a.slack
    s_n = h["sstart"] <= 1 + a.slack
    s_c = h["send"] >= h["slen"] - a.slack
    pid = h.get("pid", h["pident"])
    long_enough = h["length"] >= a.min_shared_aa
    if pid >= a.min_pident and q_full and long_enough:
        if s_n and s_c:
            extra = (h["qend"] - h["qstart"]) - (h["send"] - h["sstart"])
            if extra > a.slack:
                return "insertion", extra
            if extra < -a.slack:
                return "deletion", -extra
            return "near_identical", 0
        if s_c:
            return "N_truncated", h["sstart"] - 1
        if s_n:
            return "C_truncated", h["slen"] - h["send"]
        return "internal", h["sstart"] - 1
    if pid >= a.min_pident and s_n and s_c and long_enough:
        n_ext, c_ext = h["qstart"] - 1, h["qlen"] - h["qend"]
        if n_ext > a.slack and c_ext > a.slack:
            return "NC_extended", n_ext + c_ext
        return ("N_extended", n_ext) if n_ext > a.slack else ("C_extended", c_ext)
    if pid >= a.min_pident and long_enough:
        return "partial", h["sstart"] - 1
    if h["evalue"] <= a.homolog_evalue:
        return "homolog", ""
    return None


def cmd_classify(a):
    useq = {}
    meta = read_uniprot(a.uniprot_fasta, useq)
    by_seq = defaultdict(list)
    for acc, s in useq.items():
        by_seq[s].append(acc)
    qseq = read_fasta(a.queries)
    ids = []
    for path in a.queries:
        with open(path) as fh:
            ids += [l[1:].split()[0] for l in fh if l.startswith(">")]
    if a.bed12:   # ORFs the query step skipped still get a row
        with open(a.bed12) as fh:
            ids = [l.split("\t")[3] for l in fh if l.strip()]
    best, identical_to, bit40 = {}, defaultdict(list), set()
    unknown = 0
    pairs = defaultdict(list)   # (orf, accession) -> HSPs, chained below
    for path in a.blast:
        with open(path) as fh:
            for line in fh:
                v = line.rstrip("\n").split("\t")
                h = dict(zip(BLAST_COLS, v))
                for k in ("length", "mismatch", "gapopen", "qstart", "qend", "qlen", "sstart", "send", "slen"):
                    h[k] = int(h[k])
                for k in ("pident", "evalue", "bitscore"):
                    h[k] = float(h[k])
                acc = accession(h["sseqid"])
                if acc not in meta:
                    unknown += 1
                    continue
                if h["bitscore"] >= 40 and h["evalue"] <= 1e-4:
                    bit40.add(h["qseqid"])
                pairs[(h["qseqid"], acc)].append(h)
    for (q, acc), hsps in pairs.items():
        h = chain(hsps)
        db, entry, gene, iso, frag = meta[acc]
        c = call(h, a)
        if c is None:
            continue
        kind, off = c
        tier = (0 if kind in NEAR_ORDER else 1 if kind in EXTENDED or kind in INDEL
                else 2 if kind == "partial" else 3)
        key = (tier, kind != "near_identical", db != "sp", NEAR_ORDER.get(kind, 9),
               iso, -h["bitscore"])
        if q not in best or key < best[q][0]:
            best[q] = (key, kind, off, acc, h)
    # exact identity by sequence, ahead of any BLAST call: Swiss-Prot canonical first
    for q, s in qseq.items():
        accs = [acc for acc in by_seq.get(s, []) if acc in meta]
        if not accs:
            continue
        accs.sort(key=lambda acc: (meta[acc][0] != "sp", meta[acc][3], meta[acc][4], acc))
        identical_to[q] = accs
        L = len(s)
        h = dict(pident=100.0, length=L, mismatch=0, gapopen=0, qstart=1, qend=L, qlen=L,
                 sstart=1, send=L, slen=L, evalue=0.0, bitscore=0.0)
        best[q] = ((-1,), "identical", 0, accs[0], h)
    if a.queries and not qseq:
        sys.exit("blastoff classify: --queries gave no sequences, so exact identity can't be tested")
    if unknown:
        print(f"WARNING: {unknown:,} hits to accessions absent from {a.uniprot_fasta}",
              file=sys.stderr)

    n = Counter()
    with open(a.out, "w") as out:
        out.write("\t".join(OUT_COLS) + "\n")
        for oid in ids:
            if oid not in best:
                row = [oid, "novel"] + [""] * (len(OUT_COLS) - 3) + [str(oid in bit40)]
                n["novel"] += 1
            else:
                _, kind, off, acc, h = best[oid]
                db, entry, gene, iso, frag = meta[acc]
                cls = f"{db}_{kind}"
                n[cls] += 1
                div = diverging_ends(h, a.slack) if kind in IN_FRAME else ""
                row = [oid, cls, db, acc, entry, gene, str(iso), str(frag), str(off), div,
                       f"{h['pident']:.1f}",
                       f"{(h['qend'] - h['qstart'] + 1) / h['qlen']:.3f}",
                       f"{(h['send'] - h['sstart'] + 1) / h['slen']:.3f}",
                       f"{h['evalue']:.3g}", f"{h['bitscore']:.1f}",
                       ",".join(identical_to.get(oid, [])), str(oid in bit40)]
            out.write("\t".join(row) + "\n")
    total = sum(n.values())
    print(f"classified {total:,} ORFs; pass_bit40 {len(bit40 & set(ids)):,}", file=sys.stderr)
    for cls, k in sorted(n.items(), key=lambda x: -x[1]):
        print(f"  {cls:<18} {k:>8,}  ({100 * k / total:.1f}%)", file=sys.stderr)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("query", help="bed12 + genome -> short/long query FASTA")
    q.add_argument("--bed12", required=True)
    q.add_argument("--genome", required=True)
    q.add_argument("--short-aa", type=int, default=30)
    q.add_argument("--out-short", required=True)
    q.add_argument("--out-long", required=True)
    c = sub.add_parser("classify", help="blastp outfmt 6 -> per-ORF UniProt class")
    c.add_argument("--blast", nargs="+", required=True,
                   help="outfmt '6 " + " ".join(BLAST_COLS) + "'")
    c.add_argument("--uniprot-fasta", required=True, help="the FASTA the BLAST db was built from")
    c.add_argument("--queries", nargs="*", default=[], help="query FASTAs (ORF ids)")
    c.add_argument("--bed12", help="ORF universe; overrides --queries so skipped ORFs get a row")
    c.add_argument("--min-pident", type=float, default=95.0)
    c.add_argument("--slack", type=int, default=2, help="aa tolerated at each alignment end")
    c.add_argument("--min-shared-aa", type=int, default=10,
                   help="partial: shortest in-frame shared stretch (guards chance short-query hits)")
    c.add_argument("--homolog-evalue", type=float, default=1e-4)
    c.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    return cmd_query(a) if a.cmd == "query" else cmd_classify(a)


if __name__ == "__main__":
    raise SystemExit(main())
