#!/usr/bin/env python3
"""PRICE + RiboCode calls -> one cohort ORF catalogue (BED12 + TSV).

Three details of the merge:

1. Seqnames are harmonised BEFORE merging. GEDI/PRICE strips 'chr' and emits
   bare Ensembl names ('1', 'M'); RiboCode follows the GTF ('chr1'). Merged
   raw, no PRICE call could ever cluster with a RiboCode call.
2. aa_length is taken from the block span. RiboCode's ORF_length column
   excludes the stop while tstart..tstop includes it; subtracting 3 from
   ORF_length would undercount every RiboCode ORF by one residue.
3. The stop-codon convention is measured per caller against the genome and
   reported, not assumed. Scoring's normalize step then enforces it.

--merge exact (default) keeps every distinct ORF structure (chrom, strand, exon blocks
incl. the stop) as its own row, as Chothani et al. do ("all isoforms retained"): calls
of the same structure, from either caller or any transcript, become one row credited to
each caller that made it. ORFs that share a stop codon share `stop_key`
(final_outputs.py turns it into isoform_group). --merge overlap folds by overlap: calls
at >= --min-reciprocal-overlap with the same junctions become one row, the
representative chosen by class then length, so a start variant replaces another.

Coordinates: BED12 0-based half-open. PRICE Location blocks are already
0-based half-open; RiboCode ORF_tstart/tstop are 1-based inclusive on the
transcript and are projected onto that transcript's exons from the hybrid GTF.
"""
import argparse
import collections
import hashlib
import math
import os
import re
import sys

from pyfaidx import Fasta

ATTR = re.compile(r'(\S+) "([^"]*)"')
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CLASS_MAP = {
    "price": {"cds": "canonical_cds", "ext": "canonical_cds", "trunc": "canonical_cds",
              "variant": "canonical_cds", "uoorf": "uoORF", "uorf": "uORF", "dorf": "dORF",
              "iorf": "intORF", "ncrna": "novel_u", "orphan": "novel_u", "intronic": "other"},
    "ribocode": {"annotated": "canonical_cds", "uorf": "uORF", "overlap_uorf": "uoORF",
                 "dorf": "dORF", "overlap_dorf": "doORF", "internal": "intORF",
                 "novel": "novel_u"},
}
CLASS_RANK = {c: i for i, c in enumerate(
    ["canonical_cds", "uoORF", "uORF", "doORF", "dORF", "intORF", "novel_u", "other"])}
CALLERS = ["price", "ribocode"]


def die(msg):
    sys.exit(f"build_catalogue: {msg}")


# --------------------------------------------------------------------------- inputs

def duplicate_inputs(paths):
    """Groups of non-empty files with byte-identical content. Two samples'
    caller outputs should never match; when they do, the caller reused one
    sample's intermediate for another (RiboCode's CWD psites.hd5 cache did this).
    Empty files are skipped: several samples may legitimately have no calls."""
    by_hash = collections.defaultdict(list)
    for path in paths:
        if os.path.getsize(path) == 0:
            continue
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        by_hash[h.hexdigest()].append(path)
    return [sorted(g) for g in by_hash.values() if len(g) > 1]


def load_exons(gtf):
    """transcript_id -> (chrom, strand, [(start0, end), ...] in genomic order), plus gene_id."""
    ex = collections.defaultdict(list)
    meta = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon":
                continue
            at = dict(ATTR.findall(f[8]))
            tid = at["transcript_id"]
            ex[tid].append((int(f[3]) - 1, int(f[4])))
            meta[tid] = (f[0], f[6], at.get("gene_id", ""))
    return {t: (meta[t][0], meta[t][1], sorted(b), meta[t][2]) for t, b in ex.items()}


class SeqnameMap:
    def __init__(self, names):
        self.names = set(names)
        self.renamed = collections.Counter()
        self.unknown = collections.Counter()

    def __call__(self, name):
        if name in self.names:
            return name
        for cand in ("chr" + name, "chrM" if name in ("M", "MT") else None):
            if cand and cand in self.names:
                self.renamed[f"{name}->{cand}"] += 1
                return cand
        self.unknown[name] += 1
        return None


def map_class(caller, native):
    key = native.split(":")[0].strip().lower()
    cls = CLASS_MAP[caller].get(key)
    if cls is None:
        die(f"{caller} ORF type '{native}' has no class mapping; add it to CLASS_MAP.")
    return cls


def bh_p_cutoff(pvals, fdr):
    """Largest raw p kept by Benjamini-Hochberg at `fdr` (-1.0 if none pass)."""
    ps = sorted(pvals)
    m, cut = len(ps), -1.0
    for i, p in enumerate(ps, 1):
        if p <= i / m * fdr:
            cut = p
    return cut


def price_fdr_cutoff(path, fdr):
    """PRICE's own default is -fdr 0.1, but allsamples.orfs.tsv lists every
    candidate it tested (39% of muscle rows have p >= 0.5). This applies that FDR
    to the table: BH over all its rows, returned as the raw-p cutoff."""
    with open(path) as fh:
        i = fh.readline().rstrip("\n").split("\t").index("p value")
        return bh_p_cutoff([float(l.rstrip("\n").split("\t")[i]) for l in fh if l.strip()], fdr)


def read_price(path, seqmap, max_p, max_fdr=None):
    rows, dropped = [], 0
    if max_fdr is not None:
        cut = price_fdr_cutoff(path, max_fdr)
        print(f"  PRICE BH FDR {max_fdr}: keep p <= {cut:.3g}", file=sys.stderr)
        max_p = cut if max_p is None else min(max_p, cut)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for need in ("Id", "Location", "Type", "p value", "Gene"):
            if need not in col:
                die(f"{path}: missing PRICE column '{need}'")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            p = float(f[col["p value"]])
            if max_p is not None and p > max_p:
                dropped += 1
                continue
            m = re.match(r"^(.+)([+-]):(.+)$", f[col["Location"]])
            if not m:
                die(f"{path}: unparsable Location '{f[col['Location']]}'")
            chrom = seqmap(m.group(1))
            if chrom is None:
                continue
            blocks = tuple(sorted(tuple(map(int, b.split("-"))) for b in m.group(3).split("|")))
            pid = f[col["Id"]]
            tid = re.sub(r"_[^_]+_\d+$", "", pid)
            gene = f[col["Gene"]]
            rows.append(dict(caller="price", sample="pooled", native_id=pid, chrom=chrom,
                             strand=m.group(2), blocks=blocks, transcript_id=tid, gene_id=gene,
                             native_type=f[col["Type"]], orf_class=map_class("price", f[col["Type"]]),
                             pval=p, bed_score=1000 if p <= 0 else min(1000, int(-math.log10(p) * 100))))
    return rows, dropped


def price_gene_from_transcript(rows, tx):
    """PRICE's Gene column is unreliable: the same ORF on the same transcript got
    olfactory-receptor or Hbb genes 140 kb away in some tissues, and 37-64% of
    PRICE-only ORFs lacked their own transcript's gene. Its transcript (from Id) is
    right, so take gene_id from the GTF. Returns how many rows changed."""
    changed = 0
    for r in rows:
        t = tx.get(r["transcript_id"])
        if t is None:
            die(f"PRICE transcript {r['transcript_id']} ({r['native_id']}) is not in the GTF")
        if r["gene_id"] != t[3]:
            r["gene_id"] = t[3]
            changed += 1
    return changed


def project(exons, strand, t0, t1):
    """Transcript interval [t0, t1) (0-based, 5'->3') -> genomic blocks."""
    ordered = exons if strand == "+" else exons[::-1]
    out, pos = [], 0
    for s, e in ordered:
        n = e - s
        lo, hi = max(t0, pos), min(t1, pos + n)
        if lo < hi:
            if strand == "+":
                out.append((s + lo - pos, s + hi - pos))
            else:
                out.append((e - (hi - pos), e - (lo - pos)))
        pos += n
    return tuple(sorted(out)) if sum(b - a for a, b in out) == t1 - t0 else None


def read_ribocode(path, sample, tx, seqmap, max_p, stats):
    rows, dropped = [], 0
    if os.path.getsize(path) == 0:
        return rows, dropped
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            p = float(f[col["pval_combined"]])
            if max_p is not None and p > max_p:
                dropped += 1
                continue
            tid, strand = f[col["transcript_id"]], f[col["strand"]]
            t0, t1 = int(f[col["ORF_tstart"]]) - 1, int(f[col["ORF_tstop"]])
            blocks = None
            if tid in tx:
                blocks = project(tx[tid][2], tx[tid][1], t0, t1)
            if blocks is None:
                stats["ribocode_single_block_fallback"] += 1
                gs, ge = sorted((int(f[col["ORF_gstart"]]), int(f[col["ORF_gstop"]])))
                blocks = ((gs - 1, ge),)
            chrom = seqmap(f[col["chrom"]])
            if chrom is None:
                continue
            rows.append(dict(caller="ribocode", sample=sample, native_id=f[col["ORF_ID"]],
                             chrom=chrom, strand=strand, blocks=blocks, transcript_id=tid,
                             gene_id=f[col["gene_id"]], native_type=f[col["ORF_type"]],
                             orf_class=map_class("ribocode", f[col["ORF_type"]]),
                             pval=p, bed_score=int((1 - p) * 1000)))
    return rows, dropped


# --------------------------------------------------------------------------- genome checks

def last_codon(genome, chrom, strand, blocks):
    if strand == "+":
        s, e = blocks[-1][1] - 3, blocks[-1][1]
        if s < blocks[-1][0]:
            return None
        return str(genome[chrom][s:e]).upper()
    s, e = blocks[0][0], blocks[0][0] + 3
    if e > blocks[0][1]:
        return None
    return str(genome[chrom][s:e]).upper().translate(COMP)[::-1]


def stop_position(strand, blocks):
    return blocks[-1][1] if strand == "+" else blocks[0][0]


# --------------------------------------------------------------------------- merge

def length(blocks):
    return sum(e - s for s, e in blocks)


def overlap(a, b):
    i = j = tot = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            tot += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def rank(p):
    return (CLASS_RANK[p["orf_class"]], -p["aa_length"], p["key"])


def same_junctions(a, b):
    """True if two ORFs have identical exon blocks wherever both span. Different starts
    or stops in a shared exon pass; an extra, skipped or shifted exon, or one ORF
    running into the other's intron, does not. ORFs whose spans don't overlap share
    nothing, so they are not the same."""
    lo, hi = max(a[0][0], b[0][0]), min(a[-1][1], b[-1][1])
    if lo >= hi:
        return False

    def clip(blocks):
        return [(max(s, lo), min(e, hi)) for s, e in blocks if s < hi and e > lo]
    return clip(a) == clip(b)


def cluster_by_overlap(proxies, min_ro):
    """Greedy complete linkage at >= min_ro reciprocal overlap on exon blocks, among ORFs
    with the same splice junctions where they overlap. Merging on overlap alone
    would fold splice isoforms (an extra 14-aa exon in a 500-aa protein is 97%
    overlap; a spliced Dnajc19-ps N-isoform was 85%) into one representative."""
    def close(p, q):
        if q["blocks"][0][0] >= p["blocks"][-1][1] or p["blocks"][0][0] >= q["blocks"][-1][1]:
            return False
        if not same_junctions(p["blocks"], q["blocks"]):
            return False
        ov = overlap(p["blocks"], q["blocks"])
        return ov >= min_ro * length(p["blocks"]) and ov >= min_ro * length(q["blocks"])

    clusters = []
    for p in sorted(proxies, key=rank):
        placed = False
        for c in clusters:
            if all(close(p, q) for q in c):
                c.append(p)
                placed = True
                break
        if not placed:
            clusters.append([p])
    return clusters


def merge(rows, min_ro, mode="exact"):
    proxies = {}
    for r in rows:
        key = (r["chrom"], r["strand"], r["blocks"])
        p = proxies.get(key)
        if p is None:
            p = proxies[key] = dict(key=key, chrom=r["chrom"], strand=r["strand"],
                                    blocks=r["blocks"], members=[],
                                    aa_length=(length(r["blocks"]) - 3) // 3)
        p["members"].append(r)
    for p in proxies.values():
        best = min(p["members"], key=lambda r: (CLASS_RANK[r["orf_class"]], r["caller"], r["native_id"]))
        p.update(orf_class=best["orf_class"], transcript_id=best["transcript_id"],
                 gene_id=best["gene_id"], bed_score=max(m["bed_score"] for m in p["members"]))
    if mode == "exact":
        return [[p] for p in proxies.values()]

    groups = collections.defaultdict(list)
    for p in proxies.values():
        if p["orf_class"] == "novel_u":
            groups[("ov", p["chrom"], p["strand"])].append(p)
        elif p["orf_class"] == "canonical_cds":
            groups[("ov", p["transcript_id"], p["strand"])].append(p)
        else:
            groups[("exact", p["transcript_id"], p["strand"],
                    p["blocks"][0][0], p["blocks"][-1][1])].append(p)

    clusters = []
    for gkey, ps in groups.items():
        clusters.extend(cluster_by_overlap(ps, min_ro) if gkey[0] == "ov" else [ps])
    return clusters


# --------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--price", help="PRICE allsamples.orfs.tsv")
    ap.add_argument("--ribocode", nargs="*", default=[], help="RiboCode <sample>.txt files")
    ap.add_argument("--gtf", required=True, help="hybrid GTF the callers were given")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--smorf-max-aa", type=int, default=150)
    ap.add_argument("--min-callers", type=int, default=1)
    ap.add_argument("--min-samples", type=int, default=1)
    ap.add_argument("--price-max-p", type=float, default=None)
    ap.add_argument("--price-max-fdr", type=float, default=None,
                    help="Benjamini-Hochberg FDR over every PRICE candidate (PRICE's default is 0.1)")
    ap.add_argument("--ribocode-max-p", type=float, default=None)
    ap.add_argument("--merge", choices=["exact", "overlap"], default="exact",
                    help="exact: one row per distinct ORF structure; overlap: fold by reciprocal overlap")
    ap.add_argument("--min-reciprocal-overlap", type=float, default=0.8,
                    help="--merge overlap only")
    ap.add_argument("--out-prefix", required=True, help="writes <prefix>.bed12/.tsv/.consensus.bed12")
    ap.add_argument("--report", required=True)
    a = ap.parse_args()

    genome = Fasta(a.genome, as_raw=True, sequence_always_upper=True)
    seqmap = SeqnameMap(genome.keys())
    tx = load_exons(a.gtf)
    stats = collections.Counter()

    dups = duplicate_inputs(a.ribocode)
    if dups:
        die("identical RiboCode outputs for different samples (a caller reused one "
            f"sample's data for another): {dups}")

    rows = []
    if a.price:
        r, d = read_price(a.price, seqmap, a.price_max_p, a.price_max_fdr)
        stats["price_gene_id_changed"] = price_gene_from_transcript(r, tx)
        rows += r
        stats["price_rows"], stats["price_dropped_p"] = len(r), d
    for path in a.ribocode:
        sample = os.path.basename(path).removesuffix(".txt")
        r, d = read_ribocode(path, sample, tx, seqmap, a.ribocode_max_p, stats)
        rows += r
        stats[f"ribocode_rows:{sample}"], stats[f"ribocode_dropped_p:{sample}"] = len(r), d
    if not rows:
        die("no ORF calls from any caller")
    if seqmap.unknown:
        die(f"seqnames not in the genome: {dict(seqmap.unknown)}")

    for c in CALLERS:
        mine = [r for r in rows if r["caller"] == c]
        if not mine:
            continue
        ends = collections.Counter(
            "stop" if last_codon(genome, r["chrom"], r["strand"], r["blocks"]) in STOPS else "no_stop"
            for r in mine)
        stats[f"{c}_ends_in_stop"] = ends["stop"]
        stats[f"{c}_ends_not_stop"] = ends["no_stop"]

    clusters = merge(rows, a.min_reciprocal_overlap, a.merge)
    reps = []
    for c in clusters:
        rep = min(c, key=rank)
        members = [m for p in c for m in p["members"]]
        called = {k: int(any(m["caller"] == k for m in members)) for k in CALLERS}
        best_p = {k: min((m["pval"] for m in members if m["caller"] == k), default=None) for k in CALLERS}
        samples = sorted({m["sample"] for m in members if m["sample"] != "pooled"})
        reps.append(dict(rep=rep, called=called, best_p=best_p, samples=samples,
                         native=",".join(sorted({f"{m['caller']}:{m['native_type']}" for m in members}))))

    reps.sort(key=lambda d: (d["rep"]["chrom"], d["rep"]["blocks"][0][0], d["rep"]["blocks"][-1][1],
                             d["rep"]["strand"], d["rep"]["gene_id"], d["rep"]["transcript_id"],
                             d["rep"]["orf_class"]))
    cols = (["orf_id", "chrom", "start", "end", "strand", "gene_id", "transcript_id", "orf_class",
             "aa_length", "is_smorf"] + [f"called_by_{c}" for c in CALLERS] + ["n_callers"]
            + [f"pval_{c}" for c in CALLERS] + ["n_samples", "samples", "orf_type_native",
                                                "stop_key"])
    by_stop = collections.defaultdict(set)
    n_cons = 0
    with open(a.out_prefix + ".bed12", "w") as bed, open(a.out_prefix + ".tsv", "w") as tsv, \
            open(a.out_prefix + ".consensus.bed12", "w") as cons:
        tsv.write("\t".join(cols) + "\n")
        for i, d in enumerate(reps, 1):
            p, oid = d["rep"], f"orf_{i:08d}"
            bl = p["blocks"]
            start, end = bl[0][0], bl[-1][1]
            line = "\t".join(map(str, [
                p["chrom"], start, end, oid, p["bed_score"], p["strand"], start, end, 0, len(bl),
                ",".join(str(e - s) for s, e in bl) + ",",
                ",".join(str(s - start) for s, e in bl) + ","])) + "\n"
            bed.write(line)
            n_callers = sum(d["called"].values())
            stop_key = f"{p['chrom']}:{p['strand']}:{stop_position(p['strand'], bl)}"
            for c in CALLERS:
                if d["called"][c]:
                    by_stop[stop_key].add(c)
            tsv.write("\t".join(map(str, [
                oid, p["chrom"], start, end, p["strand"], p["gene_id"], p["transcript_id"],
                p["orf_class"], p["aa_length"], int(0 < p["aa_length"] <= a.smorf_max_aa)]
                + [d["called"][c] for c in CALLERS] + [n_callers]
                + ["" if d["best_p"][c] is None else f"{d['best_p'][c]:.6g}" for c in CALLERS]
                + [len(d["samples"]), ",".join(d["samples"]), d["native"], stop_key])) + "\n")
            # PRICE is pooled, so a PRICE call carries no per-sample count.
            if n_callers >= a.min_callers and (len(d["samples"]) >= a.min_samples or d["called"]["price"]):
                cons.write(line)
                n_cons += 1

    stats["catalogue_orfs"] = len(reps)
    stats["consensus_orfs"] = n_cons
    stats["orfs_called_by_both"] = sum(all(d["called"].values()) for d in reps)
    stats["stop_codons_shared_by_both"] = sum(len(v) == len(CALLERS) for v in by_stop.values())
    for k, v in seqmap.renamed.items():
        stats[f"seqname_renamed:{k}"] = v
    with open(a.report, "w") as rep:
        rep.write("metric\tvalue\n")
        for k in sorted(stats):
            rep.write(f"{k}\t{stats[k]}\n")
    print(f"build_catalogue: {len(reps)} ORFs ({n_cons} consensus); "
          f"{stats['orfs_called_by_both']} called by both, "
          f"{stats['stop_codons_shared_by_both']} stop codons shared; "
          f"{sum(seqmap.renamed.values())} seqnames harmonised", file=sys.stderr)


if __name__ == "__main__":
    main()
