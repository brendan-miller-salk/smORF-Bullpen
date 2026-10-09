#!/usr/bin/env python3
"""DIAMOND blastx hits (species transcripts vs ORF proteins) -> per-ORF conservation.

Input per species: outfmt 6 (gzip or plain) with columns
    qseqid sseqid pident length qstart qend qlen qframe sstart send slen evalue bitscore
where the query is an Ensembl transcript and the subject is an ORF, plus the
transcript FASTA the queries came from and, optionally, a run-species -> species ortholog
table (refs/fetch_conservation.sh). Reverse frames are dropped: Ensembl transcripts
are sense.

E-values are recomputed as tblastn would report them for the ORF as query:
    E = slen * database aa * 2^-bitscore     (Karlin-Altschul, no edge correction)
DIAMOND's own E-value uses the ORF set as the database, so it changes with the run;
the bitscore does not. Two databases per species:
    transcriptome  every transcript (nt / 3)                          -> _evalue
    ortholog       transcripts of the orthologs of the ORF's host gene -> _ortholog_evalue
A hit in an ortholog is judged on its ortholog E. A short ORF can only reach a
significant E in the small ortholog space: a perfect 15-aa match is ~33 bits,
E ~0.3 against the human transcriptome, ~5e-6 against one gene's transcripts.

Per ORF and species the best hit is the one that passes (E <= --evalue and ORF
coverage (send - sstart + 1) / slen >= --min-cov), else the highest bitscore; ties
go to the first in the file. _call is
    conserved      the best hit passes
    not_conserved  it fails, but a perfect full-length match would have passed in
                   at least one of the databases
    untestable     even a perfect match could not pass (too short for the space),
                   or the ORF has no clean protein (internal stop)
A perfect match's bitscore is from BLOSUM62 self-scores with DIAMOND's gapped
BLOSUM62 11/1 statistics (lambda 0.267, K 0.041), plus PERFECT_SLACK bits:
MEASURED on 32 perfect full-length hits (muscle vs human), DIAMOND scored up to
1.15 bits above this (median 1.5 below, composition correction).

Output (first column orf_id, one row per ORF in --bed12):
    diamond_n_species      species where conserved
    diamond_species        those species, comma-separated, in --hits order
    diamond_n_untestable   species where untestable
    diamond_<sp>_call, _ortholog (the host gene's ortholog symbols), and the best
    hit's _gene (gene_symbol, else gene id), _biotype, _pident, _cov, _evalue,
    _ortholog_evalue (blank if the hit is outside the orthologs)
"""
import argparse
import csv
import gzip
import math
import re
import sys

COLS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen", "qframe",
        "sstart", "send", "slen", "evalue", "bitscore"]
FIELD = re.compile(r"(\w+):(\S+)")
BLOSUM62_SELF = dict(A=4, R=5, N=6, D=6, C=9, Q=5, E=5, G=6, H=8, I=4, L=4, K=5, M=5,
                     F=6, P=7, S=4, T=5, W=11, Y=7, V=4)
LAMBDA, K = 0.267, 0.041
PERFECT_SLACK = 1.5
PER_SPECIES = ("call", "ortholog", "gene", "biotype", "pident", "cov", "evalue",
               "ortholog_evalue")


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _unversion(gene_id):
    return gene_id.split(".")[0]


def header_meta(title):
    """Ensembl FASTA header (without '>') -> (transcript id, gene id, symbol, biotype)."""
    tid, _, rest = title.partition(" ")
    kv = dict(FIELD.findall(rest))
    gid = _unversion(kv.get("gene", "")) or tid
    return tid, gid, kv.get("gene_symbol") or gid, kv.get("transcript_biotype", "")


def read_transcripts(path):
    """-> ({tid: (gene id, symbol, biotype)}, {gene id: nt}, database aa = nt / 3)."""
    meta, gene_nt, gid = {}, {}, None
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                tid, gid, sym, biotype = header_meta(line[1:].rstrip("\n"))
                meta[tid] = (gid, sym, biotype)
                gene_nt.setdefault(gid, 0)
            elif gid is not None:
                gene_nt[gid] += len(line.rstrip("\n"))
    return meta, gene_nt, sum(gene_nt.values()) / 3


def read_orthologs(path):
    """run-species gene -> {ortholog gene}; columns <source>_gene_id, ortholog_gene_id[, type]."""
    orth = {}
    with _open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and f[1]:
                orth.setdefault(_unversion(f[0]), set()).add(_unversion(f[1]))
    return orth


def read_proteins(path):
    seqs, oid = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                oid = line[1:].split()[0]
                seqs[oid] = []
            elif oid is not None:
                seqs[oid].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def perfect_bits(aa):
    raw = sum(BLOSUM62_SELF.get(c, 0) for c in aa.upper())
    return (LAMBDA * raw - math.log(K)) / math.log(2) + PERFECT_SLACK


def evalue(slen, bits, letters):
    return slen * letters * 2.0 ** -bits


def host_genes(catalogue, annotator):
    """orf_id -> {run-species gene}: the Annotator's host_gene_id, else the catalogue's
    gene_id when it is an Ensembl gene (StringTie MSTRG ids have no orthologs)."""
    genes = {}
    with open(catalogue) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r.get("gene_id") or ""
            if g.startswith("ENS"):
                genes[r["orf_id"]] = {_unversion(x) for x in g.split(",")}
    if annotator:
        with open(annotator) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                g = r.get("host_gene_id") or ""
                if g.startswith("ENS"):
                    genes[r["orf_id"]] = {_unversion(x) for x in g.split(",")}
    return genes


def best_hits(path, meta, letters, orth_of, evalue_max, min_cov):
    """orf_id -> (passes, bitscore, row) for the best forward-frame hit.
    orth_of: orf_id -> (ortholog gene set, ortholog database aa)."""
    best = {}
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(COLS):
                sys.exit(f"diamond_classify: {path}: expected {len(COLS)} columns, got {len(f)}")
            h = dict(zip(COLS, f))
            if int(h["qframe"]) < 0:
                continue
            oid, slen, bits = h["sseqid"], int(h["slen"]), float(h["bitscore"])
            h["cov"] = (abs(int(h["send"]) - int(h["sstart"])) + 1) / slen
            h["E"] = evalue(slen, bits, letters)
            genes, orth_aa = orth_of.get(oid, (set(), 0))
            in_orth = meta.get(h["qseqid"], ("",))[0] in genes
            h["E_orth"] = evalue(slen, bits, orth_aa) if in_orth else None
            e = h["E_orth"] if in_orth else h["E"]
            ok = e <= evalue_max and h["cov"] >= min_cov
            key = (ok, bits)
            if oid not in best or key > best[oid][:2]:
                best[oid] = (*key, h)
    return best


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hits", nargs="+", required=True, help="species=blastx.tsv[.gz], in output order")
    p.add_argument("--transcripts", nargs="+", required=True,
                   help="species=transcripts.fa[.gz], the blastx queries")
    p.add_argument("--orthologs", nargs="*", default=[],
                   help="species=orthologs.tsv (run-species gene -> species gene); species without "
                        "one get no ortholog database")
    p.add_argument("--proteins", required=True, help="the ORF proteins DIAMOND searched (orfs.faa)")
    p.add_argument("--catalogue", required=True, help="cohort.catalogue.tsv (orf_id, gene_id)")
    p.add_argument("--annotator", help="annotations/annotator.tsv (orf_id, host_gene_id)")
    p.add_argument("--bed12", required=True, help="ORF universe: every ORF gets a row")
    p.add_argument("--evalue", type=float, default=1e-3)
    p.add_argument("--min-cov", type=float, default=0.5)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    species = [s.split("=", 1) for s in a.hits]
    fasta = dict(s.split("=", 1) for s in a.transcripts)
    orth_paths = dict(s.split("=", 1) for s in a.orthologs)
    missing = [sp for sp, _ in species if sp not in fasta]
    if missing:
        sys.exit(f"diamond_classify: no --transcripts for {missing}")
    orfs = [line.split("\t")[3] for line in open(a.bed12) if line.strip()]
    prot = read_proteins(a.proteins)
    maxbits = {o: perfect_bits(s) for o, s in prot.items()}
    hosts = host_genes(a.catalogue, a.annotator)

    cols = {}      # species -> {orf_id: [per-species cells]}
    summary = {}
    for sp, path in species:
        meta, gene_nt, letters = read_transcripts(fasta[sp])
        orth = read_orthologs(orth_paths[sp]) if sp in orth_paths else {}
        sym = {gid: s for gid, s, _ in meta.values()}
        orth_of = {}
        for o in orfs:
            genes = set().union(*(orth.get(g, set()) for g in hosts.get(o, ())))
            genes &= gene_nt.keys()
            if genes:
                orth_of[o] = (genes, sum(gene_nt[g] for g in genes) / 3)
        best = best_hits(path, meta, letters, orth_of, a.evalue, a.min_cov)
        n = {"conserved": 0, "not_conserved": 0, "untestable": 0}
        cols[sp] = {}
        for o in orfs:
            genes, orth_aa = orth_of.get(o, (set(), 0))
            if o in prot:
                slen, mb = len(prot[o]), maxbits[o]
                testable = (evalue(slen, mb, letters) <= a.evalue
                            or (orth_aa > 0 and evalue(slen, mb, orth_aa) <= a.evalue))
            else:
                testable = False
            hit = best.get(o)
            call = ("conserved" if hit and hit[0] else
                    "not_conserved" if testable else "untestable")
            n[call] += 1
            cells = [call, ",".join(sorted(sym.get(g, g) for g in genes))]
            if hit:
                h = hit[2]
                _, gsym, biotype = meta.get(h["qseqid"], ("", h["qseqid"], ""))
                cells += [gsym, biotype, h["pident"], round(h["cov"], 3), f"{h['E']:.2e}",
                          "" if h["E_orth"] is None else f"{h['E_orth']:.2e}"]
            else:
                cells += [""] * 6
            cols[sp][o] = cells
        summary[sp] = (n, len(orth_of), letters)

    hdr = ["orf_id", "diamond_n_species", "diamond_species", "diamond_n_untestable"]
    for sp, _ in species:
        hdr += [f"diamond_{sp}_{c}" for c in PER_SPECIES]
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(hdr)
        for o in orfs:
            calls = [(sp, cols[sp][o][0]) for sp, _ in species]
            cons = [sp for sp, c in calls if c == "conserved"]
            unt = sum(c == "untestable" for _, c in calls)
            w.writerow([o, len(cons), ",".join(cons), unt,
                        *(x for sp, _ in species for x in cols[sp][o])])

    print(f"diamond: {len(orfs):,} ORFs ({len(orfs) - sum(o in prot for o in orfs):,} without a "
          f"clean protein); E <= {a.evalue:g}, cov >= {a.min_cov:g}", file=sys.stderr)
    for sp, (n, n_orth, letters) in summary.items():
        print(f"  {sp:<12} conserved {n['conserved']:>7,}  not_conserved {n['not_conserved']:>7,}  "
              f"untestable {n['untestable']:>7,}  with orthologs {n_orth:>7,}  "
              f"transcriptome {letters:.3g} aa", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
