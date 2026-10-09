"""UniProt classifier (incl. HSP chaining and insertions/deletions) and relabel checks. Prints PASS/FAIL, exits 1 on failure.

Synthetic BLAST rows, one per class, against a small UniProt FASTA: Swiss-Prot
canonical P1 (100 aa), its isoform P1-2 (same sequence, to test the order), TrEMBL T1
(P1's last 60 aa), T2 (= P1) and T4 (a 7-aa peptide). Exact identity is by sequence
(query FASTA), not from BLAST.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "workflow", "scripts", "blastoff")
sys.path.insert(0, SCRIPTS)
from blastoff import relabel_smorf_type  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


def hit(q, s, pident, qs, qe, ql, ss, se, sl, ev="1e-30", bit="200"):
    return "\t".join(map(str, [q, s, pident, qe - qs + 1, 0, 0, qs, qe, ql, ss, se, sl, ev, bit]))


def gapped(q, s, ident, length, mismatch, gapopen, qs, qe, ql, ss, se, sl, bit="200"):
    """A BLAST row with gaps: pident counts gap columns against identity."""
    return "\t".join(map(str, [q, s, round(100 * ident / length, 3), length, mismatch, gapopen,
                                qs, qe, ql, ss, se, sl, "1e-30", bit]))


AAS = "ACDEFGHIKLMNPQRSTVWY"
P1 = "M" + "".join(AAS[(i * 7 + 3) % 20] for i in range(99))
T1 = P1[40:]
T4 = "MKVLRWS"
FASTA = f""">sp|P1|PROT1_MOUSE Protein one OS=Mus musculus OX=10090 GN=Prot1 PE=1 SV=1
{P1[:60]}
{P1[60:]}
>sp|P1-2|PROT1_MOUSE Isoform 2 of Protein one OS=Mus musculus OX=10090 GN=Prot1
{P1}
>tr|T1|T1_MOUSE Protein one (Fragment) OS=Mus musculus OX=10090 GN=Prot1 PE=4 SV=1
{T1}
>tr|T2|T2_MOUSE Protein one OS=Mus musculus OX=10090 GN=Prot1 PE=4 SV=1
{P1}
>tr|T4|T4_MOUSE Peptide OS=Mus musculus OX=10090 GN=Pep4 PE=4 SV=1
{T4}
"""
QUERIES = {
    "orf_identical": P1,
    "orf_tr": T1,
    "orf_indel1": P1[:50] + "W" + P1[50:],   # one extra residue: near, not exact
    "orf_near": P1[:30] + "WWW" + P1[33:],   # 3 substitutions
    "orf_exact_nohit": T4,                   # BLAST found nothing; still exact
    "orf_short": P1[29:37],                  # 8 aa inside P1
}
ROWS = [
    hit("orf_identical", "sp|P1|PROT1_MOUSE", 100, 1, 100, 100, 1, 100, 100),
    hit("orf_identical", "sp|P1-2|PROT1_MOUSE", 100, 1, 100, 100, 1, 100, 100),  # loses: isoform
    hit("orf_ntrunc", "sp|P1|PROT1_MOUSE", 100, 1, 60, 60, 41, 100, 100),
    hit("orf_tr", "sp|P1|PROT1_MOUSE", 100, 1, 60, 60, 41, 100, 100),
    hit("orf_tr", "T1", 100, 1, 60, 60, 1, 60, 60),  # tr exact beats sp truncation
    hit("orf_near", "sp|P1|PROT1_MOUSE", 97, 1, 100, 100, 1, 100, 100),
    hit("orf_short", "sp|P1|PROT1_MOUSE", 100, 1, 8, 8, 30, 37, 100, "0.5", "25"),  # 8 aa: chance
    hit("orf_ctrunc", "sp|P1|PROT1_MOUSE", 99, 1, 50, 50, 1, 50, 100),
    hit("orf_internal", "sp|P1|PROT1_MOUSE", 100, 1, 30, 30, 20, 49, 100),
    hit("orf_ext", "sp|P1|PROT1_MOUSE", 100, 21, 120, 120, 1, 100, 100),
    hit("orf_cext", "sp|P1|PROT1_MOUSE", 100, 1, 100, 115, 1, 100, 100),
    hit("orf_ncext", "sp|P1|PROT1_MOUSE", 100, 6, 105, 112, 1, 100, 100),
    hit("orf_homolog", "sp|P1|PROT1_MOUSE", 60, 1, 80, 80, 1, 80, 100, "1e-10", "90"),
    hit("orf_weak", "sp|P1|PROT1_MOUSE", 50, 5, 20, 40, 5, 20, 100, "3", "20"),
    hit("orf_partial", "sp|P1|PROT1_MOUSE", 100, 1, 20, 40, 11, 30, 100, "1e-8", "45"),  # shared N-term, own C-term
    hit("orf_tiny", "sp|P1|PROT1_MOUSE", 100, 1, 6, 20, 11, 16, 100, "5", "18"),  # 6 aa: chance
    hit("orf_slack", "sp|P1|PROT1_MOUSE", 96, 2, 58, 59, 43, 99, 100),  # within 2 aa slack
    # extra 20-aa exon bridged by one gapped HSP: pident 83 but 100% over ungapped columns
    gapped("orf_ins", "sp|P1|PROT1_MOUSE", 100, 120, 0, 1, 1, 120, 120, 1, 100, 100),
    # skipped 20-aa exon
    gapped("orf_del", "sp|P1|PROT1_MOUSE", 80, 100, 0, 1, 1, 80, 80, 1, 100, 100),
    # 40-aa insertion BLAST split into two HSPs, chained back into one hit
    hit("orf_ins_split", "sp|P1|PROT1_MOUSE", 100, 1, 50, 140, 1, 50, 100, bit="100"),
    hit("orf_ins_split", "sp|P1|PROT1_MOUSE", 100, 91, 140, 140, 51, 100, 100, bit="95"),
    # HSPs in opposite order on the two sequences (a repeat): not chained; best one is partial
    hit("orf_swapped", "sp|P1|PROT1_MOUSE", 100, 1, 50, 100, 51, 100, 100, bit="100"),
    hit("orf_swapped", "sp|P1|PROT1_MOUSE", 100, 51, 100, 100, 1, 50, 100, bit="90"),
    # 1-aa indel is within slack: near-identical (not the same sequence)
    gapped("orf_indel1", "sp|P1|PROT1_MOUSE", 100, 101, 0, 1, 1, 101, 101, 1, 100, 100),
    # gaps do not rescue mismatches: 70% over ungapped columns stays a homolog
    gapped("orf_gaphomolog", "sp|P1|PROT1_MOUSE", 70, 110, 30, 2, 1, 110, 110, 1, 100, 100, bit="120"),
]
WANT = {
    "orf_identical": ("sp_identical", "0", "P1,P1-2,T2"),
    "orf_ntrunc": ("sp_N_truncated", "40", ""),
    "orf_tr": ("tr_identical", "0", "T1"),
    "orf_ctrunc": ("sp_C_truncated", "50", ""),
    "orf_internal": ("sp_internal", "19", ""),
    "orf_ext": ("sp_N_extended", "20", ""),
    "orf_cext": ("sp_C_extended", "15", ""),
    "orf_ncext": ("sp_NC_extended", "12", ""),
    "orf_homolog": ("sp_homolog", "", ""),
    "orf_weak": ("novel", "", ""),
    "orf_slack": ("sp_N_truncated", "42", ""),
    "orf_partial": ("sp_partial", "10", ""),
    "orf_tiny": ("novel", "", ""),
    "orf_nohit": ("novel", "", ""),
    "orf_ins": ("sp_insertion", "20", ""),
    "orf_del": ("sp_deletion", "20", ""),
    "orf_ins_split": ("sp_insertion", "40", ""),
    "orf_swapped": ("sp_partial", "50", ""),
    "orf_indel1": ("sp_near_identical", "0", ""),
    "orf_near": ("sp_near_identical", "0", ""),
    "orf_exact_nohit": ("tr_identical", "0", "T4"),
    "orf_short": ("novel", "", ""),
    "orf_gaphomolog": ("sp_homolog", "", ""),
}

with tempfile.TemporaryDirectory() as d:
    fa, bl, bed, out, qf = (os.path.join(d, x) for x in ("u.fa", "b.tsv", "o.bed12", "o.tsv", "q.fa"))
    open(fa, "w").write(FASTA)
    open(qf, "w").write("".join(f">{k}\n{v}\n" for k, v in QUERIES.items()))
    open(bl, "w").write("\n".join(ROWS) + "\n")
    with open(bed, "w") as fh:
        for oid in WANT:
            fh.write(f"chr1\t0\t30\t{oid}\t0\t+\t0\t30\t0\t1\t30,\t0,\n")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "blastoff.py"), "classify",
                        "--blast", bl, "--uniprot-fasta", fa, "--bed12", bed, "--queries", qf,
                        "--out", out],
                       capture_output=True, text=True)
    check("classify exit", r.returncode, 0)
    if r.returncode:
        print(r.stderr)
    rows = [l.rstrip("\n").split("\t") for l in open(out)]
    hdr, body = rows[0], {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}
    check("orf_id first", hdr[0], "orf_id")
    check("one row per ORF", sorted(body), sorted(WANT))
    for oid, (cls, off, ident) in WANT.items():
        got = body.get(oid, {})
        check(f"{oid} class", got.get("blastoff_class"), cls)
        check(f"{oid} offset", got.get("blastoff_offset_aa"), off)
        check(f"{oid} identical_to", got.get("blastoff_identical_to"), ident)
    check("canonical beats isoform", body["orf_identical"]["blastoff_acc"], "P1")
    check("exact pident", body["orf_exact_nohit"]["blastoff_pident"], "100.0")
    check("exact qcov", body["orf_exact_nohit"]["blastoff_qcov"], "1.000")
    for oid, d in [("orf_identical", ""), ("orf_ntrunc", "N"), ("orf_ctrunc", "C"),
                   ("orf_internal", "NC"), ("orf_ext", "N"), ("orf_cext", "C"),
                   ("orf_ncext", "NC"), ("orf_partial", "C"), ("orf_homolog", ""),
                   ("orf_ins", ""), ("orf_del", ""), ("orf_swapped", "C")]:
        check(f"{oid} diverges", body[oid]["blastoff_diverges"], d)
    check("gene parsed", body["orf_identical"]["blastoff_gene"], "Prot1")
    check("pass_bit40 strong", body["orf_homolog"]["blastoff_pass_bit40"], "True")
    check("pass_bit40 weak", body["orf_weak"]["blastoff_pass_bit40"], "False")

for smorf, cls, want in [
    ("iORF", "sp_identical", "annotated_CDS"),
    ("iORF", "tr_identical", "TrEMBL"),
    ("iORF", "tr_N_truncated", "iORF"),       # TrEMBL label only for 100% matches
    ("iORF", "tr_homolog", "iORF"),
    ("iORF", "sp_N_truncated", "Trunc-Iso"),
    ("iORF", "tr_C_truncated", "iORF"),
    ("iORF", "sp_internal", "Trunc-Iso"),
    ("iORF", "sp_insertion", "Ins-Iso"),
    ("iORF", "sp_deletion", "Del-Iso"),
    ("iORF", "tr_insertion", "iORF"),
    ("iORF", "tr_partial", "iORF"),
    ("iORF", "sp_N_extended", "N-Iso"),
    ("iORF", "sp_C_extended", "C-Iso"),
    ("iORF", "tr_C_extended", "iORF"),
    ("iORF", "sp_homolog", "iORF"),
    ("iORF", "novel", "iORF"),
    ("iORF", float("nan"), "iORF"),
    ("uORF", "sp_N_extended", "uORF"),      # only iORF + overlapping/ambiguous types
    ("uORF", "sp_insertion", "uORF"),
    ("uORF", "sp_partial", "uORF"),
    ("psORF", "sp_identical", "annotated_CDS"),  # identical relabels any type
    ("intergenic", "sp_identical", "annotated_CDS"),
    ("lncRNA", "tr_identical", "TrEMBL"),
    ("lncRNA", "tr_partial", "lncRNA"),
    ("iORF", "sp_near_identical", "iORF"),       # not the protein: keeps its type
    ("lncRNA", "tr_near_identical", "lncRNA"),
]:
    check(f"relabel {smorf}+{cls}", relabel_smorf_type(smorf, cls), want)

for smorf, cls, div, want in [
    ("uoORF", "sp_N_extended", "N", "N-Iso"),
    ("uoORF", "sp_partial", "N", "N-Iso"),
    ("doORF", "sp_C_extended", "C", "C-Iso"),
    ("doORF", "sp_partial", "C", "C-Iso"),
    ("daoORF", "sp_N_truncated", "N", "Trunc-Iso"),   # inside the protein, adds nothing
    ("doORF", "sp_C_truncated", "C", "Trunc-Iso"),
    ("uaoORF", "sp_NC_extended", "NC", "NC-Iso"),
    ("uaORF", "sp_internal", "NC", "Trunc-Iso"),
    ("uoORF", "sp_insertion", "", "Ins-Iso"),
    ("daoORF", "sp_deletion", "", "Del-Iso"),
    ("daORF", "sp_identical", "", "annotated_CDS"),
    ("uoORF", "sp_near_identical", "", "uoORF"),
    ("uoORF", "tr_N_extended", "N", "uoORF"),
    ("uoORF", "sp_homolog", "", "uoORF"),       # divergent: not in frame
    ("doORF", "novel", "", "doORF"),
    ("iORF", "sp_partial", "C", "C-Iso"),   # own C-terminal residues: not a truncation
    ("iORF", "sp_partial", "N", "N-Iso"),
]:
    check(f"relabel {smorf}+{cls}+{div}", relabel_smorf_type(smorf, cls, div), want)

print(f"\n{'ALL PASS' if not failed else f'{failed} FAILED'}")
sys.exit(1 if failed else 0)
