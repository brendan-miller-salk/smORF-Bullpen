"""rnaseq_strand.py (strandedness from STAR gene counts) and assembly_filter.py (novel
StringTie transcripts with TPM support). Prints PASS/FAIL, exits 1 on failure."""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
STRAND = os.path.join(HERE, "..", "workflow", "scripts", "rnaseq_strand.py")
FILTER = os.path.join(HERE, "..", "workflow", "scripts", "assembly_filter.py")
failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


def reads_per_gene(d, sample, fr, rf):
    os.makedirs(os.path.join(d, sample), exist_ok=True)
    p = os.path.join(d, sample, "ReadsPerGene.out.tab")
    open(p, "w").write("N_unmapped\t9\t9\t9\nN_noFeature\t9\t9\t9\n"
                       f"G1\t{fr + rf}\t{fr}\t{rf}\nG2\t0\t0\t0\n")
    return p


def strand(d, paths, mode="auto"):
    out, flag = os.path.join(d, "rep.tsv"), os.path.join(d, "flag.txt")
    if os.path.exists(flag):
        os.remove(flag)
    r = subprocess.run([sys.executable, STRAND, "--counts", *paths, "--mode", mode,
                        "--out", out, "--out-flag", flag], capture_output=True, text=True)
    return r.returncode, open(flag).read().strip() if os.path.exists(flag) else None


with tempfile.TemporaryDirectory() as d:
    a, b = reads_per_gene(d, "A", 50, 950), reads_per_gene(d, "B", 30, 970)
    check("dUTP libraries -> rf", strand(d, [a, b]), (0, "rf"))
    u = reads_per_gene(d, "U", 500, 520)
    check("unstranded", strand(d, [u]), (0, "unstranded"))
    check("mixed cohort fails under auto", strand(d, [a, u])[0], 1)
    m = reads_per_gene(d, "M", 300, 700)
    check("ambiguous (0.7) fails", strand(d, [m])[0], 1)
    check("explicit mode wins over the reads", strand(d, [a, u], "unstranded"), (0, "unstranded"))
    check("N_* rows ignored (fr)", strand(d, [reads_per_gene(d, "F", 990, 10)]), (0, "fr"))

    def strict(paths, mode="auto"):
        out, flag = os.path.join(d, "rep.tsv"), os.path.join(d, "flag.txt")
        return subprocess.run([sys.executable, STRAND, "--counts", *paths, "--mode", mode,
                               "--require-stranded", "--out", out, "--out-flag", flag],
                              capture_output=True, text=True).returncode
    check("x guard: stranded cohort passes", strict([a, b]), 0)
    check("x guard: unstranded cohort fails", strict([u]), 1)
    check("x guard: explicit rf does not override unstranded reads", strict([u], "rf"), 1)


def gtf(tid, feat, s, e, strand="+", extra=""):
    return f'chr1\tStringTie\t{feat}\t{s}\t{e}\t1000\t{strand}\t.\tgene_id "G"; transcript_id "{tid}";{extra}\n'


def tx(tid, exons, strand="+", tpm=None):
    attr = f' TPM "{tpm}";' if tpm is not None else ""
    rows = [gtf(tid, "transcript", exons[0][0], exons[-1][1], strand, attr)]
    return "".join(rows + [gtf(tid, "exon", s, e, strand) for s, e in exons])


with tempfile.TemporaryDirectory() as d:
    ref = os.path.join(d, "ref.gtf")
    open(ref, "w").write(tx("ENST1", [(100, 200), (300, 400)]))
    two = [(1000, 1100), (1200, 1300)]
    merged = os.path.join(d, "merged.gtf")
    open(merged, "w").write("# stringtie --merge\n" + tx("ENST1", [(100, 200), (300, 400)])
                            + tx("MSTRG.1.1", two) + tx("MSTRG.2.1", [(5000, 5600)], ".")
                            + tx("MSTRG.3.1", two))
    tpms = {"S1": {"ENST1": 50, "MSTRG.1.1": 3, "MSTRG.2.1": 2, "MSTRG.3.1": 5},
            "S2": {"ENST1": 40, "MSTRG.1.1": 1, "MSTRG.2.1": 4, "MSTRG.3.1": 0.5},
            "S3": {"ENST1": 30, "MSTRG.1.1": 0.2, "MSTRG.2.1": 0, "MSTRG.3.1": 0}}
    expr = []
    for s, t in tpms.items():
        expr.append(os.path.join(d, f"{s}.gtf"))
        open(expr[-1], "w").write("".join(tx(k, [(1, 2)], tpm=v) for k, v in t.items()))
    out = {k: os.path.join(d, k) for k in ("novel.gtf", "tpm.tsv", "report.tsv")}
    r = subprocess.run([sys.executable, FILTER, "--merged", merged, "--reference", ref,
                        "--expr", *expr, "--min-tpm", "1", "--min-samples", "2",
                        "--out-gtf", out["novel.gtf"], "--out-tpm", out["tpm.tsv"],
                        "--out-report", out["report.tsv"]], capture_output=True, text=True)
    check("filter exit 0", r.returncode, 0)
    kept = [l.split('transcript_id "')[1].split('"')[0] for l in open(out["novel.gtf"])]
    check("kept: TPM >= 1 in >= 2 samples, reference left out", sorted(set(kept)), ["MSTRG.1.1", "MSTRG.2.1"])
    check("lines copied (transcript + exons)", len(kept), 3 + 2)
    rep = dict(l.rstrip("\n").split("\t") for l in list(open(out["report.tsv"]))[1:])
    check("report novel", rep["novel_transcripts"], "3")
    check("report multi-exon kept", rep["novel_kept_multi_exon"], "1")
    check("report unstranded kept", rep["novel_kept_unstranded (dropped by the hybrid GTF)"], "1")
    row = [l.split("\t") for l in open(out["tpm.tsv"]) if l.startswith("MSTRG.3.1")][0]
    check("tpm table: MSTRG.3.1 n_samples, kept", row[-2:], ["1", "False\n"])

REPORT = os.path.join(HERE, "..", "workflow", "scripts", "fastp_report.py")
with tempfile.TemporaryDirectory() as d:
    import json

    def fastp_json(name, r1, r2, before, cut):
        p = os.path.join(d, f"{name}.fastp.json")
        json.dump({"summary": {"before_filtering": {"total_reads": before},
                               "after_filtering": {"total_reads": before - 10}},
                   "adapter_cutting": {"read1_adapter_sequence": r1, "read2_adapter_sequence": r2,
                                       "adapter_trimmed_reads": cut}}, open(p, "w"))
        return p
    T1, T2 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    out = os.path.join(d, "rep.tsv")
    r = subprocess.run([sys.executable, REPORT, "--json", fastp_json("S_L0", T1, T2, 1000, 50),
                        fastp_json("S_L1", T1, T2, 2000, 0), "--out", out], capture_output=True, text=True)
    rows = [l.rstrip("\n").split("\t") for l in open(out)][1:]
    check("adapter report: one row per lane", [x[0] for x in rows], ["S_L0", "S_L1"])
    check("adapter report: detected adapter + pct", (rows[0][1], rows[0][6]), (T1, "5.00"))
    check("adapter report: same adapters, no warning", "WARNING" in r.stderr, False)
    r = subprocess.run([sys.executable, REPORT, "--json", fastp_json("S_L0", T1, T2, 1000, 50),
                        fastp_json("S_L2", "unspecified", "unspecified", 10, 0), "--out", out],
                       capture_output=True, text=True)
    check("adapter report: differing lanes warn, exit 0", ("WARNING" in r.stderr, r.returncode), (True, 0))

sys.exit(1 if failed else 0)
