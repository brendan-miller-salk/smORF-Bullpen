#!/usr/bin/env python3
"""Novel transcripts of a StringTie merge that the RNA-seq supports.

Novel = a merged transcript whose transcript_id the reference GTF lacks (StringTie
--merge -G keeps reference ids; new models are MSTRG.*). One is kept if its -e TPM
is >= --min-tpm in >= --min-samples of the per-sample -e GTFs. Kept transcripts are
written with their merged GTF lines unchanged, in merged order; gffcompare then
assigns class codes and novel_class_codes picks which reach the hybrid GTFs.

  --out-gtf     the kept novel transcripts (the run's novel_gtf)
  --out-tpm     every novel transcript: locus, exons, TPM per sample, kept
  --out-report  counts at each step; unstranded (strand '.') models are kept here
                but the hybrid GTF build drops them
"""
import argparse
import collections
import os
import re
import sys

TID = re.compile(r'transcript_id "([^"]*)"')
TPM = re.compile(r'TPM "([^"]*)"')


def gtf_rows(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            t = TID.search(f[8])
            if t:
                yield t.group(1), f, line


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--expr", nargs="+", required=True, help="<sample>.gtf from stringtie -e")
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--min-samples", type=int, default=2)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-tpm", required=True)
    ap.add_argument("--out-report", required=True)
    a = ap.parse_args()

    ref = {tid for tid, _, _ in gtf_rows(a.reference)}
    lines, info, exons = collections.OrderedDict(), {}, collections.Counter()
    for tid, f, line in gtf_rows(a.merged):
        lines.setdefault(tid, []).append(line)
        if f[2] == "transcript":
            info[tid] = (f[0], f[3], f[4], f[6])
        elif f[2] == "exon":
            exons[tid] += 1
    novel = [t for t in lines if t not in ref]

    samples, tpm = [], {}
    for path in a.expr:
        s = os.path.basename(path).removesuffix(".gtf")
        samples.append(s)
        for tid, f, _ in gtf_rows(path):
            m = TPM.search(f[8]) if f[2] == "transcript" else None
            if m:
                tpm[(tid, s)] = float(m.group(1))

    kept = []
    with open(a.out_tpm, "w") as fh:
        fh.write("\t".join(["transcript_id", "chrom", "start", "end", "strand", "exons"]
                           + [f"TPM_{s}" for s in samples] + ["n_samples", "kept"]) + "\n")
        for tid in novel:
            vals = [tpm.get((tid, s), 0.0) for s in samples]
            n = sum(v >= a.min_tpm for v in vals)
            ok = n >= a.min_samples
            if ok:
                kept.append(tid)
            fh.write("\t".join([tid, *info.get(tid, (".",) * 4), str(exons[tid])]
                               + [f"{v:.4g}" for v in vals] + [str(n), str(ok)]) + "\n")
    with open(a.out_gtf, "w") as fh:
        for tid in kept:
            fh.writelines(lines[tid])

    report = [
        ("merged_transcripts", len(lines)),
        ("reference_transcripts", len(lines) - len(novel)),
        ("novel_transcripts", len(novel)),
        (f"novel_kept (TPM >= {a.min_tpm:g} in >= {a.min_samples} of {len(samples)})", len(kept)),
        ("novel_kept_multi_exon", sum(exons[t] > 1 for t in kept)),
        ("novel_kept_unstranded (dropped by the hybrid GTF)", sum(info.get(t, (".",) * 4)[3] == "." for t in kept)),
    ]
    with open(a.out_report, "w") as fh:
        fh.write("step\ttranscripts\n")
        fh.writelines(f"{k}\t{v}\n" for k, v in report)
    print(open(a.out_report).read(), file=sys.stderr, end="")
    if not kept:
        print("WARNING: no novel transcript passed; the hybrid GTFs will equal the reference.", file=sys.stderr)


if __name__ == "__main__":
    main()
