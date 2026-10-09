#!/usr/bin/env python3
"""Convert scored ORFs (BED12) into the GTF the Annotator expects.

The Annotator is used unmodified, so this has to meet its input contract
exactly. Three things it requires, none of which a caller's BED12 provides:

1. **CDS features.** It skips any row where `parts[2] != 'CDS'`.

2. **Coordinates embedded in `gene_id`.** It recovers the true ORF span with
   `COORD_RE = r'(chr[\\dXYM]+):(\\d+)-(\\d+)'` applied to the gene_id string,
   because a single GTF row is one block, not the whole ORF. With no match it
   silently falls back to that row's own coordinates -- the failure mode the
   Annotator's own notes blame for ~50k inflated uoORF calls. So every ORF is
   renamed to `<orf_id>-chrN:start-end`, matching the ShortStop convention.

3. **`chr`-prefixed seqnames.** It hard-codes `line.startswith('chr')`.

That last point drops ORFs on scaffolds (GL*/JH*/MU*) -- `COORD_RE` cannot
match them either. Rather than let those vanish, they are written to a
separate file and counted, so they can be rejoined as `not_annotated_scaffold`
instead of quietly disappearing from the results.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys

COORD_RE = re.compile(r'(chr[\dXYM]+):(\d+)-(\d+)')  # the Annotator's own pattern
ANNOTATABLE = re.compile(r'^chr[\dXYM]+$')           # what it can actually handle


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed12", required=True)
    p.add_argument("--out-gtf", required=True)
    p.add_argument("--out-skipped", required=True,
                   help="TSV of ORFs the Annotator cannot handle (scaffolds).")
    p.add_argument("--out-idmap", required=True,
                   help="TSV: original orf_id -> renamed gene_id, to join results back.")
    p.add_argument("--source", default="riboseq")
    args = p.parse_args(argv)

    n_in = n_out = n_skip = n_preexisting = 0

    with _open(args.bed12) as fh, \
         open(args.out_gtf, "w") as out, \
         open(args.out_skipped, "w") as skipped, \
         open(args.out_idmap, "w") as idmap:

        skipped.write("orf_id\tchrom\treason\n")
        idmap.write("orf_id\tannotator_gene_id\n")

        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 12:
                raise ValueError(f"expected 12+ BED columns, got {len(f)}")
            n_in += 1

            chrom, cs, ce, name, strand = f[0], int(f[1]), int(f[2]), f[3], f[5]

            if not ANNOTATABLE.match(chrom):
                skipped.write(f"{name}\t{chrom}\tseqname not matched by the Annotator's chr[0-9XYM]+ pattern\n")
                n_skip += 1
                continue

            # Our coordinate must be the FIRST COORD_RE hit in the id, since the
            # Annotator uses .search(). Warn if the caller's own id already has one.
            if COORD_RE.search(name):
                n_preexisting += 1
            gene_id = f"{name}-{chrom}:{cs + 1}-{ce}"   # 1-based inclusive span
            idmap.write(f"{name}\t{gene_id}\n")

            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offs = [int(x) for x in f[11].rstrip(",").split(",")]
            blocks = sorted((cs + o, cs + o + n) for o, n in zip(offs, sizes))
            if strand == "-":
                blocks = blocks[::-1]

            attrs = f'gene_id "{gene_id}"; transcript_id "{gene_id}";'
            cumulative = 0
            for bs, be in blocks:
                phase = (3 - cumulative % 3) % 3
                out.write("\t".join([
                    chrom, args.source, "CDS", str(bs + 1), str(be), ".", strand, str(phase), attrs
                ]) + "\n")
                cumulative += be - bs
            n_out += 1

    print(f"read        {n_in:,} ORFs", file=sys.stderr)
    print(f"written     {n_out:,} ORFs as CDS features -> {args.out_gtf}", file=sys.stderr)
    if n_skip:
        print(f"skipped     {n_skip:,} on scaffolds -> {args.out_skipped}", file=sys.stderr)
        print(f"            (the Annotator drops these; they rejoin as "
              f"'not_annotated_scaffold')", file=sys.stderr)
    if n_preexisting:
        print(f"  WARNING: {n_preexisting:,} ORF ids already contained a chrN:start-end\n"
              f"           pattern. The Annotator reads the FIRST match, which is now\n"
              f"           the appended one, so this is handled -- but check the ids.",
              file=sys.stderr)
    if n_out == 0:
        print("ERROR: nothing written; every ORF was on an unhandled seqname.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
