"""collapse_catalogue.py: MMseqs2 clusters fold only within one locus. Prints PASS/FAIL,
exits 1 on failure."""
import os
import subprocess
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "workflow", "scripts")
sys.path.insert(0, SCRIPTS)
from collapse_catalogue import same_locus_groups  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


cols = ["orf_id", "chrom", "start", "end", "strand", "orf_class", "aa_length", "is_smorf",
        "called_by_price", "called_by_ribocode", "pval_price", "pval_ribocode", "n_callers",
        "n_samples", "samples", "orf_type_native"]
rows = [
    # parent gene (chr3) and its processed pseudogene (chr7): same protein, two loci
    ["parent", "chr3", "1000", "1339", "-", "canonical_cds", "112", "1", "0", "1", "", "0.01", "1", "1", "s1", "a"],
    ["pseudo", "chr7", "5000", "5351", "-", "novel_u", "116", "1", "0", "1", "", "0.02", "1", "1", "s1", "b"],
    # two near-identical smORFs at one locus
    ["loc_a", "chr1", "100", "190", "+", "novel_u", "29", "1", "1", "0", "0.001", "", "1", "1", "s1", "c"],
    ["loc_b", "chr1", "103", "190", "+", "novel_u", "28", "1", "0", "1", "", "0.03", "1", "1", "s2", "d"],
    # same locus and cluster, but spliced (intron 150-200 inside the shared span)
    ["iso", "chr1", "100", "290", "+", "novel_u", "46", "1", "0", "1", "", "0.04", "1", "1", "s1", "e"],
]
rows += [   # a second foldable pair, so the fold loop runs more than once
    ["loc_c", "chr2", "500", "590", "+", "novel_u", "29", "1", "1", "0", "0.001", "", "1", "1", "s1", "f"],
    ["loc_d", "chr2", "503", "590", "+", "novel_u", "28", "1", "0", "1", "", "0.02", "1", "1", "s2", "g"],
]
BLOCKS = {"iso": ((0, 50), (100, 190))}   # BED12 relative blocks: 100-150 and 200-290
df = pd.DataFrame(rows, columns=cols)
check("cross-locus members split", sorted(sorted(df.loc[g, "orf_id"]) for g in same_locus_groups(df, [0, 1])),
      [["parent"], ["pseudo"]])
check("same-locus members stay together", [sorted(df.loc[g, "orf_id"]) for g in same_locus_groups(df, [2, 3, 4])],
      [["iso", "loc_a", "loc_b"]])

with tempfile.TemporaryDirectory() as d:
    tsv, bed, clu, out = (os.path.join(d, x) for x in ("c.tsv", "c.bed12", "clu.tsv", "out"))
    df.to_csv(tsv, sep="\t", index=False)
    with open(bed, "w") as fh:
        for r in rows:
            b = BLOCKS.get(r[0], ((0, int(r[3]) - int(r[2])),))
            sizes = ",".join(str(e - st) for st, e in b) + ","
            starts = ",".join(str(st) for st, _ in b) + ","
            fh.write(f"{r[1]}\t{r[2]}\t{r[3]}\t{r[0]}\t0\t{r[4]}\t{r[2]}\t{r[3]}\t0\t{len(b)}\t{sizes}\t{starts}\n")
    with open(clu, "w") as fh:
        fh.write("parent(-)\tparent(-)\nparent(-)\tpseudo(-)\nloc_a(+)\tloc_a(+)\nloc_a(+)\tloc_b(+)\n"
                 "loc_a(+)\tiso(+)\nloc_c(+)\tloc_c(+)\nloc_c(+)\tloc_d(+)\n")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "collapse_catalogue.py"), "--tsv", tsv,
                        "--bed12", bed, "--clusters", clu, "--smorf-max-aa", "150", "--out-prefix", out],
                       capture_output=True, text=True)
    check("exit 0", r.returncode, 0)
    o = pd.read_csv(out + ".tsv", sep="\t", dtype=str, keep_default_na=False).set_index("orf_id")
    check("pseudogene and spliced isoform keep their own rows", sorted(o.index), ["iso", "loc_a", "loc_c", "parent", "pseudo"])
    check("second group folds too", o.at["loc_c", "collapsed_orf_ids"], "loc_d")
    check("same-locus pair folded, traceable", o.at["loc_a", "collapsed_orf_ids"], "loc_b")
    check("callers OR-ed on fold", (o.at["loc_a", "called_by_price"], o.at["loc_a", "called_by_ribocode"]), ("1", "1"))

sys.exit(1 if failed else 0)
