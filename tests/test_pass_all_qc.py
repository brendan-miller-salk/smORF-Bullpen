"""pass_all_qc checks against Chothani et al. Prints PASS/FAIL, exits 1 on failure."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workflow", "scripts"))
from pass_all_qc import run  # noqa: E402

failed = 0


def check(name, got, want):
    global failed
    ok = got == want
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got={got} want={want}")


def status(rows):
    return {r["check"]: r["status"] for r in rows}


thr_in = pd.DataFrame([dict(pif=0.6, uniformity=0.8, dropoff=0.9, tail="lower",
                            percentile=0.95, min_psites=100)])
thr_out = thr_in.assign(pif=0.3, uniformity=0.5)
cds = pd.DataFrame(dict(psites_total=[200] * 20, pif=[0.59] + [0.7] * 19,
                        uniformity=[0.79] + [0.9] * 19))

# Paper-like: 36% primary, short ORFs pass PIF more, overlapping ORFs pass less.
n = 100
scores = pd.DataFrame(dict(
    orf_id=[f"o{i}" for i in range(n)],
    orf_class="novel_u",
    length_nt=[60 if i < 50 else 300 for i in range(n)],
    pass_pif=[str(i < 40 or 50 <= i < 70) for i in range(n)],
    pass_all=[str(i < 36) for i in range(n)],
    smorf_type=["uORF" if i < 50 else "iORF" for i in range(n)],
))
st = status(run(scores, thr_in, cds))
for c in ("threshold_pif", "threshold_uniformity", "threshold_dropoff", "calibration_pif",
          "calibration_uniformity", "primary_fraction", "length_bias", "overlap_bias"):
    check(f"paper-like {c}", st.get(c), "OK")

bad = scores.assign(length_nt=[300 if i < 50 else 60 for i in range(n)],
                    pass_all=[str(i < 90) for i in range(n)],
                    smorf_type=["iORF" if i < 50 else "uORF" for i in range(n)])
st = status(run(bad, thr_out, cds))
for c in ("threshold_pif", "threshold_uniformity", "length_bias", "primary_fraction", "overlap_bias"):
    check(f"off-paper {c}", st.get(c), "WARN")
check("off-paper dropoff still OK", st.get("threshold_dropoff"), "OK")

# orf_set restricts to the Comprehensive set; canonical CDS are excluded.
tagged = scores.assign(orf_set=["primary" if i < 36 else "not_supported" if i >= 90 else
                                "comprehensive" for i in range(n)],
                       orf_class=["canonical_cds" if i == 0 else "novel_u" for i in range(n)])
rows = {r["check"]: r for r in run(tagged, thr_in, cds)}
check("orf_set denominators", rows["primary_fraction"]["observed"], "0.393 (35/89)")

no_cls = scores.drop(columns=["orf_class"])
check("missing orf_class -> SKIP", status(run(no_cls, thr_in, cds)).get("primary_fraction"), "SKIP")

print(f"\n{'ALL PASS' if not failed else f'{failed} FAILED'}")
sys.exit(1 if failed else 0)
