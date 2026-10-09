#!/usr/bin/env python3
"""Score distributions: annotated CDS vs called ORFs, with derived thresholds.

This is the plot to look at first. If the called-ORF distribution sits far left
of the CDS distribution, the thresholds are doing real work. If the CDS
distribution itself is broad or low, the libraries are the problem, not the ORFs.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

METRICS = [("pif", "PIF"), ("uniformity", "Uniformity"), ("dropoff", "Drop-off")]

orfs = pd.read_csv(snakemake.input.orfs, sep="\t")
cds = pd.read_csv(snakemake.input.cds, sep="\t")
thr = pd.read_csv(snakemake.input.thresholds, sep="\t").iloc[0]

fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)

for ax, (col, label) in zip(axes, METRICS):
    c = cds[col].dropna()
    o = orfs[col].dropna()
    bins = 50
    ax.hist(c, bins=bins, range=(0, 1), density=True, alpha=0.55,
            label=f"annotated CDS (n={len(c):,})", color="#4C72B0")
    ax.hist(o, bins=bins, range=(0, 1), density=True, alpha=0.55,
            label=f"called ORFs (n={len(o):,})", color="#DD8452")
    ax.axvline(thr[col], color="k", ls="--", lw=1.4,
               label=f"threshold = {thr[col]:.3f}")
    ax.set_xlabel(label)
    ax.set_ylabel("density")
    ax.legend(fontsize=7, frameon=False)

n_pass = int(orfs["pass_all"].sum())
fig.suptitle(
    f"Translation signature scores  |  primary set: {n_pass:,} / {len(orfs):,} ORFs  "
    f"|  calibrated on {int(thr['n_cds_used']):,} CDS, tail={thr['tail']}",
    fontsize=10,
)
fig.savefig(snakemake.output[0], dpi=160)
