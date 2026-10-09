#!/usr/bin/env python3
"""Translation signature scores for Ribo-seq ORFs.

Implements PIF, Uniformity and Drop-off as defined in Chothani et al.,
Nucleic Acids Research 2026;54:gkag234 (and Chothani et al., Mol Cell 2022;
Menon, Rackham & Chothani, Methods Mol Biol 2026;2992:91-112).

Inputs are all genomic:
  - ORF catalogue BED12
  - P-site bedGraphs      (fwd + rev)

Scores are computed per ORF against pooled P-sites. Thresholds are not
hardcoded: they are derived by fitting a normal to the score distribution of
annotated CDS in the same data. See --threshold-tail.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass

import io
import numpy as np
import pandas as pd
from scipy import stats

IN_FRAME = 0  # the ORF's own reading frame, "Frame 1" in the paper's wording
UNIFORMITY_CODON_FRAC = 1.0 / 3.0  # paper: codon is translated if >33.33% in frame
DROPOFF_FLANK_NT = 15  # 15 nt each side of the stop codon = 30 nt window


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# --------------------------------------------------------------------------
# P-site tracks
# --------------------------------------------------------------------------


def _pool_intervals(starts: np.ndarray, ends: np.ndarray, values: np.ndarray):
    """Sum possibly overlapping [start, end) intervals into disjoint ones (value > 0)."""
    pts = np.concatenate([starts, ends])
    dv = np.concatenate([values, -values])
    order = np.argsort(pts, kind="stable")
    pts, dv = pts[order], dv[order]
    uniq, first = np.unique(pts, return_index=True)
    level = np.round(np.cumsum(np.add.reduceat(dv, first)), 9)
    keep = level[:-1] > 0
    return uniq[:-1][keep], uniq[1:][keep], level[:-1][keep]


class PsiteIndex:
    """Per-strand, per-chromosome P-site counts from bedGraph, queried by position.

    bedGraph intervals are half-open [start, end). We keep sorted start arrays
    per chromosome and resolve positions with searchsorted, which avoids
    materialising a genome-length array.
    """

    def __init__(self) -> None:
        self._idx: dict[tuple[str, str], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    @classmethod
    def from_bedgraphs(cls, forward: list[str], reverse: list[str]) -> "PsiteIndex":
        self = cls()
        frames = {"+": [], "-": []}
        for strand, paths in (("+", forward), ("-", reverse)):
            for path in paths:
                # plastid's make_wiggle writes ~21 lines of '##' provenance
                # headers BEFORE the UCSC 'track' line. The old comment="t"
                # skipped neither -- and worse, it made "t" a comment character,
                # truncating every line at its first "t". Filter explicitly.
                with _open(path) as fh:
                    body = io.StringIO("".join(
                        ln for ln in fh
                        if not ln.startswith(("#", "track", "browser")) and ln.strip()
                    ))
                    df = pd.read_csv(
                        body,
                        sep="\t",
                        header=None,
                        names=["chrom", "start", "end", "value"],
                        dtype={"chrom": str, "start": np.int64, "end": np.int64, "value": np.float64},
                    )
                # plastid writes minus-strand coverage as negative values
                df["value"] = df["value"].abs()
                frames[strand].append(df)

        for strand, dfs in frames.items():
            if not dfs:
                continue
            allrows = pd.concat(dfs, ignore_index=True)
            # Pool replicates exactly: plastid merges runs of equal values into wider
            # intervals, so libraries' intervals overlap without matching. Sum them over
            # every breakpoint (a sweep), giving disjoint intervals counts_at can read.
            for chrom, sub in allrows.groupby("chrom", sort=False):
                self._idx[(chrom, strand)] = _pool_intervals(
                    sub["start"].to_numpy(np.int64), sub["end"].to_numpy(np.int64),
                    sub["value"].to_numpy(np.float64))
        return self

    def counts_at(self, chrom: str, strand: str, positions: np.ndarray) -> np.ndarray:
        """P-site count at each 0-based genomic position. Zero where uncovered."""
        key = (chrom, strand)
        out = np.zeros(positions.size, dtype=np.float64)
        if key not in self._idx or positions.size == 0:
            return out
        starts, ends, values = self._idx[key]
        # Index of the last interval starting at or before each position.
        i = np.searchsorted(starts, positions, side="right") - 1
        valid = (i >= 0) & (positions < ends[np.clip(i, 0, None)])
        out[valid] = values[i[valid]]
        return out


# --------------------------------------------------------------------------
# ORF models
# --------------------------------------------------------------------------


@dataclass(slots=True)
class Orf:
    name: str
    chrom: str
    strand: str
    positions: np.ndarray  # 0-based genomic positions, 5'->3' in transcript order
    extra: dict

    @property
    def n_codons(self) -> int:
        return self.positions.size // 3


def parse_bed12(path: str, keep_extra_cols: bool = True):
    """Yield Orf objects from a BED12 file, blocks spliced in transcript order."""
    with _open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 12:
                raise ValueError(f"{path}:{lineno}: expected 12+ BED columns, got {len(f)}")

            chrom, chrom_start, name, strand = f[0], int(f[1]), f[3], f[5]
            block_sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            block_starts = [int(x) for x in f[11].rstrip(",").split(",")]

            pos = np.concatenate(
                [np.arange(chrom_start + s, chrom_start + s + n, dtype=np.int64)
                 for s, n in zip(block_starts, block_sizes)]
            )
            pos.sort()
            if strand == "-":
                pos = pos[::-1]  # transcript order is descending in genomic coords

            extra = {"score": f[4], "thick_start": f[6], "thick_end": f[7]} if keep_extra_cols else {}
            if len(f) > 12:
                extra["extra_fields"] = "\t".join(f[12:])
            yield Orf(name=name, chrom=chrom, strand=strand, positions=pos, extra=extra)


def cds_from_gtf(path: str, one_per_gene: bool = True):
    """Yield Orf objects for annotated CDS, used to calibrate thresholds.

    CDS segments of a transcript are concatenated in transcript order, WITH the stop
    codon (Ensembl stop_codon rows), because catalogue ORFs include theirs: drop-off's
    window is placed from the ORF end. Incomplete CDS (cds_start_NF / cds_end_NF tags,
    no stop codon, first-row phase != 0, or length not a multiple of 3) are skipped:
    their frame or ends are unknown. With one_per_gene, the transcript with the
    longest CDS wins, which approximates the paper's representative transcript.
    """
    by_tx: dict[str, dict] = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("CDS", "stop_codon"):
                continue
            attrs = f[8]
            tx = _gtf_attr(attrs, "transcript_id")
            gene = _gtf_attr(attrs, "gene_id")
            if tx is None or gene is None:
                continue
            rec = by_tx.setdefault(tx, {"chrom": f[0], "strand": f[6], "gene": gene, "segs": [],
                                        "stop": [], "incomplete": False, "phase0": None})
            if "cds_start_NF" in attrs or "cds_end_NF" in attrs:
                rec["incomplete"] = True
            iv = (int(f[3]) - 1, int(f[4]))  # GTF is 1-based inclusive
            if f[2] == "stop_codon":
                rec["stop"].append(iv)
            else:
                rec["segs"].append(iv)
                first = rec["phase0"] is None or (iv[0] < rec["phase0"][0] if f[6] == "+"
                                                  else iv[0] > rec["phase0"][0])
                if first:
                    rec["phase0"] = (iv[0], f[7])

    for tx in [t for t, r in by_tx.items()
               if r["incomplete"] or not r["segs"] or sum(e - s for s, e in r["stop"]) != 3
               or r["phase0"][1] not in ("0", ".") or sum(e - s for s, e in r["segs"]) % 3]:
        del by_tx[tx]
    for rec in by_tx.values():
        rec["segs"] = rec["segs"] + rec["stop"]

    best: dict[str, tuple[int, str]] = {}
    for tx, rec in by_tx.items():
        length = sum(e - s for s, e in rec["segs"])
        if not one_per_gene:
            best[tx] = (length, tx)
        else:
            g = rec["gene"]
            if g not in best or length > best[g][0]:
                best[g] = (length, tx)

    for _, tx in best.values():
        rec = by_tx[tx]
        pos = np.concatenate([np.arange(s, e, dtype=np.int64) for s, e in sorted(rec["segs"])])
        if rec["strand"] == "-":
            pos = pos[::-1]
        yield Orf(name=tx, chrom=rec["chrom"], strand=rec["strand"], positions=pos,
                  extra={"gene_id": rec["gene"]})


def _merge(segs):
    segs = sorted(segs)
    merged = [list(segs[0])]
    for a, b in segs[1:]:
        if a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    m = np.array(merged, dtype=np.int64)
    return m[:, 0], m[:, 1]


def _in(ivs, positions):
    starts, ends = ivs
    i = np.searchsorted(starts, positions, side="right") - 1
    return (i >= 0) & (positions < ends[np.clip(i, 0, None)])


class CdsMask:
    """Every annotated CDS (all transcripts, stop codons included), per strand, to find
    the nucleotides an ORF does not share with any known protein.

    frame_aware (default): a nucleotide counts as shared only where some CDS reads it in
    the ORF's own frame. Each CDS segment is stored with its genomic frame key, which
    is constant along a segment: (g - codon position) % 3 on +, (g + codon position) % 3
    on -. An ORF nucleotide at g with codon position c is shared if a CDS segment with
    the same key covers g. frame_aware=False: any CDS nucleotide counts (frame-blind)."""

    def __init__(self, path: str, frame_aware: bool = True):
        self.frame_aware = frame_aware
        by_tx: dict[str, dict] = {}
        with _open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 8)
                if len(f) < 9 or f[2] not in ("CDS", "stop_codon"):
                    continue
                tx = _gtf_attr(f[8], "transcript_id") or f"{f[0]}:{f[3]}"
                rec = by_tx.setdefault(tx, {"key": (f[0], f[6]), "rows": []})
                rec["rows"].append((int(f[3]) - 1, int(f[4]), f[2], f[7]))
        blind: dict = {}
        framed: dict = {}
        for rec in by_tx.values():
            chrom, strand = rec["key"]
            rows = sorted(rec["rows"], key=lambda r: r[0], reverse=(strand == "-"))
            cds = [r for r in rows if r[2] == "CDS"]
            ph = cds[0][3] if cds else "0"
            c = (3 - int(ph)) % 3 if ph in ("0", "1", "2") else 0   # codon position of the first nt
            for s, e, _, _ in rows:                               # CDS then stop, transcript order
                blind.setdefault(rec["key"], []).append((s, e))
                k = (s - c) % 3 if strand == "+" else (e - 1 + c) % 3
                framed.setdefault((chrom, strand, k), []).append((s, e))
                c = (c + e - s) % 3
        self._blind = {k: _merge(v) for k, v in blind.items()}
        self._framed = {k: _merge(v) for k, v in framed.items()}

    def covered(self, chrom: str, strand: str, positions: np.ndarray) -> np.ndarray:
        """Frame-blind: positions any annotated CDS covers."""
        ivs = self._blind.get((chrom, strand))
        return _in(ivs, positions) if ivs is not None else np.zeros(positions.size, dtype=bool)

    def shared(self, chrom: str, strand: str, positions: np.ndarray) -> np.ndarray:
        """ORF positions (transcript order) whose nucleotide an annotated CDS reads in
        the ORF's frame (frame_aware), or covers at all (frame-blind)."""
        if not self.frame_aware:
            return self.covered(chrom, strand, positions)
        cp = np.arange(positions.size) % 3
        key = (positions - cp) % 3 if strand == "+" else (positions + cp) % 3
        out = np.zeros(positions.size, dtype=bool)
        for k in (0, 1, 2):
            ivs = self._framed.get((chrom, strand, k))
            m = key == k
            if ivs is not None and m.any():
                out[m] = _in(ivs, positions[m])
        return out


def _gtf_attr(attrs: str, key: str) -> str | None:
    tag = f'{key} "'
    i = attrs.find(tag)
    if i < 0:
        return None
    i += len(tag)
    j = attrs.find('"', i)
    return attrs[i:j] if j > 0 else None


# --------------------------------------------------------------------------
# The three scores
# --------------------------------------------------------------------------


def score_orf(orf: Orf, psites: PsiteIndex, uniformity_zero_codons: str = "fail",
              cds_mask: CdsMask | None = None, psites_unique: PsiteIndex | None = None) -> dict:
    """PIF, Uniformity and Drop-off for one ORF.

    counts[i] is the P-site count at the i-th nucleotide of the ORF in
    transcript order, so frame = i % 3 and frame 0 is the ORF's own frame.
    """
    counts = psites.counts_at(orf.chrom, orf.strand, orf.positions)
    total = counts.sum()
    n_codons = orf.n_codons

    result = {
        "orf_id": orf.name,
        "chrom": orf.chrom,
        "strand": orf.strand,
        "length_nt": int(orf.positions.size),
        "n_codons": int(n_codons),
        "psites_total": float(total),
        "psites_inframe": float(counts[IN_FRAME::3].sum()),
        "pif": np.nan,
        "uniformity": np.nan,
        "dropoff": np.nan,
    }
    if cds_mask is not None:
        # Evidence of its own: P-sites in nucleotides no annotated CDS covers (an extra
        # exon, an N-extension, an intergenic ORF's whole length). An isoform can pass
        # the whole-ORF scores on the P-sites it shares with a highly translated CDS.
        own = ~cds_mask.shared(orf.chrom, orf.strand, orf.positions)
        frame0 = (np.arange(orf.positions.size) % 3) == IN_FRAME
        result["unshared_nt"] = int(own.sum())
        result["unshared_psites"] = float(counts[own].sum())
        result["unshared_inframe_psites"] = float(counts[own & frame0].sum())
    if psites_unique is not None:
        # uniquely mapped reads only: a pseudogene's parent-gene reads drop out here
        cu = psites_unique.counts_at(orf.chrom, orf.strand, orf.positions)
        result["psites_unique"] = float(cu.sum())
        result["psites_inframe_unique"] = float(cu[IN_FRAME::3].sum())
        if cds_mask is not None:
            result["unshared_inframe_psites_unique"] = float(cu[own & frame0].sum())
        # Density for the above-CDS-median flag (final_outputs): in-frame P-sites over the body
        # codons, between start and stop, so initiation / termination peaks don't inflate short
        # ORFs. free_*: only body codons no annotated CDS touches in any frame (host-free).
        nb = max(n_codons - 2, 0)
        body = cu[3 + IN_FRAME:3 * (n_codons - 1):3] if nb else np.zeros(0)
        result["body_codons"] = int(nb)
        result["body_inframe_psites_unique"] = float(body.sum())
        if cds_mask is not None:
            hit = cds_mask.covered(orf.chrom, orf.strand, orf.positions)[:3 * n_codons]
            free = ~hit.reshape(n_codons, 3).any(axis=1)[1:n_codons - 1] if nb else np.zeros(0, bool)
            result["free_body_codons"] = int(free.sum())
            result["free_body_inframe_psites_unique"] = float(body[free].sum())
    if uniformity_zero_codons not in ("fail", "ignore"):
        raise ValueError(f"uniformity_zero_codons must be 'fail' or 'ignore', got {uniformity_zero_codons!r}")

    # PIF and Uniformity are 0/0 with no P-sites in the ORF, so they stay NaN.
    # Drop-off is computed below regardless: it only needs the stop-codon window,
    # and 0.0 (reads only after the stop) is a real, informative result.
    if total > 0 and n_codons > 0:
        # --- PIF: in-frame P-sites over all P-sites in the ORF
        result["pif"] = float(counts[IN_FRAME::3].sum() / total)

        # --- Uniformity: codons whose in-frame share exceeds 33.33%
        usable = n_codons * 3
        by_codon = counts[:usable].reshape(n_codons, 3)
        codon_totals = by_codon.sum(axis=1)
        covered = codon_totals > 0
        frac_inframe = np.zeros(n_codons)
        frac_inframe[covered] = by_codon[covered, IN_FRAME] / codon_totals[covered]
        translated = covered & (frac_inframe > UNIFORMITY_CODON_FRAC)

        # zero-coverage codons: 'fail' counts them against the ORF (default),
        # 'ignore' drops them from the denominator
        denom = int(covered.sum()) if uniformity_zero_codons == "ignore" else n_codons
        result["uniformity"] = float(translated.sum() / denom) if denom else np.nan

    if orf.positions.size < 3:
        return result

    # --- Drop-off: in-frame P-sites 15 nt before vs after the stop codon.
    # The ORF's last 3 nt are the stop codon; "before" is the 15 nt preceding
    # it, "after" is the 15 nt of downstream transcript sequence. Downstream
    # positions are extrapolated along the genome, which is exact unless the
    # stop sits within 15 nt of a splice junction.
    before_pos = orf.positions[-(DROPOFF_FLANK_NT + 3):-3]
    step = 1 if orf.strand == "+" else -1
    last = int(orf.positions[-1])
    after_pos = np.arange(last + step, last + step * (DROPOFF_FLANK_NT + 1), step, dtype=np.int64)

    # Keep frame relative to the ORF so "in-frame" means the same thing on both sides.
    before_counts = psites.counts_at(orf.chrom, orf.strand, before_pos)
    after_counts = psites.counts_at(orf.chrom, orf.strand, after_pos)
    n_before = orf.positions.size - 3 - before_pos.size  # index of before_pos[0] in the ORF
    before_in = before_counts[(IN_FRAME - n_before) % 3::3].sum()
    after_in = after_counts[(IN_FRAME - orf.positions.size) % 3::3].sum()

    denom = before_in + after_in
    result["dropoff"] = float(before_in / denom) if denom > 0 else np.nan
    return result


def score_many(orfs, psites: PsiteIndex, uniformity_zero_codons: str = "fail",
               cds_mask: CdsMask | None = None, psites_unique: PsiteIndex | None = None) -> pd.DataFrame:
    rows = [score_orf(o, psites, uniformity_zero_codons, cds_mask, psites_unique) for o in orfs]
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# PPM and thresholds
# --------------------------------------------------------------------------


def add_ppm(df: pd.DataFrame) -> pd.DataFrame:
    """P-sites per million: length-normalise, then scale to per-million."""
    density = df["psites_total"] / df["length_nt"].replace(0, np.nan)
    total = density.sum(skipna=True)
    df["ppm"] = density / total * 1e6 if total and total > 0 else np.nan
    return df


def calibrate(cds_scores: pd.DataFrame, tail: str = "lower", percentile: float = 0.95,
              min_psites: float = 100.0, method: str = "empirical") -> dict:
    """Derive thresholds from the annotated-CDS score distribution.

    PIF and Uniformity: the `percentile` point of the CDS distribution.
    Drop-off: the CDS mean, per the paper.

    Two things here were measured, not assumed (45M-read library):

    1. min_psites=100, not 1. With min_psites=1 the CDS set keeps thousands of
       barely-covered CDS that score exactly 0, producing a spike at zero:
           all CDS (n=22,178)        PIF p5 = 0.000   median = 0.711
           psites_total >= 100       PIF p5 = 0.497   median = 0.712
       The MEDIAN does not move -- only the zero spike is removed. Chothani et
       al. never hit this because their pooled data is 1.3e9 P-sites and covers
       essentially every CDS; on a single library the filter has to be explicit.
       At >=100 the derived thresholds (PIF 0.497, Uniformity 0.525) land in the
       paper's published range (RiboCrypt PIF 0.5112, body-map Uniformity 0.713).

    2. method='empirical', not a Gaussian fit. PIF/Uniformity/Drop-off are
       bounded [0,1] and left-skewed, so stats.norm.fit is invalid:
           mu + 1.645*sigma = 1.062 (PIF), 1.050 (Uniformity)  -- impossible.
       That impossibility is what originally motivated tail='lower'; it was really a symptom of the bad fit. With an empirical
       quantile on a coverage-filtered set, the tail choice stops being load
       bearing. method='normal' reproduces the old behaviour for comparison.

    tail='lower' puts the threshold where `percentile` of CDS score ABOVE it.
    """
    usable = cds_scores[cds_scores["psites_total"] >= min_psites]
    if usable.empty:
        raise ValueError(
            f"No annotated CDS had >= {min_psites} P-sites; cannot calibrate. "
            "Library may be too shallow -- lower min_psites and treat the "
            "resulting thresholds with suspicion."
        )

    q = (1.0 - percentile) if tail == "lower" else percentile
    thresholds = {"n_cds_used": int(len(usable)), "tail": tail,
                  "percentile": percentile, "min_psites": float(min_psites),
                  "method": method}

    for metric in ("pif", "uniformity"):
        vals = usable[metric].dropna().to_numpy()
        if vals.size < 2:
            raise ValueError(f"Too few CDS values to calibrate {metric}.")
        if method == "normal":
            mu, sigma = stats.norm.fit(vals)
            thr = float(np.clip(stats.norm.ppf(q, mu, sigma), 0.0, 1.0))
        else:
            mu, sigma = float(np.mean(vals)), float(np.std(vals))
            thr = float(np.clip(np.quantile(vals, q), 0.0, 1.0))
        thresholds[metric] = thr
        thresholds[f"{metric}_mu"] = float(mu)
        thresholds[f"{metric}_sigma"] = float(sigma)
        thresholds[f"{metric}_median"] = float(np.median(vals))

    vals = usable["dropoff"].dropna().to_numpy()
    if vals.size >= 2:
        mu, sigma = float(np.mean(vals)), float(np.std(vals))
    else:
        mu = sigma = float("nan")
    thresholds["dropoff"] = float(np.clip(mu, 0.0, 1.0))  # paper uses the mean
    thresholds["dropoff_mu"] = mu
    thresholds["dropoff_sigma"] = sigma
    thresholds["n_cds_dropoff"] = int(vals.size)
    return thresholds


def cds_median_density(cds: pd.DataFrame, min_psites: float) -> dict:
    """Reference for final_outputs' above-CDS-median flag, from CDS scored with unique P-sites:
    the median over calibration CDS (>= min_psites pooled P-sites) of in-frame unique P-sites per
    body codon, per million in-frame unique P-sites in all CDS (the run's depth)."""
    tot = float(cds["psites_inframe_unique"].sum())
    cal = cds[(cds["psites_total"] >= min_psites) & (cds["body_codons"] > 0)]
    dens = cal["body_inframe_psites_unique"] / cal["body_codons"] / tot * 1e6 if tot > 0 else pd.Series(dtype=float)
    return {"cds_median_density_unique": float(dens.median()) if len(dens) else np.nan,
            "inframe_cds_psites_unique": tot, "n_cds_density": int(len(cal))}


def apply_unshared(orfs: pd.DataFrame, mode: str, min_inframe: int) -> int:
    """pass_unshared: no unshared sequence, or >= min_inframe in-frame P-sites in it.
    mode 'drop' also sets exclude_reason and clears pass_all. Returns ORFs failing it."""
    orfs["pass_unshared"] = (orfs["unshared_nt"] == 0) | (orfs["unshared_inframe_psites"] >= min_inframe)
    orfs["exclude_reason"] = ""
    if mode == "drop":
        orfs.loc[~orfs["pass_unshared"], "exclude_reason"] = (
            f"fewer than {min_inframe} in-frame P-sites outside annotated CDS")
        orfs["pass_all"] = orfs["pass_all"] & orfs["pass_unshared"]
    return int((~orfs["pass_unshared"]).sum())


def apply_thresholds(df: pd.DataFrame, thr: dict, require: list[str]) -> pd.DataFrame:
    for metric in ("pif", "uniformity", "dropoff"):
        df[f"pass_{metric}"] = df[metric] >= thr[metric]
    df["pass_all"] = np.logical_and.reduce([df[f"pass_{m}"] for m in require])
    return df


# --------------------------------------------------------------------------


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orf-bed12", required=True, help="ORF catalogue in BED12 (genomic).")
    p.add_argument("--gtf", required=True, help="GENCODE GTF, for CDS threshold calibration.")
    p.add_argument("--psite-forward", required=True, nargs="+", help="Forward-strand P-site bedGraph(s).")
    p.add_argument("--psite-reverse", required=True, nargs="+", help="Reverse-strand P-site bedGraph(s).")
    p.add_argument("--psite-forward-unique", nargs="+", help="Forward P-sites, uniquely mapped reads only.")
    p.add_argument("--psite-reverse-unique", nargs="+", help="Reverse P-sites, uniquely mapped reads only.")
    p.add_argument("--out-orfs", required=True, help="Per-ORF scores TSV.")
    p.add_argument("--out-cds", required=True, help="Per-CDS calibration scores TSV.")
    p.add_argument("--out-thresholds", required=True, help="Derived thresholds TSV.")
    p.add_argument("--threshold-tail", choices=("lower", "upper"), default="lower")
    p.add_argument("--threshold-percentile", type=float, default=0.95)
    p.add_argument("--threshold-method", choices=("empirical", "normal"), default="empirical",
                   help="empirical quantile (default) or legacy Gaussian fit.")
    p.add_argument("--calibration-min-psites", type=float, default=100.0,
                   help="Minimum P-sites for a CDS to enter threshold calibration.")
    p.add_argument("--uniformity-zero-codons", choices=("fail", "ignore"), default="fail")
    p.add_argument("--require", nargs="+", default=["pif", "uniformity", "dropoff"],
                   choices=("pif", "uniformity", "dropoff"),
                   help="Metrics an ORF must pass to enter the primary set.")
    p.add_argument("--min-ppm", type=float, default=None,
                   help="If set, also require PPM above this (paper uses 1).")
    p.add_argument("--unshared-evidence", choices=("drop", "report"), default="report",
                   help="drop: an ORF with nucleotides no annotated CDS covers, but no in-frame "
                        "P-site there, leaves the primary set and gets exclude_reason "
                        "(final_outputs moves it out of Comprehensive too). report: columns only.")
    p.add_argument("--unshared-min-inframe", type=int, default=1,
                   help="in-frame P-sites an ORF needs in its unshared sequence (default 1)")
    p.add_argument("--unshared-frame", choices=("aware", "blind"), default="aware",
                   help="aware (default): CDS nucleotides count as shared only in the ORF's own "
                        "frame (an other-frame ORF is a different protein). blind: any CDS nt.")
    args = p.parse_args(argv)

    print("Loading P-site tracks...", file=sys.stderr)
    psites = PsiteIndex.from_bedgraphs(args.psite_forward, args.psite_reverse)
    uniq = (PsiteIndex.from_bedgraphs(args.psite_forward_unique, args.psite_reverse_unique)
            if args.psite_forward_unique and args.psite_reverse_unique else None)

    print("Scoring annotated CDS for calibration...", file=sys.stderr)
    cds = add_ppm(score_many(cds_from_gtf(args.gtf), psites, args.uniformity_zero_codons, psites_unique=uniq))
    cds.to_csv(args.out_cds, sep="\t", index=False)

    thr = calibrate(cds, tail=args.threshold_tail, percentile=args.threshold_percentile,
                    min_psites=args.calibration_min_psites, method=args.threshold_method)
    if uniq is not None:
        thr.update(cds_median_density(cds, args.calibration_min_psites))
    pd.DataFrame([thr]).to_csv(args.out_thresholds, sep="\t", index=False)
    print(
        f"  thresholds from {thr['n_cds_used']} CDS "
        f"(tail={thr['tail']}): PIF={thr['pif']:.4f} "
        f"Uniformity={thr['uniformity']:.4f} Drop-off={thr['dropoff']:.4f}",
        file=sys.stderr,
    )

    print("Scoring called ORFs...", file=sys.stderr)
    mask = CdsMask(args.gtf, frame_aware=(args.unshared_frame == "aware"))
    orfs = add_ppm(score_many(parse_bed12(args.orf_bed12), psites, args.uniformity_zero_codons, mask, uniq))
    if uniq is not None:
        orfs["ppm_unique"] = add_ppm(orfs[["psites_unique", "length_nt"]].rename(
            columns={"psites_unique": "psites_total"}))["ppm"]
    orfs = apply_thresholds(orfs, thr, args.require)
    if args.min_ppm is not None:
        orfs["pass_ppm"] = orfs["ppm"] >= args.min_ppm
        orfs["pass_all"] = orfs["pass_all"] & orfs["pass_ppm"]
    n_out = apply_unshared(orfs, args.unshared_evidence, args.unshared_min_inframe)
    if args.unshared_evidence == "drop":
        print(f"  {n_out} ORFs excluded: fewer than {args.unshared_min_inframe} in-frame P-sites in "
              "the sequence they don't share with an annotated CDS", file=sys.stderr)
    orfs.to_csv(args.out_orfs, sep="\t", index=False)

    n_pass = int(orfs["pass_all"].sum())
    print(f"  {n_pass} / {len(orfs)} ORFs in the primary set", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
