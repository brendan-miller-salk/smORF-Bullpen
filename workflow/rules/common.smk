import csv
import re
import os

# Provided mode (config `orf_gtf`): annotate the GTF's ORFs; no Ribo-seq, callers,
# scoring or counting, and no samplesheet.
PROVIDED = bool(config.get("orf_gtf"))
# blastoff classes exempt from the own-sequence and RPKM rules: Reviewed, plus exact TrEMBL
# matches with trembl_identical_exempt (they must still pass the scores)
EXEMPT_CLASSES = ["sp_identical"] + (["tr_identical"] if config.get("trembl_identical_exempt") else [])
REQUIRED = ["outdir", "genome", "gtf", "index_dir", "containers"] + ([] if PROVIDED else ["samplesheet", "sortmerna_dbs"])
_missing = [k for k in REQUIRED if k not in config]
if _missing:
    raise WorkflowError(
        f"Config is missing {_missing}. Pass the run's config with --configfile "
        f"(template: {workflow.basedir}/../config/config.template.yaml)."
    )

OUT = config["outdir"].rstrip("/")
IDX = config["index_dir"].rstrip("/")
IMG = config["containers"]
SCRIPTS = os.path.join(workflow.basedir, "scripts")
# Tool registry. Each tool is one module under rules/tools/<name>.smk.
#   callers    : BAMs in -> <OUT>/callers/<name>/... read by build_catalogue.py
#   annotators : scored ORFs in -> <OUT>/annotations/<name>.tsv, first column
#                orf_id, exactly one row per ORF in scoring/orf_candidates.bed12
#   counters   : BAMs + annotated ORFs in -> <OUT>/counts.tsv
#   visualizers: P-site tracks + ORFs in -> <OUT>/tracks/ (browser files)
#   conservation: ORF proteins + other species in -> <OUT>/annotations/<name>.tsv,
#                same contract as an annotator (join_annotations merges it)
#   assembly   : bulk RNA-seq in -> <OUT>/assembly/novel.gtf, used as novel_gtf
KNOWN_TOOLS = {
    "callers": ["price", "ribocode", "ribotish"],
    "annotators": ["annotator", "shortstop", "blastoff"],
    "counters": ["featurecounts"],
    "visualizers": ["frame_tracks"],
    "conservation": ["diamond"],
    "assembly": ["stringtie"],
}
# An absent kind means all of it, except visualizers, conservation and assembly: configs
# written before they existed should not start building tracks, searching species or assembling.
# Absent callers means PRICE + RiboCode: Ribo-TISH needs TIS libraries in the samplesheet.
OPT_IN = ("visualizers", "conservation", "assembly")
DEFAULTS = {**KNOWN_TOOLS, "callers": ["price", "ribocode"]}
TOOLS = {kind: list(config.get("tools", {}).get(kind, [] if kind in OPT_IN else DEFAULTS[kind]))
         for kind in KNOWN_TOOLS}
for _kind, _names in TOOLS.items():
    _unknown = set(_names) - set(KNOWN_TOOLS[_kind])
    if _unknown:
        raise WorkflowError(f"Unknown {_kind} {sorted(_unknown)}; available: {KNOWN_TOOLS[_kind]}.")
if PROVIDED:
    _set = [k for k in ("callers", "counters", "visualizers", "assembly") if config.get("tools", {}).get(k)]
    if _set:
        raise WorkflowError(f"orf_gtf is set (provided mode: annotation only), so tools.{_set[0]} "
                            "must be [] or absent; they need Ribo-seq.")
    TOOLS["callers"], TOOLS["counters"], TOOLS["visualizers"], TOOLS["assembly"] = [], [], [], []
elif not TOOLS["callers"]:
    raise WorkflowError("tools.callers is empty; at least one ORF caller is required.")
CALLERS = TOOLS["callers"]
ANNOTATORS = TOOLS["annotators"]
COUNTERS = TOOLS["counters"]
VISUALIZERS = TOOLS["visualizers"]
CONSERVATION = TOOLS["conservation"]
ASSEMBLY = TOOLS["assembly"]
if ASSEMBLY and config.get("novel_gtf"):
    raise WorkflowError("Set novel_gtf or tools.assembly, not both: the assembly is the novel GTF.")
# Novel transcript models for the hybrid GTFs (rules/reference.smk), or None.
NOVEL_GTF = f"{config['outdir'].rstrip('/')}/assembly/novel.gtf" if ASSEMBLY else config.get("novel_gtf")
# Antisense (x) models are only real if the novel GTF came from strand-specific RNA-seq.
# Assembly: rnaseq_strand must measure every library rf or fr. Supplied GTF: declared.
NEED_STRANDED = bool(NOVEL_GTF) and "x" in (set(config.get("novel_class_codes") or [])
                                            | set(config.get("ribocode_novel_class_codes") or []))
if NEED_STRANDED and not ASSEMBLY and config.get("novel_gtf_stranded") is not True:
    raise WorkflowError("novel_class_codes has 'x' (antisense), which needs a strand-specific novel_gtf. "
                        "Set novel_gtf_stranded: true if its RNA-seq was stranded, else drop 'x'.")
if NEED_STRANDED and ASSEMBLY and str(config.get("rnaseq_strandedness", "auto")) == "unstranded":
    raise WorkflowError("novel_class_codes has 'x' (antisense) but rnaseq_strandedness is unstranded; drop 'x'.")
if COUNTERS and "blastoff" not in ANNOTATORS:
    raise WorkflowError("tools.counters needs blastoff in tools.annotators: novel smORFs are "
                        "the ones blastoff finds no Swiss-Prot/TrEMBL match for.")


# Library types (optional samplesheet column `type`, default elongation):
#   elongation  ordinary Ribo-seq: P-sites, PRICE, RiboCode, scoring, counting, tracks
#   tis_harr / tis_ltm  initiation-enriched (harringtonine / lactimidomycin): trimmed,
#       aligned and strand-checked, then used only by Ribo-TISH, paired with the
#       elongation libraries of the same `condition` (optional column, default all).
RIBO_TYPES = ("elongation", "tis_harr", "tis_ltm")
# A name that suggests a TIS library, with `type` left blank, stops the run: an
# unlabelled TIS library would be pooled with elongation in every caller. An explicit
# `elongation` is accepted; psite_report's start-codon peak check still applies.
TIS_NAME_HINT = re.compile(r"harr|ltm|lactimidomycin|(?<![a-z])(lact|tis)(?![a-z])", re.IGNORECASE)


def _read_samples(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or not {"sample", "fastq"} <= set(rows[0]):
        raise WorkflowError(f"{path}: need tab-separated columns 'sample' and 'fastq'.")
    samples, types, conds, unlabelled = {}, {}, {}, []
    for r in rows:
        s, fq = r["sample"].strip(), r["fastq"].strip()
        if s in samples:
            raise WorkflowError(f"{path}: duplicate sample {s}.")
        if not os.path.isabs(fq):
            raise WorkflowError(f"{path}: fastq for {s} must be an absolute path.")
        t = (r.get("type") or "").strip()
        if not t and (TIS_NAME_HINT.search(s) or TIS_NAME_HINT.search(os.path.basename(fq))):
            unlabelled.append(s)
        t = t or "elongation"
        if t not in RIBO_TYPES:
            raise WorkflowError(f"{path}: type of {s} is '{t}'; use one of {list(RIBO_TYPES)}.")
        samples[s], types[s], conds[s] = fq, t, (r.get("condition") or "all").strip()
    if unlabelled:
        raise WorkflowError(
            f"{path}: {', '.join(unlabelled)} look like TIS libraries by name (harr / ltm / lact / tis) "
            "but have no type. Set type to tis_harr or tis_ltm (with a condition shared with "
            "elongation libraries), or to elongation if they are ordinary Ribo-seq.")
    return samples, types, conds


SAMPLE_FASTQ, SAMPLE_TYPE, SAMPLE_COND = ({}, {}, {}) if PROVIDED else _read_samples(config["samplesheet"])
ALL_SAMPLES = sorted(SAMPLE_FASTQ)          # trimming, alignment, adapter + strand checks
SAMPLES = [s for s in ALL_SAMPLES if SAMPLE_TYPE[s] == "elongation"]
TIS_SAMPLES = [s for s in ALL_SAMPLES if SAMPLE_TYPE[s] != "elongation"]
if not PROVIDED:
    if not SAMPLES:
        raise WorkflowError("The samplesheet has no elongation libraries; every caller but Ribo-TISH needs them.")
    if TIS_SAMPLES and "ribotish" not in CALLERS:
        raise WorkflowError(f"TIS libraries {TIS_SAMPLES} are only used by Ribo-TISH; add ribotish to "
                            "tools.callers or drop them from the samplesheet.")
    if "ribotish" in CALLERS and not TIS_SAMPLES:
        raise WorkflowError("tools.callers has ribotish but the samplesheet has no TIS libraries "
                            "(type tis_harr or tis_ltm).")
    _orphan = [t for t in TIS_SAMPLES if not any(SAMPLE_COND[s] == SAMPLE_COND[t] for s in SAMPLES)]
    if _orphan:
        raise WorkflowError(f"TIS libraries {_orphan} have no elongation library with the same condition.")

wildcard_constraints:
    sample="|".join(map(re.escape, ALL_SAMPLES)),
    strand="fw|rc",


def scaled(mb):
    """Memory request that doubles on the profile's one retry."""
    return lambda wc, attempt: mb * attempt
