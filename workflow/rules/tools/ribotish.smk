# Caller: Ribo-TISH (Zhang et al. 2017), for initiation-enriched (TIS) libraries:
# harringtonine (type tis_harr, --harr) or lactimidomycin (tis_ltm). One predict per TIS
# library, with the elongation libraries of its condition as ordinary Ribo-seq (-b), so
# each call carries a TIS test (start-codon peak) and a frame test (elongation body).
# P-site offsets come from `ribotish quality` (-t for TIS) on each reference genome BAM.
# Candidates: every isoform + novel transcript (FULL_HYBRID, as PRICE); AUG plus
# near-cognate starts with ribotish_alt (--alt). TIS background and the known-CDS labels
# (-a): the canonical GTF, one transcript per gene; on the full human hybrid the
# background estimate alone ran > 11 min on chr22 (user, 2026-10-09). Calls are Ribo-TISH's own output, Fisher
# q (TIS + frame) <= ribotish_max_q (its --fsqth, default 0.05).
# Contract: <OUT>/callers/ribotish/<tis_sample>.txt, read by build_catalogue.py.

RT = f"{OUT}/callers/ribotish"
if "ribotish" not in IMG:
    raise WorkflowError("tools.callers has ribotish but containers has no 'ribotish' "
                        "(template: docker://quay.io/biocontainers/ribotish:0.2.8--pyhdfd78af_0).")


def ribotish_partners(tis):
    """The elongation libraries a TIS library is paired with: same condition."""
    return [s for s in SAMPLES if SAMPLE_COND[s] == SAMPLE_COND[tis]]


wildcard_constraints:
    tis="|".join(map(re.escape, TIS_SAMPLES)),


rule ribotish_quality:
    """Read-length offsets (.para.py). Canonical GTF: offsets come from annotated
    start/stop codons, one transcript per gene is enough and much faster. Fails if no
    read length passes Ribo-TISH's quality threshold, rather than predicting without."""
    input:
        bam=f"{OUT}/align/{{sample}}.genome.sorted.bam",
        bai=f"{OUT}/align/{{sample}}.genome.sorted.bam.bai",
        gtf=f"{IDX}/reference.canonical.gtf",
        strand=f"{OUT}/qc/strand/{{sample}}.tsv",   # gate: the library is sense-stranded
    output:
        para=f"{RT}/quality/{{sample}}.para.py",
        txt=f"{RT}/quality/{{sample}}_qual.txt",
        pdf=f"{RT}/quality/{{sample}}_qual.pdf",
    params:
        tis=lambda wc: "" if SAMPLE_TYPE[wc.sample] == "elongation" else "-t",
    threads: 4
    resources:
        mem_mb=scaled(16000),
        runtime=240,
    container:
        IMG["ribotish"]
    log:
        f"{OUT}/logs/ribotish/quality_{{sample}}.log",
    shell:
        r"""
        ribotish quality {params.tis} -b {input.bam} -g {input.gtf} -o {output.txt} \
            -f {output.pdf} -r {output.para} -p {threads} > {log} 2>&1
        grep -qE "offdict = \{{[0-9]" {output.para} || {{
            echo "ribotish quality: no read length passed for {wildcards.sample} ({output.txt})" >&2; exit 1; }}
        """


rule ribotish_predict:
    input:
        tis=f"{OUT}/align/{{tis}}.genome.sorted.bam",
        tis_bai=f"{OUT}/align/{{tis}}.genome.sorted.bam.bai",
        tispara=f"{RT}/quality/{{tis}}.para.py",
        ribo=lambda wc: [f"{OUT}/align/{s}.genome.sorted.bam" for s in ribotish_partners(wc.tis)],
        ribo_bai=lambda wc: [f"{OUT}/align/{s}.genome.sorted.bam.bai" for s in ribotish_partners(wc.tis)],
        ribopara=lambda wc: [f"{RT}/quality/{s}.para.py" for s in ribotish_partners(wc.tis)],
        gtf=FULL_HYBRID,
        canonical=f"{IDX}/reference.canonical.gtf",
        genome=config["genome"],
    output:
        f"{RT}/{{tis}}.txt",
    params:
        harr=lambda wc: "--harr" if SAMPLE_TYPE[wc.tis] == "tis_harr" else "",
        alt="--alt" if config.get("ribotish_alt", True) else "",
        q=config.get("ribotish_max_q", 0.05),
        ribo=lambda wc, input: ",".join(input.ribo),
        ribopara=lambda wc, input: ",".join(input.ribopara),
        bg=f"{RT}/{{tis}}.tisBackground.txt",
    # MEASURED (chr22, 1 process, --alt): background 45 s on the canonical GTF, prediction
    # > 20 min at ~0.5 GB, i.e. ~15-20 CPU-h genome-wide per TIS library.
    threads: 16
    resources:
        mem_mb=scaled(64000),
        runtime=1440,
    container:
        IMG["ribotish"]
    log:
        f"{OUT}/logs/ribotish/predict_{{tis}}.log",
    shell:
        r"""
        # Ribo-TISH's TIS-background pool (NoDaemonProcess, -p > 1) predates Python 3.8,
        # where Pool passes its context as the first Process argument; it lands in `group`
        # and fails an assertion (0.2.7 and 0.2.8 alike). The shim drops that context
        # and then runs the ribotish command unchanged.
        python - predict -t {input.tis} -b {params.ribo} -g {input.gtf} -a {input.canonical} -f {input.genome} \
            --tispara {input.tispara} --ribopara {params.ribopara} {params.harr} {params.alt} \
            --fsqth {params.q} --blocks --aaseq --allresult off -e {params.bg} \
            -p {threads} -o {output} > {log} 2>&1 <<'SHIM'
import multiprocessing.context as ctx, multiprocessing.process as proc, runpy, shutil, sys
_init = proc.BaseProcess.__init__
def _init_fixed(self, group=None, *a, **k):
    _init(self, None if isinstance(group, ctx.BaseContext) else group, *a, **k)
proc.BaseProcess.__init__ = _init_fixed
sys.argv[0] = shutil.which("ribotish")
runpy.run_path(sys.argv[0], run_name="__main__")
SHIM
        """
