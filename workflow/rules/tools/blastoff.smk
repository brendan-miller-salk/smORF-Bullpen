# blastoff: BLASTP of each ORF's protein against the UniProt mouse reference
# proteome (canonical + isoforms). Labels ORFs as Swiss-Prot/TrEMBL identical,
# N-/C-truncated, internal, extended, insertion / deletion (an extra or skipped
# internal segment), partial, homolog, or novel. When the annotator is also enabled,
# join_annotations relabels its iORF and overlap calls to annotated_CDS / Trunc-Iso /
# N-, C-, NC-Iso / Ins-Iso / Del-Iso / TrEMBL (original kept in smorf_type_annotator).
# Needs config `blastoff_fasta` (written by refs/fetch_references.sh).
# Contract: <OUT>/annotations/blastoff.tsv

BO = f"{OUT}/annotations/blastoff"
BO_SCRIPTS = f"{SCRIPTS}/blastoff"
BLAST_FMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore"

if "blastoff_fasta" not in config:
    raise WorkflowError("tools.annotators has 'blastoff' but config has no blastoff_fasta.")
if "blast" not in IMG:
    raise WorkflowError("tools.annotators has 'blastoff' but containers has no 'blast' "
                        "(template: docker://quay.io/biocontainers/blast:2.14.1--pl5321h6f7f691_0).")


rule blastoff_db:
    input:
        config["blastoff_fasta"],
    output:
        multiext(f"{BO}/db/blastoff", ".pin", ".phr", ".psq"),
    params:
        prefix=f"{BO}/db/blastoff",
    container:
        IMG["blast"]
    log:
        f"{OUT}/logs/blastoff/db.log",
    shell:
        "makeblastdb -in {input} -dbtype prot -parse_seqids -out {params.prefix} > {log} 2>&1"


rule blastoff_query:
    input:
        bed=f"{OUT}/scoring/orf_candidates.bed12",
        genome=config["genome"],
    output:
        short=f"{BO}/query_short.fa",
        long=f"{BO}/query_long.fa",
    params:
        short_aa=config.get("blastoff_short_aa", 30),
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/blastoff/query.log",
    shell:
        r"""
        python {BO_SCRIPTS}/blastoff.py query --bed12 {input.bed} --genome {input.genome} \
            --short-aa {params.short_aa} --out-short {output.short} --out-long {output.long} 2> {log}
        cat {log} >&2
        """


rule blastoff_blast:
    """Short queries need blastp-short and a loose search E-value. A perfect 15-aa
    match has E ~ 0.4 against the proteome, so near-identity is judged from the
    alignment in classify, not from the E-value. Long queries use PANDA's settings,
    but keep up to 10 HSPs per protein: an inserted or skipped exon splits the
    alignment, and classify chains the pieces back into one hit."""
    input:
        short=f"{BO}/query_short.fa",
        long=f"{BO}/query_long.fa",
        db=rules.blastoff_db.output,
    output:
        short=f"{BO}/blast_short.tsv",
        long=f"{BO}/blast_long.tsv",
    params:
        db=f"{BO}/db/blastoff",
        fmt=BLAST_FMT,
    threads: 16
    resources:
        mem_mb=scaled(16000),
        # 3 h, doubled on the profile's retry. MEASURED: ~30 min for 77k ORFs on 16 cpu.
        runtime=lambda wc, attempt: 180 * attempt,
    container:
        IMG["blast"]
    log:
        f"{OUT}/logs/blastoff/blast.log",
    shell:
        r"""
        : > {log}
        : > {output.short}; : > {output.long}
        if [[ -s {input.short} ]]; then
            blastp -task blastp-short -query {input.short} -db {params.db} -evalue 10 \
                -outfmt '{params.fmt}' -max_target_seqs 50 -max_hsps 1 -seg no \
                -num_threads {threads} -out {output.short} >> {log} 2>&1
        fi
        if [[ -s {input.long} ]]; then
            blastp -query {input.long} -db {params.db} -evalue 1e-4 -word_size 6 \
                -soft_masking false -outfmt '{params.fmt}' -max_target_seqs 50 -max_hsps 10 \
                -seg no -num_threads {threads} -out {output.long} >> {log} 2>&1
        fi
        """


rule blastoff:
    input:
        blast=[f"{BO}/blast_short.tsv", f"{BO}/blast_long.tsv"],
        queries=[f"{BO}/query_short.fa", f"{BO}/query_long.fa"],
        fasta=config["blastoff_fasta"],
        bed=f"{OUT}/scoring/orf_candidates.bed12",
    output:
        f"{OUT}/annotations/blastoff.tsv",
    params:
        pident=config.get("blastoff_min_pident", 95),
        slack=config.get("blastoff_slack_aa", 2),
        shared=config.get("blastoff_min_shared_aa", 10),
    conda:
        "../../envs/score.yaml"
    log:
        f"{OUT}/logs/blastoff/classify.log",
    shell:
        r"""
        python {BO_SCRIPTS}/blastoff.py classify --blast {input.blast} \
            --uniprot-fasta {input.fasta} --bed12 {input.bed} --queries {input.queries} \
            --min-pident {params.pident} --slack {params.slack} --min-shared-aa {params.shared} --out {output} 2> {log}
        cat {log} >&2
        """
