#!/usr/bin/env bash
# Fetch mouse GRCm39 references: Ensembl release 110, renamed to chr-convention.
#
# Why Ensembl 110 rather than GENCODE vM33: they are the same annotation on the
# same assembly, but only the Ensembl GTF carries the two things the Annotator
# needs -- `five_prime_utr` / `three_prime_utr` features (GENCODE has a single
# undifferentiated `UTR`) and `gene_biotype` / `transcript_biotype` attributes
# (GENCODE spells them `gene_type` / `transcript_type`). Fed a GENCODE GTF the
# Annotator does not fail, it silently emits zero uORFs, dORFs, lncRNA and
# psORF and dumps everything into `UA`.
#
# Ensembl 110 == GENCODE vM33. A novel GTF (e.g. a StringTie assembly) should be
# reference-guided with the same release so provenance stays consistent.
#
# Why rename to chr-convention: Ensembl names chromosomes `1`, `X`, `MT`; the
# workflow uses UCSC-style `chr1`, `chrX`, `chrM`. The Annotator additionally
# hard-codes `chr` in two places
# (`line.startswith('chr')` and `COORD_RE`). Scaffolds (GL*/JH*/MU*) already
# share INSDC accessions between the two, so only main chromosomes change.
#
# Usage: bash refs/fetch_references.sh <out_dir> [novel.gtf]
# Also fetches the UniProt mouse proteome for the blastoff annotator (FETCH_UNIPROT=0 skips it).
# ~3 GB into <out_dir>; point `genome`/`gtf` in the run config at it. Give the
# optional novel GTF (your `novel_gtf`) to check its seqnames match the genome.
set -euo pipefail

REF_DIR=${1:?usage: fetch_references.sh <out_dir> [novel.gtf]}
STRINGTIE_GTF=${2:-}
[[ -z "$STRINGTIE_GTF" || -f "$STRINGTIE_GTF" ]] || { echo "not found: $STRINGTIE_GTF" >&2; exit 1; }
[[ -z "$STRINGTIE_GTF" ]] || STRINGTIE_GTF=$(cd "$(dirname "$STRINGTIE_GTF")" && pwd)/$(basename "$STRINGTIE_GTF")
mkdir -p "$REF_DIR"
REF_DIR=$(cd "$REF_DIR" && pwd)
ENSEMBL="https://ftp.ensembl.org/pub/release-110"

cd "$REF_DIR"

# 1 .. 19, X, Y -> chrN ; MT -> chrM ; scaffolds untouched.
rename_gtf() {
  awk 'BEGIN{FS=OFS="\t"}
       /^#/ {print; next}
       { if ($1 ~ /^([0-9]+|X|Y)$/) $1="chr" $1; else if ($1=="MT") $1="chrM"; print }'
}

rename_fasta() {
  awk '/^>/ { n=substr($1,2)
              if (n ~ /^([0-9]+|X|Y)$/) n="chr" n; else if (n=="MT") n="chrM"
              print ">" n; next }
       { print }'
}

# --- Genome ------------------------------------------------------------------
if [[ ! -f GRCm39.primary_assembly.genome.fa ]]; then
  echo "Downloading genome ..."
  curl -fL --retry 3 -o ens_genome.fa.gz \
    "$ENSEMBL/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
  echo "Renaming seqnames to chr-convention ..."
  gunzip -c ens_genome.fa.gz | rename_fasta > GRCm39.primary_assembly.genome.fa
  rm -f ens_genome.fa.gz
fi

# --- Annotation --------------------------------------------------------------
if [[ ! -f ensembl110.GRCm39.chr.gtf ]]; then
  echo "Downloading annotation ..."
  curl -fL --retry 3 -o ens.gtf.gz \
    "$ENSEMBL/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz"
  gunzip -c ens.gtf.gz | rename_gtf > ensembl110.GRCm39.chr.gtf
  rm -f ens.gtf.gz
fi

# --- Transcript FASTA (only needed if you add matched RNA-seq) ---------------
if [[ ! -f ensembl110.GRCm39.cdna.fa ]]; then
  curl -fL --retry 3 -o ens_cdna.fa.gz \
    "$ENSEMBL/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
  gunzip -c ens_cdna.fa.gz > ensembl110.GRCm39.cdna.fa
  rm -f ens_cdna.fa.gz
fi

# --- UniProt mouse reference proteome (skip with FETCH_UNIPROT=0) -------------
# For the blastoff annotator (workflow/scripts/blastoff/blastoff.py). Canonical + _additional (Swiss-Prot
# and TrEMBL isoforms -- truncations often match an isoform, not the canonical).
# Pinned: UniProt only serves the current release here, so a different release
# aborts rather than silently changing the database under an existing run.
UNIPROT_RELEASE="${UNIPROT_RELEASE:-2026_03}"
UNIPROT_DIR="uniprot_UP000000589_${UNIPROT_RELEASE}"
if [[ "${FETCH_UNIPROT:-1}" == 1 && ! -f "$UNIPROT_DIR/uniprot_mouse.fasta" ]]; then
  UP="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589"
  mkdir -p "$UNIPROT_DIR"
  curl -fsL --retry 3 -o "$UNIPROT_DIR/RELEASE.metalink" "$UP/RELEASE.metalink"
  served=$(sed -n 's:.*<version>\(.*\)</version>.*:\1:p' "$UNIPROT_DIR/RELEASE.metalink")
  [[ "$served" == "$UNIPROT_RELEASE" ]] || {
    echo "FAIL: UniProt now serves $served, pinned $UNIPROT_RELEASE. Set UNIPROT_RELEASE=$served to move."; exit 1; }
  for f in UP000000589_10090.fasta.gz UP000000589_10090_additional.fasta.gz; do
    echo "Downloading UniProt $f ($UNIPROT_RELEASE) ..."
    curl -fL --retry 3 -o "$UNIPROT_DIR/$f" "$UP/$f"
    want=$(awk -v f="$f" '$0 ~ "name=\""f"\"" {g=1} g && /hash type="md5"/ {gsub(/.*">|<.*/,""); print; exit}' "$UNIPROT_DIR/RELEASE.metalink")
    got=$(md5sum "$UNIPROT_DIR/$f" | cut -d' ' -f1)
    [[ -n "$want" && "$want" == "$got" ]] || { echo "FAIL: md5 mismatch for $f ($got vs '$want')"; exit 1; }
  done
  gunzip -c "$UNIPROT_DIR"/UP000000589_10090.fasta.gz "$UNIPROT_DIR"/UP000000589_10090_additional.fasta.gz \
    > "$UNIPROT_DIR/uniprot_mouse.fasta"
  # The workflow builds its own BLAST db from the FASTA; this one is for standalone use.
  if command -v makeblastdb >/dev/null; then
    makeblastdb -in "$UNIPROT_DIR/uniprot_mouse.fasta" -dbtype prot -parse_seqids \
      -out "$UNIPROT_DIR/uniprot_mouse" -title "UniProt UP000000589 $UNIPROT_RELEASE canonical+additional"
  else
    echo "  makeblastdb not on PATH: skipped the standalone BLAST db (the workflow builds its own)"
  fi
  echo "  UniProt entries: $(grep -c '^>' "$UNIPROT_DIR/uniprot_mouse.fasta")"
fi

# --- Verification ------------------------------------------------------------
echo
echo "=== Annotator requirements ==="
for feat in five_prime_utr three_prime_utr; do
  n=$(awk -F'\t' -v f="$feat" '$3==f' ensembl110.GRCm39.chr.gtf | wc -l | tr -d ' ')
  printf "  %-20s %s\n" "$feat" "$n"
  [[ "$n" -gt 0 ]] || { echo "  FAIL: $feat missing - the Annotator cannot classify uORFs/dORFs"; exit 1; }
done
for attr in gene_biotype transcript_biotype; do
  grep -q "$attr \"" ensembl110.GRCm39.chr.gtf \
    && printf "  %-20s present\n" "$attr" \
    || { echo "  FAIL: $attr missing"; exit 1; }
done

echo
echo "=== Seqname agreement with the genome ==="
awk '/^>/{print substr($1,2)}' GRCm39.primary_assembly.genome.fa | sort -u > .fa_names
cut -f1 ensembl110.GRCm39.chr.gtf | grep -v '^#' | sort -u > .gtf_names
missing=$(comm -13 .fa_names .gtf_names | wc -l | tr -d ' ')
echo "  GTF seqnames absent from genome FASTA: $missing"
[[ "$missing" -eq 0 ]] || { echo "  FAIL:"; comm -13 .fa_names .gtf_names | head; exit 1; }

if [[ -n "$STRINGTIE_GTF" ]]; then
  echo
  echo "=== Seqname agreement with the StringTie assembly ==="
  cut -f1 "$STRINGTIE_GTF" | sort -u > .st_names
  bad=$(comm -13 .fa_names .st_names | wc -l | tr -d ' ')
  echo "  StringTie seqnames absent from genome FASTA: $bad"
  [[ "$bad" -eq 0 ]] || { echo "  FAIL - the workflow will error on these:"; comm -13 .fa_names .st_names | head; exit 1; }
fi

rm -f .fa_names .gtf_names .st_names
echo
echo "All checks passed. References in $REF_DIR:"
ls -lh "$REF_DIR" | grep -vE '^total'
