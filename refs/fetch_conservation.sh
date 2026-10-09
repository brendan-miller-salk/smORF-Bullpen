#!/usr/bin/env bash
# Fetch other species' transcripts and orthologs for the diamond conservation tool,
# all one Ensembl release (the run annotation's release; default 110, mouse):
#   <name>.cdna_ncrna.fa.gz   cDNA + ncRNA, checked against Ensembl's CHECKSUMS (BSD `sum`)
#   <name>.orthologs.tsv      source gene -> <name> gene, orthology type (Compara, via a
#                             BioMart of that release; checked for BioMart's [success]
#                             and that the mart serves RELEASE)
#
# Usage: bash refs/fetch_conservation.sh <out_dir>
# ~1 GB into <out_dir>/ensembl<RELEASE>_conservation/. Prints the diamond_species and
# diamond_orthologs blocks for the run config.
#   SOURCE="name:biomart_prefix"  the run's species (default mouse:mmusculus)
#   RELEASE=<n>                   Ensembl release (default 110)
#   MART=<martservice url>        a BioMart serving RELEASE (default: the 110 archive)
#   SPECIES="name:ensembl_dir:biomart_prefix ..." overrides the panel (default: human,
#                                 mouse, rat, cow, chicken, zebrafish minus the source)
# Human run on Ensembl 116, e.g.:
#   SOURCE=human:hsapiens RELEASE=116 MART=https://www.ensembl.org/biomart/martservice \
#     bash refs/fetch_conservation.sh <out_dir>
set -euo pipefail

REF_DIR=${1:?usage: fetch_conservation.sh <out_dir>}
RELEASE=${RELEASE:-110}
IFS=: read -r SRC SRC_MART <<< "${SOURCE:-mouse:mmusculus}"
MART=${MART:-https://jul2023.archive.ensembl.org/biomart/martservice}   # = Ensembl 110
mkdir -p "$REF_DIR/ensembl${RELEASE}_conservation"
DIR=$(cd "$REF_DIR/ensembl${RELEASE}_conservation" && pwd)
ENSEMBL="https://ftp.ensembl.org/pub/release-${RELEASE}/fasta"
PANEL="human:homo_sapiens:hsapiens mouse:mus_musculus:mmusculus rat:rattus_norvegicus:rnorvegicus cow:bos_taurus:btaurus chicken:gallus_gallus:ggallus zebrafish:danio_rerio:drerio"
SPECIES=${SPECIES:-$(for e in $PANEL; do [[ ${e%%:*} == "$SRC" ]] || printf '%s ' "$e"; done)}

served=$(curl -fsSL --retry 3 "$MART?type=registry" | grep -o 'displayName="Ensembl Genes [0-9]*"' | grep -o '[0-9]*' | head -1)
[[ "$served" == "$RELEASE" ]] || { echo "FAIL: $MART serves Ensembl ${served:-?}, not $RELEASE; set MART"; exit 1; }

cd "$DIR"
for entry in $SPECIES; do
  IFS=: read -r name sp mart <<< "$entry"
  orth="$name.orthologs.tsv"
  if [[ ! -s "$orth" ]]; then
    echo "Fetching $SRC -> $name orthologs ..."
    q='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" completionStamp="1"><Dataset name="'$SRC_MART'_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="'$mart'_homolog_ensembl_gene"/><Attribute name="'$mart'_homolog_orthology_type"/></Dataset></Query>'
    curl -fsSL --retry 3 --max-time 1800 -G "$MART" --data-urlencode "query=$q" -o "$orth.raw"
    [[ "$(tail -1 "$orth.raw")" == "[success]" ]] || { echo "FAIL: BioMart response for $name is incomplete"; exit 1; }
    { printf '%s_gene_id\t%s_gene_id\torthology_type\n' "$SRC" "$name"
      awk -F'\t' '$2 != ""' "$orth.raw" | grep -v '^\[success\]$'; } > "$orth"
    rm -f "$orth.raw"
    echo "  $name: $(($(wc -l < "$orth") - 1)) ortholog pairs, $(tail -n +2 "$orth" | cut -f1 | sort -u | wc -l) $SRC genes"
  fi
  out="$name.cdna_ncrna.fa.gz"
  [[ -s "$out" ]] && { echo "have $out"; continue; }
  parts=()
  for kind in cdna ncrna; do
    url="$ENSEMBL/$sp/$kind"
    f=$(curl -fsL "$url/" | grep -o 'href="[^"]*\.\(cdna\.all\|ncrna\)\.fa\.gz"' | sed 's/href="//; s/"$//' | head -1)
    [[ -n "$f" ]] || { echo "FAIL: no $kind FASTA under $url"; exit 1; }
    echo "Downloading $f ..."
    curl -fL --retry 3 -o "$f" "$url/$f"
    want=$(curl -fsL "$url/CHECKSUMS" | awk -v f="$f" '$3 == f {print $1, $2}')
    got=$(sum "$f" | awk '{print $1, $2}')
    [[ -n "$want" && "$want" == "$got" ]] || { echo "FAIL: checksum mismatch for $f ($got vs '$want')"; exit 1; }
    parts+=("$f")
  done
  cat "${parts[@]}" > "$out.tmp" && mv "$out.tmp" "$out"   # concatenated gzip is valid gzip
  rm -f "${parts[@]}"
  echo "  $name: $(gzip -cd "$out" | grep -c '^>') transcripts"
done

echo
echo "diamond_species:"
for entry in $SPECIES; do
  name=${entry%%:*}
  echo "  $name: $DIR/$name.cdna_ncrna.fa.gz"
done
echo "diamond_orthologs:"
for entry in $SPECIES; do
  name=${entry%%:*}
  echo "  $name: $DIR/$name.orthologs.tsv"
done
