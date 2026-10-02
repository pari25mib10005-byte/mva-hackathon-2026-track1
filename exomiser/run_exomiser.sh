#!/usr/bin/env bash
# Gene-agnostic Exomiser run for Track 1 (report §4.1).
# Requires: Java 17+, Exomiser CLI 15.1.0, data release 2512 (hg38 + phenotype)
# unzipped into $EXOMISER_DIR/data. ~9 GB RAM. Runtime ~1 min.
set -euo pipefail
EXOMISER_DIR=${EXOMISER_DIR:-$HOME/exomiser/exomiser-cli-15.1.0}
VCF=${1:?usage: run_exomiser.sh <proband.vcf.gz> <proband_phenopacket.yml>}
SAMPLE=${2:?usage: run_exomiser.sh <proband.vcf.gz> <proband_phenopacket.yml>}
OUT=${OUT:-results}

cd "$EXOMISER_DIR"
java -Xmx9g -jar exomiser-cli-15.1.0.jar analyse \
  --sample="$SAMPLE" \
  --vcf="$VCF" \
  --assembly=GRCh38 \
  --preset=exome \
  --output-directory="$OUT" \
  --output-filename=PROBAND01 \
  --output-format=HTML,JSON,TSV_GENE,TSV_VARIANT

cut -f1-6 "$OUT/PROBAND01.genes.tsv" | head -11
