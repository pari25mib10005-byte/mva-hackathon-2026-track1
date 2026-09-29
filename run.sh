#!/usr/bin/env bash
# Reproduces the Track 1 analysis end to end. Requires: conda env with bcftools + python requests.
# Usage: bash run.sh /path/to/WGS_EX2312012_HGWCNDSX7.vcf.gz
set -euo pipefail
VCF="$1"

# 1. Build check (GRCh38: contig 15 length 101991189)
bcftools view -h "$VCF" | grep -w "contig=<ID=15"

# 2. Extract the three established MVA genes (+/-10 kb), quality filter, split multi-allelics
bcftools view -R mva_genes.bed "$VCF" -Oz -o candidates.vcf.gz
bcftools view -i 'QUAL>=30 && FMT/DP>=10 && FMT/GQ>=20' candidates.vcf.gz -Oz -o candidates.pass.vcf.gz
bcftools norm -m -any candidates.pass.vcf.gz -Oz -o candidates.norm.vcf.gz
bcftools index -f candidates.norm.vcf.gz

# 3. Annotate via Ensembl VEP REST and rank by impact then rarity
python annotate.py
