#!/usr/bin/env bash

# 23_vcf_to_tsv.sh — Export filtered VCFs to TSV format
# Input:  IND_VCF_FILT_DIR/{species}.QUAL_100-150.vcf.gz  (from 22_filter_qual_individuals)
# Output: TSV_DIR/{species}.GT.tsv    — genotypes only (CHROM, POS, REF, ALT, GT per individual)
#         TSV_DIR/{species}.FULL.tsv  — full info (GT, PL, DP, ADF, ADR, AD per individual)
#
# Usage: bash 21_vcf_to_tsv.sh

set -euo pipefail
shopt -s nullglob
source "$(dirname "$0")/../config.sh"

mkdir -p "$TSV_DIR"

for VCF in "$IND_VCF_FILT_DIR"/*.vcf.gz; do
    SPECIES=$(basename "$VCF" .vcf.gz)
    SAMPLES=$(bcftools query -l "$VCF" | paste -sd $'\t' -)

    # Genotypes only
    {
        echo -e "CHROM\tPOS\tREF\tALT\t${SAMPLES}"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$VCF"
    } > "$TSV_DIR/${SPECIES}.GT.tsv"

    # Full info per individual
    {
        echo -e "CHROM\tPOS\tREF\tALT\tFORMAT\t${SAMPLES}"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT[\t%TGT:%PL:%DP:%ADF:%ADR:%AD]\n' "$VCF"
    } > "$TSV_DIR/${SPECIES}.FULL.tsv"

    echo "OK: $SPECIES -> ${SPECIES}.GT.tsv | ${SPECIES}.FULL.tsv"
done

echo "Done. TSVs saved in: $TSV_DIR"
