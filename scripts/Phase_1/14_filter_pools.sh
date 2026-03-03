#!/usr/bin/env bash
# 14_filter_pools.sh — Filter variants by quality (QUAL) and recurrence (>=3 pools per country)
# Input:  VCF_RAW_DIR/{species}_Spain|Norway|Portugal.vcf.gz  (from 13_indels_snps_VCF.sh)
# Output: VCF_FILTERED_DIR/{species}.FINAL_pool3.vcf.gz
#
# Filters applied:
#   SNPs:   QUAL >= 100
#   INDELs: QUAL >= 150
#   Recurrence: GT=1 in >= 3 of 10 pools per country
#
# Usage: bash 14_filter_pools.sh

set -euo pipefail
source "$(dirname "$0")/config.sh"

REF_DIR="${RESULTS_DIR}/ref_per_species"
QUAL_SNP=100
QUAL_INDEL=150
FILTER="((TYPE=\"snp\" && QUAL>=${QUAL_SNP}) || (TYPE=\"indel\" && QUAL>=${QUAL_INDEL})) && COUNT(GT=\"1\")>=3"

mkdir -p "$VCF_FILTERED_DIR"

filter_country() {
    local in_vcf="$1" ref="$2" out_vcf="$3"
    [[ -s "$in_vcf" ]] || return 0
    bcftools norm -f "$ref" -m -both "$in_vcf" -Ou \
    | bcftools view -i "$FILTER" -Oz -o "$out_vcf"
    tabix -f -p vcf "$out_vcf"
}

for VCF in "$VCF_RAW_DIR"/*_Spain.vcf.gz; do
    SPECIES=$(basename "$VCF" _Spain.vcf.gz)
    REF="$REF_DIR/${SPECIES}.fa"
    [[ -s "$REF" ]] || { echo "WARNING: no reference for $SPECIES, skipping"; continue; }

    echo "Filtering $SPECIES"

    TMP_ES="$VCF_FILTERED_DIR/${SPECIES}_Spain.pool3.vcf.gz"
    TMP_NO="$VCF_FILTERED_DIR/${SPECIES}_Norway.pool3.vcf.gz"
    TMP_PT="$VCF_FILTERED_DIR/${SPECIES}_Portugal.pool3.vcf.gz"

    filter_country "$VCF_RAW_DIR/${SPECIES}_Spain.vcf.gz"   "$REF" "$TMP_ES"
    filter_country "$VCF_RAW_DIR/${SPECIES}_Norway.vcf.gz"  "$REF" "$TMP_NO"
    filter_country "$VCF_RAW_DIR/${SPECIES}_Portugal.vcf.gz" "$REF" "$TMP_PT"

    MERGE=()
    [[ -s "$TMP_ES" ]] && MERGE+=("$TMP_ES")
    [[ -s "$TMP_NO" ]] && MERGE+=("$TMP_NO")
    [[ -s "$TMP_PT" ]] && MERGE+=("$TMP_PT")

    OUT="$VCF_FILTERED_DIR/${SPECIES}.FINAL_pool3.vcf.gz"

    if [[ ${#MERGE[@]} -eq 0 ]]; then
        echo "WARNING: no variants passed filters for $SPECIES"
    else
        bcftools merge -Oz -o "$OUT" "${MERGE[@]}"
        tabix -f -p vcf "$OUT"
        rm -f "$TMP_ES"* "$TMP_NO"* "$TMP_PT"*
        echo "Done: $SPECIES → $(basename "$OUT")"
    fi
done

echo "Done. Filtered VCFs saved in: $VCF_FILTERED_DIR"
