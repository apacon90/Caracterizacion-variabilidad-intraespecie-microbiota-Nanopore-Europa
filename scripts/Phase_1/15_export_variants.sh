#!/usr/bin/env bash

# 15_export_variants.sh — Export filtered variants to TSV with genotypes per pool
# Input:  VCF_FILTERED_DIR/{species}.FINAL_pool3.vcf.gz  (from 14_filter_pools.sh)
# Output: RESULTS_DIR/ALL_species_pool3_variants_with_pools.tsv
#
# Usage: bash 15_export_variants.sh

set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUT="${TABLES_DIR}/ALL_species_pool3_variants_with_pools.tsv"

# Get pool order from first VCF
first_vcf=$(ls "$VCF_FILTERED_DIR"/*.FINAL_pool3.vcf.gz | head -n 1)
[[ -z "$first_vcf" ]] && { echo "ERROR: no VCFs found in $VCF_FILTERED_DIR" >&2; exit 1; }
mapfile -t SAMPLES < <(bcftools query -l "$first_vcf")

# Header
{
    echo -ne "Species\tCHROM\tPOS\tREF\tALT\tTYPE\tQUAL\tINFO\tFORMAT"
    for S in "${SAMPLES[@]}"; do echo -ne "\t${S}"; done
    echo
} > "$OUT"

# One row per variant per species
for VCF in "$VCF_FILTERED_DIR"/*.FINAL_pool3.vcf.gz; do
    [[ -s "$VCF" ]] || continue
    SPECIES=$(basename "$VCF" .FINAL_pool3.vcf.gz)
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%QUAL\t%INFO\t%FORMAT[\t%GT]\n' \
        "$VCF" \
    | awk -v sp="$SPECIES" 'BEGIN{OFS="\t"} {print sp,$0}' >> "$OUT"
done

echo "Done. Variants table saved: $OUT"
