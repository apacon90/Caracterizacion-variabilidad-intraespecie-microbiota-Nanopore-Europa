#!/usr/bin/env bash
# 13_indels_snps_VCF.sh — Variant calling with bcftools mpileup+call per country, merge by species
# Input:  OUT_BAM/{species}/{country}_pool{01..10}.bam  (from 12_mapp_call.sh)
#         REF_DIR/{species}.fa
# Output: VCF_RAW_DIR/{species}_Spain|Norway|Portugal.vcf.gz
#         VCF_RAW_DIR/{species}.vcf.gz  (merged)
# Usage:  bash 13_indels_snps_VCF.sh

set -euo pipefail
source "$(dirname "$0")/config.sh"

REF_DIR="${RESULTS_DIR}/ref_per_species"
BAM_BASE="${RESULTS_DIR}/pool_mappings"
OUT_VCF="${VCF_RAW_DIR}"

THREADS_PER_JOB=8
MAPQ_MIN=20
BQ_MIN=10
MAX_BQ=35
MAX_DEPTH=500
PLOIDY=1
P_PRIOR="1e-3"
ANNOT="FORMAT/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR"

mkdir -p "$OUT_VCF"
shopt -s nullglob

call_one_country() {
    local out_vcf="$1"; shift
    local -a bams=("$@")
    [[ ${#bams[@]} -gt 0 ]] || return 0

    bcftools mpileup -Ou -f "$REF" \
        -B \
        -q "$MAPQ_MIN" -Q "$BQ_MIN" \
        --max-BQ "$MAX_BQ" \
        --max-depth "$MAX_DEPTH" \
        -a "$ANNOT" \
        --threads "$THREADS_PER_JOB" \
        "${bams[@]}" \
    | bcftools call -mv \
        --ploidy "$PLOIDY" -P "$P_PRIOR" \
        --threads "$THREADS_PER_JOB" \
        -Oz -o "$out_vcf"
    tabix -f -p vcf "$out_vcf"
}

find "$BAM_BASE" -maxdepth 1 -mindepth 1 -type d | sort | while read -r SPEC_DIR; do
    SPECIES=$(basename "$SPEC_DIR")
    REF="$REF_DIR/${SPECIES}.fa"
    [[ -s "$REF" ]] || { echo "WARNING: no reference for $SPECIES, skipping"; continue; }

    BAMS_ES=("$SPEC_DIR"/Spain_pool*.bam)
    BAMS_NO=("$SPEC_DIR"/Norway_pool*.bam)
    BAMS_PT=("$SPEC_DIR"/Portugal_pool*.bam)

    VCF_ES="$OUT_VCF/${SPECIES}_Spain.vcf.gz"
    VCF_NO="$OUT_VCF/${SPECIES}_Norway.vcf.gz"
    VCF_PT="$OUT_VCF/${SPECIES}_Portugal.vcf.gz"

    # Call variants per country in parallel
    call_one_country "$VCF_ES" "${BAMS_ES[@]}" &
    call_one_country "$VCF_NO" "${BAMS_NO[@]}" &
    call_one_country "$VCF_PT" "${BAMS_PT[@]}" &
    wait

    # Merge country VCFs into single species VCF
    MERGE=()
    [[ -s "$VCF_ES" ]] && MERGE+=("$VCF_ES")
    [[ -s "$VCF_NO" ]] && MERGE+=("$VCF_NO")
    [[ -s "$VCF_PT" ]] && MERGE+=("$VCF_PT")

    OUT="$OUT_VCF/${SPECIES}.vcf.gz"
    if [[ ${#MERGE[@]} -eq 0 ]]; then
        echo "WARNING: no VCFs for $SPECIES"
    elif [[ ${#MERGE[@]} -eq 1 ]]; then
        cp -f "${MERGE[0]}" "$OUT"; tabix -f -p vcf "$OUT"
    else
        bcftools merge -Oz "${MERGE[@]}" -o "$OUT"
        tabix -f -p vcf "$OUT"
    fi
    echo "Done: $SPECIES"
done

echo "Done. VCFs saved in: $OUT_VCF"
