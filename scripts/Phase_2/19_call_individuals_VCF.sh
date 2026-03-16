#!/usr/bin/env bash

# 19_call_individuals_VCF.sh — Joint variant calling across all individuals per species
# Input:  IND_BAM_DIR/{species}/{country}/{sample}.bam  (from 18_map_individuals)
#         REF_DIR/{species}.fa  (from Phase 1)
# Output: IND_VCF_DIR/{species}.vcf.gz + .tbi
#
# Joint call: all individuals of a species (all countries) in a single bcftools call.
# Same parameters as Phase 1 pool calling.
# Usage: bash 19_call_individuals_VCF.sh

set -euo pipefail
shopt -s nullglob
source "$(dirname "$0")/../config.sh"


mkdir -p "$IND_VCF_DIR"

wait_for_slot() {
    while [[ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]]; do
        sleep 1
    done
}

for SPEC_DIR in "$IND_BAM_DIR"/*/; do
    SPECIES=$(basename "$SPEC_DIR")
    REF="$REF_DIR/${SPECIES}.fa"

    # Collect all indexed BAMs across countries
    BAMS=()
    for BAM in "$SPEC_DIR"/*/*.bam; do
        [[ -s "$BAM" && -s "${BAM}.bai" ]] && BAMS+=("$BAM")
    done
    [[ ${#BAMS[@]} -gt 0 ]] || continue

    OUT="$IND_VCF_DIR/${SPECIES}.vcf.gz"
    if [[ -s "$OUT" && -s "${OUT}.tbi" ]]; then
        echo "Skipping $SPECIES (VCF already exists)"
        continue
    fi

    echo "Joint call: $SPECIES (n_individuals=${#BAMS[@]})"
    wait_for_slot
    (
        bcftools mpileup -Ou -f "$REF" \
            -B \
            -q "$MAPQ_MIN" -Q "$BQ_MIN" \
            --max-BQ "$MAX_BQ" \
            --max-depth "$MAX_DEPTH" \
            -a "$ANNOT" \
            --threads "$THREADS_PER_JOB" \
            "${BAMS[@]}" \
        | bcftools call -mv \
            --ploidy "$PLOIDY" -P "$P_PRIOR" \
            --threads "$THREADS_PER_JOB" \
            -Oz -o "$OUT"

        tabix -f -p vcf "$OUT"
    ) &
done

wait
echo "Done. VCFs saved in: $IND_VCF_DIR"
