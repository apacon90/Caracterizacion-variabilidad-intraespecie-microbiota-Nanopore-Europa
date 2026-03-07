#!/usr/bin/env bash
# 05_map_individuals.sh — Map individual reads against species-specific references
# Input:  FASTA_N/{species}/{country}/{species}.{sample}.n100.fa  (from 04_extract_reads)
#         REF_DIR/{species}.fa + .mmi  (from Phase 1)
# Output: IND_BAM_DIR/{species}/{country}/{sample}.bam + .bai
#
# Runs up to MAX_JOBS mappings in parallel.
# Usage: bash 05_map_individuals.sh

set -euo pipefail
shopt -s nullglob
source "$(dirname "$0")/config.sh"

FASTA_N="${PHASE2_DIR}/fasta_n100_or_max"
REF_DIR="${RESULTS_DIR}/ref_per_species"
IND_BAM_DIR="${PHASE2_DIR}/results_refdb/individual_mappings"
THREADS_PER_MAP=8
MAX_JOBS=4

mkdir -p "$IND_BAM_DIR"

wait_for_slot() {
    while [[ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]]; do
        sleep 1
    done
}

# Build minimap2 index per reference if missing
for REF in "$REF_DIR"/*.fa; do
    MMI="${REF%.fa}.mmi"
    [[ -f "$MMI" ]] && continue
    minimap2 -d "$MMI" "$REF"
    echo "Index created: $(basename "$MMI")"
done

# Map each individual
for SPEC_DIR in "$FASTA_N"/*/; do
    SPECIES=$(basename "$SPEC_DIR")
    REF="$REF_DIR/${SPECIES}.fa"
    MMI="${REF%.fa}.mmi"

    for COUNTRY_DIR in "$SPEC_DIR"/*/; do
        COUNTRY=$(basename "$COUNTRY_DIR")
        BAM_DIR="$IND_BAM_DIR/$SPECIES/$COUNTRY"
        mkdir -p "$BAM_DIR"

        for FA in "$COUNTRY_DIR"/*.n100.fa "$COUNTRY_DIR"/*.fa; do
            [[ -s "$FA" ]] || continue

            BASE=$(basename "$FA")
            NOEXT="${BASE%.fa}"; NOEXT="${NOEXT%.n100}"
            SAMPLE="${NOEXT#${SPECIES}.}"
            [[ "$SAMPLE" != "$NOEXT" ]] || SAMPLE="$NOEXT"

            BAM="$BAM_DIR/${SAMPLE}.bam"
            if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
                echo "Skipping $SPECIES / $COUNTRY / $SAMPLE (already done)"
                continue
            fi

            echo "Mapping $SPECIES / $COUNTRY / $SAMPLE"
            wait_for_slot
            (
                minimap2 -a -t "$THREADS_PER_MAP" -x map-ont --secondary=no \
                    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT\tLB:${SPECIES}\tPU:${COUNTRY}" \
                    "$MMI" "$FA" \
                | samtools sort -@ "$THREADS_PER_MAP" -o "$BAM" -
                samtools index "$BAM"
            ) &
        done
    done
done

wait
echo "Done. BAMs saved in: $IND_BAM_DIR"
