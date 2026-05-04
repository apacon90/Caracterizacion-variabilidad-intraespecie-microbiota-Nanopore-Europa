#!/usr/bin/env bash

# 12_mapp_call.sh — Extract species references, index and map pools with minimap2
# Input:  RESULTS_DIR/select_sp.tsv, POOLS_DIR/{species}/{country}_pool{01..10}.fa
# Output: REF_DIR/{species}.fa, OUT_BAM/{species}/{country}_pool{XX}.bam + .bai
# Usage:  bash 13_mapp_call.sh

set -euo pipefail
source "$(dirname "$0")/../config.sh"


mkdir -p "$REF_DIR" "$POOLS_DIR"

# Step 1: Extract reference FASTA per species from RefSeq207nr
tail -n +2 "$SELECT" | while IFS=$'\t' read -r QUERY CODE _ _ _ _ _; do
    REF="$REF_DIR/${QUERY}.fa"
    [[ -s "$REF" ]] && continue
    seqkit grep -r -p "^${CODE}$" "$DB_REF_FASTA" > "$REF"
    echo "$QUERY -> $CODE"
done
echo "References saved in: $REF_DIR"

# Step 2: Build minimap2 index per reference if missing
for REF in "$REF_DIR"/*.fa; do
    [[ -e "$REF" ]] || continue
    MMI="${REF%.fa}.mmi"
    [[ -f "$MMI" ]] && continue
    minimap2 -d "$MMI" "$REF"
    echo "Index created: $(basename "$MMI")"
done

# Step 3: Map each pool against its species-specific reference
for SPEC_DIR in "$POOLS_DIR"/*/; do
    [[ -d "$SPEC_DIR" ]] || continue
    SPECIES=$(basename "$SPEC_DIR")
    REF="$REF_DIR/${SPECIES}.fa"
    [[ -s "$REF" ]] || { echo "WARNING: no reference for $SPECIES, skipping"; continue; }
    MMI="${REF%.fa}.mmi"
    BAM_DIR="$POOLS_DIR/$SPECIES"
    mkdir -p "$BAM_DIR"

    for POOL in "$SPEC_DIR"/*.fa; do
        BASENAME=$(basename "$POOL" .fa)
        COUNTRY="${BASENAME%_pool*}"
        POOL_NO="${BASENAME##*_pool}"
        SAMPLE="${COUNTRY}_pool${POOL_NO}"
        BAM="$BAM_DIR/${SAMPLE}.bam"

        if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
            echo "Skipping $SPECIES / $SAMPLE (already done)"
            continue
        fi

        echo "Mapping $SPECIES / $SAMPLE"
        minimap2 -a -t "$THREADS" -x map-ont --secondary=no \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT\tLB:${SPECIES}\tPU:${COUNTRY}" \
            "$MMI" "$POOL" \
        | samtools sort -@ "$THREADS" -o "$BAM" -
        samtools index "$BAM"
    done
done

echo "Done. BAM files saved in: $POOLS_DIR"
