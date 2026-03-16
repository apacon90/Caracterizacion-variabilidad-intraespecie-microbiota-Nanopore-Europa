#!/bin/bash
# 01_trimm.sh — Trim first N nucleotides from 5' end with Cutadapt
# Input:  BASE_DIR/{country}_dataset/*.F.fasta.gz + *.R.fasta.gz
# Output: OUT_TRIMMED/{country}/{sample}.trim{N}.fa

set -euo pipefail
source "$(dirname "$0")/../config.sh"

mkdir -p "$OUT_TRIMMED"

for COUNTRY_DIR in "$BASE_DIR"/*_dataset; do
    COUNTRY=$(basename "$COUNTRY_DIR")
    mkdir -p "$OUT_TRIMMED/$COUNTRY"

    for F in "$COUNTRY_DIR"/*.F.fasta.gz; do
        [[ -e "$F" ]] || continue

        R="${F/.F.fasta.gz/.R.fasta.gz}"
        SAMPLE=$(basename "${F%.F.fasta.gz}")

        if [[ ! -f "$R" ]]; then
            echo "WARNING: missing reverse file for $COUNTRY/$SAMPLE, skipping"
            continue
        fi

        echo "Trimming $COUNTRY / $SAMPLE"
        zcat -f "$F" "$R" \
            | cutadapt -u "$TRIM5" -j "$THREADS" - \
                -o "$OUT_TRIMMED/$COUNTRY/${SAMPLE}.trim${TRIM5}.fa"
    done
done

echo "Done. Trimmed files saved in: $OUT_TRIMMED"
