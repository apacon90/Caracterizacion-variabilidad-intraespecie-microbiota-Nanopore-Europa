#!/usr/bin/env bash
# 06_sample_sp_1500_reads.sh — Subsample 500 reads per species per country (1500 total)
# Input:  SPECIES_READS_DIR/{country}/{species}.reads.fa  (from 05_reads_species.sh)
# Output: SPECIES_READS_DIR/combined_1500/{species}.sampled.fa
#
# Uses reservoir sampling to randomly select N reads per country file (seed fixed).
# Usage: bash 06_sample_sp_1500_reads.sh

set -euo pipefail
source "$(dirname "$0")/config.sh"

OUT_DIR="${SPECIES_READS_DIR}/combined_1500"
mkdir -p "$OUT_DIR"

N=500      # reads to sample per country
SEED=42    # fixed seed for reproducibility

for SPEC_FILE in $(find "$SPECIES_READS_DIR" -type f -name "*.reads.fa" -printf "%f\n" | sort -u); do
    SPEC="${SPEC_FILE%.reads.fa}"
    OUT_SPEC="${OUT_DIR}/${SPEC}.sampled.fa"
    : > "$OUT_SPEC"

    for INFILE in $(find "$SPECIES_READS_DIR" -type f -name "${SPEC}.reads.fa" | sort); do
        awk -v N="$N" -v SEED="$SEED" '
            BEGIN { RS=">"; ORS=""; srand(SEED); k=0 }
            NR==1 { next }
            {
                rec=">"$0; body=$0; sub(/^[^\n]*\n?/,"",body); gsub(/\n/,"",body);
                if (length(body)==0) next;
                k++;
                if (k<=N) { res[k]=rec }
                else { j=int(rand()*k)+1; if (j<=N) res[j]=rec }
            }
            END { m=(k<N?k:N); for(i=1;i<=m;i++) printf "%s", res[i] }
        ' "$INFILE" >> "$OUT_SPEC" || true
    done

    [[ -s "$OUT_SPEC" ]] || rm -f "$OUT_SPEC"
done

echo "Done. Sampled FASTAs saved in: $OUT_DIR"
