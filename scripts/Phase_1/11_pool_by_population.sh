#!/usr/bin/env bash
# 11_pool_by_population.sh — Split reads into N pools of fixed size per species and country
# Input:  SPECIES_READS_DIR/{country}/{species}.reads.fa  (from 05_reads_species.sh)
# Output: POOLS_DIR/{species}/{country}_pool{01..10}.fa
#
# Reads are shuffled (fixed seed) then split into consecutive blocks of POOL_SIZE.
# Usage: bash 11_pool_by_population.sh

set -euo pipefail
source "$(dirname "$0")/config.sh"

N_POOLS=10
POOL_SIZE=200
SEED=42

mkdir -p "$POOLS_DIR"

find "$SPECIES_READS_DIR" -mindepth 2 -maxdepth 2 -type f -name "*.reads.fa" | sort | while read -r INFA; do
    COUNTRY=$(basename "$(dirname "$INFA")")
    SPECIES=$(basename "$INFA" .reads.fa)
    OUT_DIR="$POOLS_DIR/$SPECIES"
    mkdir -p "$OUT_DIR"

    echo "Pooling $SPECIES ($COUNTRY)"

    TMP=$(mktemp)
    trap 'rm -f "$TMP"' EXIT

    seqkit shuffle -s "$SEED" "$INFA" > "$TMP"

    awk -v RS=">" -v ORS="" \
        -v n="$POOL_SIZE" -v p="$N_POOLS" \
        -v country="$COUNTRY" -v od="$OUT_DIR" '
        NR==1 { next }
        {
            rec=">"$0; k++; pool=int((k-1)/n)+1;
            if (pool<=p) {
                fn=sprintf("%s/%s_pool%02d.fa", od, country, pool);
                print rec >> fn
            }
        }' "$TMP"

    rm -f "$TMP"
done

echo "Done. Pools saved in: $POOLS_DIR"
