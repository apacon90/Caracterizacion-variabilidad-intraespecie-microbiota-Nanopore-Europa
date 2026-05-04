#!/usr/bin/env bash
# 04_reads_species_matrix.sh — Count unique reads per individual × species from PAF files
# Input:  MAPPINGS_DIR/*.keep.paf, CODE2NAME_FILE
# Output: RESULTS_DIR/reads_by_species_matrix_unique.tsv
# Usage:  bash 04_reads_species_matrix.sh

set -euo pipefail
source "$(dirname "$0")/../config.sh"

mkdir -p "$TABLES_DIR"
OUT="${TABLES_DIR}/reads_by_species_matrix_unique.tsv"

awk -F'\t' '
NR==FNR { m[$1]=$2; next }
{
    f=FILENAME; sub(/^.*\//,"",f); sub(/\.keep\.paf$/,"",f); samp[f]=1
    sp = (m[$6] ? m[$6] : "UNKNOWN")
    key = f SUBSEP sp SUBSEP $1
    if (!seen[key]++) { c[sp,f]++ }
    species[sp]=1
}
END {
    printf "Species"
    for (s in samp) printf "\t%s", s
    print ""
    for (sp in species) {
        printf "%s", sp
        for (s in samp) printf "\t%d", c[sp,s]+0
        print ""
    }
}' "$CODE2NAME_FILE" "$MAPPINGS_DIR"/*.keep.paf > "$OUT"

echo "Done -> $OUT"
