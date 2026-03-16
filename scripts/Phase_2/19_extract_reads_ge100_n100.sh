#!/usr/bin/env bash

# 19_extract_reads_ge100_n100.sh — Extract reads per individual and species, subsample to N=100
# Input:  PAIRS_FILE (pairs_14bacteria_samples_ge100.tsv, from Phase 0)
#         MAPPINGS_DIR/*.keep.paf
#         OUT_TRIMMED/{country}/{sample}.trim{N}.fa
# Output: FASTA_ALL/{species}/{country}/{species}.{sample}.all.fa
#         FASTA_N/{species}/{country}/{species}.{sample}.n100.fa
#
# Usage: bash 19_extract_reads_ge100_n100.sh

set -euo pipefail
source "$(dirname "$0")/../config.sh"

N=100
SEED=42

mkdir -p "$FASTA_ALL" "$FASTA_N" "$AUX_DIR"

command -v seqkit >/dev/null 2>&1 || { echo "ERROR: seqkit is not installed" >&2; exit 1; }

# Reuse ID -> species dictionary from Phase 1
ID2SP="${SPECIES_READS_DIR}/aux/id_to_species.tsv"
[[ -f "$ID2SP" ]] || { echo "ERROR: id_to_species.tsv not found at $ID2SP — run Phase 1 first" >&2; exit 1; }

find_fasta() {
    local sample="$1"
    ls "$OUT_TRIMMED"/*/"${sample}${READS_SUFFIX}" 2>/dev/null | head -n 1 || true
}

tail -n +2 "$PAIRS_FILE" | while IFS=$'\t' read -r SPEC SAMPLE COUNT; do
    SAMPLE=${SAMPLE//./-}

    PAF="$MAPPINGS_DIR/${SAMPLE}.keep.paf"
    RFILE=$(find_fasta "$SAMPLE")
    COUNTRY=$(basename "$(dirname "$RFILE")")

    mkdir -p "$FASTA_ALL/$SPEC/$COUNTRY" "$FASTA_N/$SPEC/$COUNTRY"

    ALLFA="$FASTA_ALL/$SPEC/$COUNTRY/${SPEC}.${SAMPLE}.all.fa"
    OUTFA="$FASTA_N/$SPEC/$COUNTRY/${SPEC}.${SAMPLE}.n${N}.fa"
    QN="$AUX_DIR/${SPEC}.${SAMPLE}.qnames"

    # Get RefSeq IDs for this species and extract matching read names from PAF
    awk -F'\t' -v s="$SPEC" '$2==s{print $1}' "$ID2SP" | \
    awk -v FS='\t' 'NR==FNR{ keep[$1]=1; next } ($6 in keep){ print $1 }' - "$PAF" \
    | sort -u > "$QN"

    # Extract reads from FASTA
    awk 'NR==FNR { want[$1]=1; next }
         /^>/ { h=$0; sub(/^>/,"",h); split(h,a,/ /); p=(a[1] in want) }
         p { print }' "$QN" "$RFILE" > "$ALLFA"
    rm -f "$QN"

    NREADS=$(grep -c '^>' "$ALLFA" || true)
    if [[ "$NREADS" -eq 0 ]]; then
        echo "ERROR: $SPEC / $SAMPLE -> 0 reads extracted" >&2
        rm -f "$ALLFA"
        continue
    fi

    # Subsample to exactly N (all individuals in PAIRS_FILE have >= N reads by construction)
    set +o pipefail
    seqkit shuffle -s "$SEED" "$ALLFA" | seqkit head -n "$N" > "$OUTFA"
    set -o pipefail

    NOUT=$(grep -c '^>' "$OUTFA" || true)
    echo "OK: $SPEC / $SAMPLE (OTU=$COUNT) -> all=$NREADS out=$NOUT"
done

echo "Done."
echo "All reads:  $FASTA_ALL"
echo "Subsampled: $FASTA_N"
