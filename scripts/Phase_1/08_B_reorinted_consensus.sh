#!/usr/bin/env bash

# 08_B_reoriented_consensus.sh — Select correct contig and reorient by 16S strand
# Input:  RRN_CONSENSUS_DIR/{species}/assembly.fasta + annotation.gff
# Output: RRN_CONSENSUS_DIR/{species}/assembly_reoriented.fasta
#
# Special cases (manually resolved multi-contig assemblies):
#   Alistipes_shahii      → contig_2
#   Bacteroides_ovatus    → contig_1
#   Streptococcus_salivarius → contig_2 (reverse complement)
#
# Usage: bash 09_B_reoriented_consensus.sh

set -euo pipefail
source "$(dirname "$0")/../config.sh"

for D in "$RRN_CONSENSUS_DIR"/*/; do
    [[ -d "$D" ]] || continue
    GFF="$D/annotation.gff"
    FA="$D/assembly.fasta"
    OUT="$D/assembly_reoriented.fasta"
    [[ -f "$GFF" && -f "$FA" ]] || continue

    SP=$(basename "$D")

    # Special cases: manually selected contigs
    if [[ "$SP" == "Alistipes_shahii" ]]; then
        seqkit grep -p contig_2 "$FA" > "$OUT"
        echo "$SP → contig_2"
        continue
    elif [[ "$SP" == "Bacteroides_ovatus" ]]; then
        seqkit grep -p contig_1 "$FA" > "$OUT"
        echo "$SP → contig_1"
        continue
    elif [[ "$SP" == "Streptococcus_salivarius" ]]; then
        seqkit grep -p contig_2 "$FA" | seqkit seq -t DNA -r -p - > "$OUT"
        echo "$SP → contig_2 (reverse complement)"
        continue
    fi

    # General case: select contig containing 16S_rRNA and reorient if on minus strand
    CONTIG=$(awk '/16S_rRNA/ {print $1; exit}' "$GFF")
    STRAND=$(awk '/16S_rRNA/ {print $7; exit}' "$GFF")

    seqkit grep -p "$CONTIG" "$FA" > "$D/tmp.fa"
    if [[ "$STRAND" == "-" ]]; then
        seqkit seq -t DNA -r -p "$D/tmp.fa" > "$OUT"
    else
        mv "$D/tmp.fa" "$OUT"
    fi
    echo "$SP → $CONTIG ($STRAND)"
done

echo "Done."
