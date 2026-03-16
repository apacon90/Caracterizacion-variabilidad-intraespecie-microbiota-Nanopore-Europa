#!/usr/bin/env bash
# 08_rrn_consensus_novo.sh — De novo assembly of ribosomal operon with Flye + Barrnap annotation
# Input:  SPECIES_READS_DIR/combined_1500/{species}.sampled.fa  (from 07_sample_sp_1500_reads.sh)
# Output: RRN_CONSENSUS_DIR/{species}/assembly.fasta
#         RRN_CONSENSUS_DIR/{species}/annotation.gff
#
# Usage: bash 08_rrn_consensus_novo.sh

set -euo pipefail
source "$(dirname "$0")/../config.sh"

IN_DIR="${SPECIES_READS_DIR}/combined_1500"
OUT_DIR="${RRN_CONSENSUS_DIR}"
mkdir -p "$OUT_DIR"

for FA in "$IN_DIR"/*.sampled.fa; do
    SP=$(basename "$FA" .sampled.fa)
    SP_DIR="${OUT_DIR}/${SP}"
    rm -rf "$SP_DIR"
    mkdir -p "$SP_DIR"

    echo "Assembling $SP"
    flye \
        --nano-raw "$FA" \
        --out-dir "$SP_DIR" \
        --threads "$THREADS" \
        --genome-size 6k \
        --asm-coverage 200 \
        --min-overlap 1200 \
        > "${SP_DIR}/flye.log" 2>&1

    if [[ -s "${SP_DIR}/assembly.fasta" ]]; then
        echo "Annotating $SP"
        "$BARRNAP" --kingdom bac --threads "$THREADS" "${SP_DIR}/assembly.fasta" \
            > "${SP_DIR}/annotation.gff"
    else
        echo "WARNING: no assembly.fasta for $SP — check ${SP_DIR}/flye.log"
    fi
done

echo "Done. Assemblies and annotations saved in: $OUT_DIR"
