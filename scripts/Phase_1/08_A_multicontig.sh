#!/usr/bin/env bash

# 08_A_multicontig.sh — Report assemblies with more than one contig
# Input:  RRN_CONSENSUS_DIR/{species}/assembly_info.txt  (from 07_rrn_consensus_novo.sh)
# Output: printed to stdout
#
# Usage: bash 09_A_multicontig.sh

set -euo pipefail
source "$(dirname "$0")/../config.sh"

for F in "$RRN_CONSENSUS_DIR"/*/assembly_info.txt; do
    [[ -f "$F" ]] || continue
    N=$(grep -vc '^#' "$F")
    (( N > 1 )) && echo "$(dirname "$F"): ${N} contigs"
done
