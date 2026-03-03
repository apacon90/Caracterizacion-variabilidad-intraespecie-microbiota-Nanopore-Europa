#!/usr/bin/env bash
# 09_blast_match_top5.sh — BLASTn top5 hits for each assembled consensus vs RefSeq207nr
# Input:  RRN_CONSENSUS_DIR/{species}/assembly_reoriented.fasta  (from 08_B)
# Output: RESULTS_DIR/blast_top5_taxrep.tsv
#
# Builds BLAST DB if missing, then runs megablast (top5, evalue<1e-10).
# Usage: bash 09_blast_match_top5.sh

set -euo pipefail
source "$(dirname "$0")/config.sh"

BLAST_DB="${RESULTS_DIR}/blast_db/refseq207nr_nrRep_db"
OUT="${RESULTS_DIR}/blast_top5_taxrep.tsv"

mkdir -p "$(dirname "$BLAST_DB")"

# Build BLAST database if missing
if [[ ! -f "${BLAST_DB}.nin" ]]; then
    echo "Building BLAST database..."
    makeblastdb -in "$DB_REF_FASTA" -dbtype nucl -out "$BLAST_DB" -parse_seqids
fi

echo -e "query\tcode\tpident\tlength\tqcovs\ttaxon_full" > "$OUT"

for F in "$RRN_CONSENSUS_DIR"/*/assembly_reoriented.fasta; do
    SP=$(basename "$(dirname "$F")")
    echo "BLASTing $SP"

    blastn -query "$F" -db "$BLAST_DB" \
           -task megablast \
           -outfmt '6 sseqid pident length qcovs' \
           -max_target_seqs 5 \
           -max_hsps 1 \
           -evalue 1e-10 \
           -num_threads "$THREADS" \
    | while read -r CODE PIDENT LENGTH QCOVS; do
        TAX=$(grep -m1 "^${CODE}"$'\t' "$TAX_FILE" | cut -f2)
        [[ -z "$TAX" ]] && TAX="UNKNOWN"
        echo -e "${SP}\t${CODE}\t${PIDENT}\t${LENGTH}\t${QCOVS}\t${TAX}" >> "$OUT"
    done
done

echo "Done. Results saved in: $OUT"
