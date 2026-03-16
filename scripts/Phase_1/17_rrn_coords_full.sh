#!/bin/bash
# 17_rrn_coords_full.sh — Extract ribosomal region coordinates from RefSeq references
# Input:  REF_DIR/{species}.fa           (from 13_mapp_call.sh)
# Output: RESULTS_DIR/master_rrn_coords_completo.tsv
#
# Runs barrnap on each species-specific RefSeq reference and extracts
# coordinates of 16S_rRNA, 23S_rRNA and ITS regions.
# Contig column = RefSeq code (e.g. C_00025396), matching CHROM in VCF files.
#
# Usage: bash 17_rrn_coords_full.sh

set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUTPUT="${RESULTS_DIR}/master_rrn_coords_completo.tsv"
TMP_GFF="${RESULTS_DIR}/tmp_barrnap.gff"

echo -e "Species\tRegion\tContig\tStart\tEnd" > "$OUTPUT"
echo "Procesando referencias en: $REF_DIR"

for fa in "$REF_DIR"/*.fa; do
    [[ -f "$fa" ]] || continue
    species=$(basename "$fa" .fa)

    # Get RefSeq contig code from FASTA header
    contig=$(grep "^>" "$fa" | head -1 | sed 's/^>//' | awk '{print $1}')

    # Run barrnap
    "$BARRNAP" --kingdom bac --threads "$THREADS" "$fa" > "$TMP_GFF" 2>/dev/null

    # Extract coordinates
    awk -v sp="$species" -v ctg="$contig" -v out="$OUTPUT" '
        BEGIN { OFS="\t"; end16S=0; start23S=0 }
        $3 == "rRNA" {
            match($9, /Name=([^;]+)/, arr);
            gene = arr[1];
            print sp, gene, ctg, $4, $5 >> out;
            if (gene == "16S_rRNA") {
                if ($5 > end16S) { end16S = $5 }
            }
            else if (gene == "23S_rRNA") {
                if (start23S == 0 || $4 < start23S) { start23S = $4 }
            }
        }
        END {
            if (end16S > 0 && start23S > 0 && start23S > end16S) {
                print sp, "ITS", ctg, end16S+1, start23S-1 >> out;
            }
        }
    ' "$TMP_GFF"

    echo "Done: $species -> $contig"
done

rm -f "$TMP_GFF"
echo "Listo. Archivo generado: $OUTPUT"
