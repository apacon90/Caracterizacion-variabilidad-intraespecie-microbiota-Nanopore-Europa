#!/usr/bin/env bash
# 02_map_filter.sh — Map trimmed reads against RefSeq207nr and filter alignments
# Input:  OUT_TRIMMED/{country}/{sample}.trim{N}.fa  (from 01_trimm.sh)
# Output: MAPPINGS_DIR/{sample}.keep.paf
# Filters: alignment length >= MIN_ALN, ref coverage >= MIN_COV, identity >= MIN_PID

set -euo pipefail
export LC_NUMERIC=C
source "$(dirname "$0")/config.sh"

mkdir -p "$MAPPINGS_DIR"

VALID_SAMPLES=$(tail -n +2 "$META_FILE" | cut -f1)

# Build minimap2 index if missing
if [[ ! -f "$DB_MMI" ]]; then
    echo "Building minimap2 index: $DB_MMI"
    minimap2 -d "$DB_MMI" "$DB_REF"
fi

for FILE in $(find "$OUT_TRIMMED" -type f -name "*.trim${TRIM5}.fa" | sort); do
    SAMPLE=$(basename "$FILE" .trim${TRIM5}.fa)
    OUT_KEEP="$MAPPINGS_DIR/${SAMPLE}.keep.paf"

    if ! echo "$VALID_SAMPLES" | grep -qx "$SAMPLE"; then
        echo "Skipping $SAMPLE (not in metadata)"
        continue
    fi

    if [[ -f "$OUT_KEEP" ]]; then
        echo "Skipping $SAMPLE (already done)"
        continue
    fi

    echo "Mapping $SAMPLE"
    START=$(date +%s)

    # Map with minimap2, filter on-the-fly with awk, append coverage/identity tags
    minimap2 -t "$THREADS" -x map-ont -c --secondary=no "$DB_MMI" "$FILE" \
        | awk -v min_aln="$MIN_ALN" -v min_cov="$MIN_COV" -v min_pid="$MIN_PID" \
            'BEGIN { FS = OFS = "\t" }
            {
                pid = $10 / $11
                cov = $11 / $7
                if ($11 >= min_aln && cov >= min_cov && pid >= min_pid) {
                    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tcr:f:%.6f\tpi:f:%.6f\n",
                           $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,cov,pid
                }
            }' > "${OUT_KEEP}.tmp"

    mv -f "${OUT_KEEP}.tmp" "$OUT_KEEP"
    echo "Done $SAMPLE ($(( $(date +%s) - START ))s)"
done

echo "All .keep.paf files saved in: $MAPPINGS_DIR"
