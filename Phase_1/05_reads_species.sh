#!/usr/bin/env bash
# 04_reads_species.sh — Extract reads per species and country from PAF alignments
# Input:  MAPPINGS_DIR/*.keep.paf + OUT_TRIMMED/{country}/{sample}.trim{N}.fa
# Output: SPECIES_READS_DIR/{country}/{species}.reads.fa
#
# Requires: GNU parallel
# Usage: bash 05_reads_species.sh

set -euo pipefail
source "$(dirname "$0")/config.sh"

AUX_DIR="${SPECIES_READS_DIR}/aux"
mkdir -p "$SPECIES_READS_DIR" "$AUX_DIR"

command -v parallel >/dev/null 2>&1 || { echo "ERROR: GNU parallel is not installed" >&2; exit 1; }

# Step 1: Build ID -> species dictionary for Top50 species only
awk -F'\t' '
    FNR==1 && NR==FNR { for(i=1;i<=NF;i++) if($i=="Species") c=i; next }
    NR==FNR { if($c!="") keep[$c]=1; next }
    BEGIN { OFS="\t" }
    {
        id=$1; tax=$2; sp="";
        n=split(tax,a,"|");
        for(i=1;i<=n;i++) if(a[i]~/^s__/){ sub(/^s__/,"",a[i]); sp=a[i]; break }
        gsub(/\[|\]/,"",sp);
        gsub(/[[:space:]]+|[\/]/,"_",sp);
        gsub(/_+/,"_",sp);
        sub(/^_+|_+$/,"",sp);
        if(sp!="" && (sp in keep)) print id, sp
    }
' "$TOP50_TSV" "$TAX_FILE" > "$AUX_DIR/id_to_species.tsv"

# Step 2+3: For each species, find read names in PAFs and extract FASTAs by country
process_one() {
    SPEC="$1"
    [ -n "$SPEC" ] || return 0

    # Get RefSeq IDs belonging to this species
    IDS=$(awk -F'\t' -v s="$SPEC" '$2==s{print $1}' "$AUX_DIR/id_to_species.tsv") || true
    [ -n "$IDS" ] || return 0

    # Collect read names (qnames) that map to this species across all PAF files
    QNAMES="$AUX_DIR/.${SPEC}.qnames.txt"
    : > "$QNAMES"
    for PAF in "$MAPPINGS_DIR"/*.keep.paf; do
        [ -f "$PAF" ] || continue
        awk -v FS='\t' -v ids="$(printf "%s\n" "$IDS")" '
            BEGIN{ split(ids,a,"\n"); for(i in a) if(a[i]!="") keep[a[i]]=1 }
            ($6 in keep){ print $1 }
        ' "$PAF" >> "$QNAMES"
    done
    [ -s "$QNAMES" ] || { rm -f "$QNAMES"; return 0; }

    # Extract FASTA sequences organized by country
    for PAF in "$MAPPINGS_DIR"/*.keep.paf; do
        [ -f "$PAF" ] || continue
        SAMPLE=$(basename "$PAF" .keep.paf)
        for RFILE in "$OUT_TRIMMED"/*/"${SAMPLE}${READS_SUFFIX}"; do
            [ -f "$RFILE" ] || continue
            COUNTRY=$(basename "$(dirname "$RFILE")")
            mkdir -p "$SPECIES_READS_DIR/$COUNTRY"
            awk 'NR==FNR { want[$1]=1; next }
                 /^>/    { h=$0; sub(/^>/,"",h); split(h,a,/ /); p=(a[1] in want) }
                 p       { print }' \
                "$QNAMES" "$RFILE" >> "$SPECIES_READS_DIR/$COUNTRY/${SPEC}.reads.fa"
        done
    done

    rm -f "$QNAMES"
}

export -f process_one
export AUX_DIR MAPPINGS_DIR OUT_TRIMMED READS_SUFFIX SPECIES_READS_DIR

awk -F'\t' 'NR==1{ for(i=1;i<=NF;i++) if($i=="Species") c=i; next } { if($c!="") print $c }' "$TOP50_TSV" \
    | parallel -j "$JOBS" process_one

echo "Done. FASTAs saved in: $SPECIES_READS_DIR"
