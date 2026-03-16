#!/usr/bin/env bash

# 22_filter_qual_individuals.sh — Filter variants by quality (SNPs QUAL>=100, INDELs QUAL>=150)
# Input:  IND_VCF_DIR/{species}.vcf.gz  (from 20_call_individuals_VCF)
#         REF_DIR/{species}.fa  (from Phase 1)
# Output: IND_VCF_FILT_DIR/{species}.QUAL_100-150.vcf.gz + .tbi
#
# Usage: bash 22_filter_qual_individuals.sh

set -euo pipefail
shopt -s nullglob
source "$(dirname "$0")/../config.sh"


mkdir -p "$IND_VCF_FILT_DIR"

for VCF in "$IND_VCF_DIR"/*.vcf.gz; do
    SPECIES=$(basename "$VCF" .vcf.gz)
    REF="$REF_DIR/${SPECIES}.fa"
    OUT="$IND_VCF_FILT_DIR/${SPECIES}.QUAL_${QUAL_SNP}-${QUAL_INDEL}.vcf.gz"

    echo "Filtering $SPECIES"

    bcftools norm -f "$REF" -m -both "$VCF" -Ou \
    | bcftools view -i "(TYPE=\"snp\" && QUAL>=$QUAL_SNP) || (TYPE=\"indel\" && QUAL>=$QUAL_INDEL)" \
        -Oz -o "$OUT"

    tabix -f -p vcf "$OUT"
done

echo "Done. Filtered VCFs saved in: $IND_VCF_FILT_DIR"
