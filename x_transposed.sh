#!/bin/bash
set -euo pipefail

#########################################
# This script identifies X-transposed regions by aligning chrY to chrX using LASTZ, 
# filtering for high-identity (>94%) and long (>1 kb) alignments, then removes overlaps 
# with known PAR and ampliconic regions using bedtools
#########################################

# Base directory
BASE_DIR="$(pwd)"

# Input files
CHR_Y_FA="$BASE_DIR/PAR_annotation/combined_mask_fastas/softmasked_chrY.fa"
CHR_X_FA="$BASE_DIR/PAR_annotation/combined_mask_fastas/softmasked_chrX.fa"
non_PAR_ANNOT_DIR="$BASE_DIR/non_PAR_annotation"
XTRAN_OUTPUT_DIR="$non_PAR_ANNOT_DIR/x_transposed_output"

# Create output directory
mkdir -p "$XTRAN_OUTPUT_DIR"

# Filtering parameters
MIN_IDENTITY=94
MIN_LENGTH=1000

#########################################
# LASTZ Alignment (Y vs X)
#########################################

echo "Running LASTZ alignment between X and Y chromosomes"
RAW_OUTPUT="${XTRAN_OUTPUT_DIR}/marmX_to_marmY_alignments.tsv"

lastz "$CHR_Y_FA" "$CHR_X_FA" \
    --format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,id%,cigarx \
    --identity="$MIN_IDENTITY" \
    --output="$RAW_OUTPUT"

#########################################
# Process Both Chromosomes
#########################################

for CHROM in X Y; do
    echo "Processing chromosome $CHROM"
    
    # Generate chromosome-specific paths
    PAR_BED="$BASE_DIR/PAR_annotation/PAR_region/${CHROM}_PAR_region.bed"
    AMP_BED="$non_PAR_ANNOT_DIR/amplicon_output/${CHROM}_pri_amplicons_final.bed"
    UNION_REGIONS="${XTRAN_OUTPUT_DIR}/chr${CHROM}_union_regions.bed"
    FINAL_BED="${XTRAN_OUTPUT_DIR}/chr${CHROM}_filtered_alignments.bed"
    
    # Extract chromosome-specific alignments with proper formatting
    FILTERED_INITIAL="${XTRAN_OUTPUT_DIR}/chr${CHROM}_unfiltered.bed"
    
    if [ "$CHROM" = "Y" ]; then
        # Process Y chromosome coordinates (columns 1-3)
        awk -F'\t' -v min_len="$MIN_LENGTH" '
            BEGIN {OFS="\t"}
            /^#/ {next}  # Skip header lines
            NF >= 3 && ($3 - $2) >= min_len {print $1, $2, $3}
        ' "$RAW_OUTPUT" > "$FILTERED_INITIAL"
    else
        # Process X chromosome coordinates (columns 4,6,7)
        awk -F'\t' -v min_len="$MIN_LENGTH" '
            BEGIN {OFS="\t"}
            /^#/ {next}  # Skip header lines
            NF >= 7 && ($3 - $2) >= min_len {print $4, $6, $7}
        ' "$RAW_OUTPUT" > "$FILTERED_INITIAL"
    fi

    # Create union of exclusion regions with validation
    EXCLUSIONS=()
    [ -f "$PAR_BED" ] && EXCLUSIONS+=("$PAR_BED")
    
    # Add ampliconic regions only for Y chromosome
    if [ "$CHROM" = "Y" ] && [ -f "$AMP_BED" ]; then
        EXCLUSIONS+=("$AMP_BED")
    fi

    if [ ${#EXCLUSIONS[@]} -eq 0 ]; then
        echo "No exclusion files found for chr$CHROM, using raw alignments"
        cp "$FILTERED_INITIAL" "$FINAL_BED"
    else
        echo "Creating union of ${#EXCLUSIONS[@]} exclusion files for chr$CHROM"
        # Process exclusion files with header skipping
        for file in "${EXCLUSIONS[@]}"; do
            awk 'BEGIN {OFS="\t"} /^[^#]/ {print $1, $2, $3}' "$file"
        done | bedtools sort | bedtools merge > "$UNION_REGIONS"
        
        # Run bedtools subtract with sorted inputs
        bedtools subtract \
            -a <(bedtools sort -i "$FILTERED_INITIAL") \
            -b "$UNION_REGIONS" \
            > "$FINAL_BED"
    fi

    echo "chr$CHROM filtered regions saved to: $FINAL_BED"
done

echo "All chromosome processing complete"
