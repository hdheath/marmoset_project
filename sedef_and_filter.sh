#!/bin/bash
# --------------------------------------------------------------------
# Sedef Analysis and Filtering Script
# This script runs Sedef on a combined softmasked FASTA and then filters
# the Sedef results for candidate PAR regions.
#
# Usage: ./sedef_and_filter.sh
# Ensure that the combined softmasked FASTA file is available in the expected directory.
# --------------------------------------------------------------------

# Base directory and directory definitions (should match those used in mask_and_merge.sh)
BASE_DIR="$(pwd)"
PAR_ANNOT_DIR="$BASE_DIR/PAR_annotation"
COMBINED_DIR="$PAR_ANNOT_DIR/combined_mask_fastas"
SEDEF_OUTPUT_DIR="$PAR_ANNOT_DIR/sedef_output"
LOG_DIR="$PAR_ANNOT_DIR/logs"
SEDEF_FINAL_DIR="$SEDEF_OUTPUT_DIR/filtered_results"

# Check if the combined softmasked FASTA file exists
if [ ! -f "$COMBINED_DIR/combined_softmasked.fa" ]; then
  echo "Error: Combined softmasked FASTA file not found at $COMBINED_DIR/combined_softmasked.fa"
  exit 1
fi

# Prompt user for the Sedef script path
read -p "Enter the full path to sedef.sh: " SEDEF_PATH

# Check if the file exists and is executable
if [ ! -f "$SEDEF_PATH" ] || [ ! -x "$SEDEF_PATH" ]; then
  echo "Error: The specified sedef.sh script does not exist or is not executable."
  exit 1
fi

# Create necessary directories for Sedef output and filtered results
mkdir -p "$SEDEF_OUTPUT_DIR" "$SEDEF_FINAL_DIR"

# --------------------------------------------------------------------
# Run Sedef on the Combined Softmasked FASTA
# --------------------------------------------------------------------
echo "Running Sedef with 40 threads..."
"$SEDEF_PATH" -o "$SEDEF_OUTPUT_DIR" -j 40 "$COMBINED_DIR/combined_softmasked.fa"
echo "Sedef analysis completed. Results saved in $SEDEF_OUTPUT_DIR."

# --------------------------------------------------------------------
# Filter Sedef Results for Candidate PAR Regions
# --------------------------------------------------------------------
echo "Filtering Sedef results for candidate PAR regions..."

# Extract only chrX and chrY alignments
SEDEF_PAR_RAW="$SEDEF_FINAL_DIR/chrX_chrY_raw_alignments.bed"
grep 'chrX' "$SEDEF_OUTPUT_DIR/final.bed" | grep 'chrY' > "$SEDEF_PAR_RAW"
echo "Filtered chrX and chrY alignments saved as: $SEDEF_PAR_RAW"

# Process the filtered file to extract relevant columns
SEDEF_PAR_PROCESSED="$SEDEF_FINAL_DIR/chrX_chrY_parsed_alignments.bed"
awk '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $11 "\t" $12 "\t" $21 "\t" $22 }' \
  "$SEDEF_PAR_RAW" | sort -k1,1 -k2,2n > "$SEDEF_PAR_PROCESSED"
echo "Processed PAR candidate alignments saved as: $SEDEF_PAR_PROCESSED"

# Sort by column 9 (descending order) and save
SEDEF_PAR_SORTED="$SEDEF_FINAL_DIR/chrX_chrY_sorted_by_similarity.bed"
sort -k9,9nr "$SEDEF_PAR_PROCESSED" > "$SEDEF_PAR_SORTED"
echo "Sorted PAR candidates (by sequence similarity) saved as: $SEDEF_PAR_SORTED"

# Save a copy specifically for plotting
SEDEF_PAR_PLOT="$SEDEF_FINAL_DIR/chrX_chrY_plot_ready.bed"
cp "$SEDEF_PAR_SORTED" "$SEDEF_PAR_PLOT"
echo "Plot-ready PAR candidate file saved as: $SEDEF_PAR_PLOT"

# Display the top 5 PAR candidates
echo -e "\nTop 5 PAR Candidates (Sorted by Similarity Score):"
head -n 5 "$SEDEF_PAR_SORTED"

echo -e "\nProcessing completed. Final PAR candidate results are available in: $SEDEF_FINAL_DIR"
