#!/bin/bash
# --------------------------------------------------------------------
# Palindrome Detection Script
# This script performs self-alignments using LASTZ on a masked FASTA
# and then runs the palindrover.py tool to identify palindromes.
#
# Usage: ./palindrome_regions.sh
# You will be prompted for the full path to palindrover.py.
# --------------------------------------------------------------------

# Base directory is the current working directory
BASE_DIR="$(pwd)"

# Define the directory for non-PAR annotations and the palindrome output directory
non_PAR_ANNOT_DIR="$BASE_DIR/non_PAR_annotation"
PAL_OUTPUT_DIR="$non_PAR_ANNOT_DIR/palindrome_output"

# Create the palindrome output directory if needed
mkdir -p "$non_PAR_ANNOT_DIR" "$PAL_OUTPUT_DIR"

# Prompt user for the full path to palindrover.py
read -p "Enter the full path to palindrover.py: " PAL_PATH

# Check if the palindrover script exists
if [ ! -f "$PAL_PATH" ]; then
  echo "Error: The specified palindrover.py script does not exist."
  exit 1
fi

# Define masked FASTA file paths (adjust if processing chrX as well)
MASKED_Y_CHROMOSOME="$BASE_DIR/PAR_annotation/combined_mask_fastas/softmasked_chrY.fa"

# Check that the masked FASTA file exists
if [ ! -f "$MASKED_Y_CHROMOSOME" ]; then
  echo "Error: Masked Y chromosome FASTA not found at:"
  echo "       $MASKED_Y_CHROMOSOME"
  exit 1
fi

# Define temporary alignment output (for LASTZ)
ALIGN_OUT_Y="$PAL_OUTPUT_DIR/chrY_alignments.tsv"

# Run LASTZ self-alignment on the Y chromosome
echo "Running LASTZ self-alignment on Y chromosome..."
lastz "$MASKED_Y_CHROMOSOME" "$MASKED_Y_CHROMOSOME" \
  --format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,id%,cigarx \
  --identity=90 \
  --output="$ALIGN_OUT_Y"

if [ ! -s "$ALIGN_OUT_Y" ]; then
  echo "Error: LASTZ did not produce any output for the Y chromosome."
  exit 1
fi

# Run palindrover.py to detect palindromes
echo "Detecting palindromes using palindrover.py..."
python "$PAL_PATH" \
  --minidentity=0.9 \
  --minlength=100000 \
  --output="$PAL_OUTPUT_DIR/chrY_palindromes.bed" < "$ALIGN_OUT_Y"

if [ -s "$PAL_OUTPUT_DIR/chrY_palindromes.bed" ]; then
  echo "Palindrome detection complete. Output:"
  echo "$PAL_OUTPUT_DIR/chrY_palindromes.bed"
else
  echo "Warning: No palindromes detected on the Y chromosome."
fi

# Repeat similar steps for chrX if needed.
