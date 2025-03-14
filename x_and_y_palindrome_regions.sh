#!/bin/bash
# --------------------------------------------------------------------
# Palindrome Detection Script
# This script performs self-alignments using LASTZ on masked FASTA
# files for both the Y and X chromosomes and then runs the palindrover.py tool
# to identify palindromes.
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

# Define masked FASTA file paths for Y and X chromosomes
MASKED_Y_CHROMOSOME="$BASE_DIR/PAR_annotation/combined_mask_fastas/softmasked_chrY.fa"
MASKED_X_CHROMOSOME="$BASE_DIR/PAR_annotation/combined_mask_fastas/softmasked_chrX.fa"

# Check that the masked Y FASTA file exists
if [ ! -f "$MASKED_Y_CHROMOSOME" ]; then
  echo "Error: Masked Y chromosome FASTA not found at:"
  echo "       $MASKED_Y_CHROMOSOME"
  exit 1
fi

# Check that the masked X FASTA file exists (if processing chrX)
if [ ! -f "$MASKED_X_CHROMOSOME" ]; then
  echo "Warning: Masked X chromosome FASTA not found at:"
  echo "       $MASKED_X_CHROMOSOME"
  echo "       Skipping X chromosome palindrome detection."
  PROCESS_X=false
else
  PROCESS_X=true
fi

# Define temporary alignment output for the Y chromosome
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

# Run palindrover.py to detect palindromes on Y chromosome alignments
echo "Detecting palindromes on Y chromosome using palindrover.py..."
python "$PAL_PATH" \
  --minidentity=0.9 \
  --minlength=100000 \
  --output="$PAL_OUTPUT_DIR/chrY_palindromes.bed" < "$ALIGN_OUT_Y"

if [ -s "$PAL_OUTPUT_DIR/chrY_palindromes.bed" ]; then
  echo "Palindrome detection complete for Y chromosome. Output:"
  echo "$PAL_OUTPUT_DIR/chrY_palindromes.bed"
else
  echo "Warning: No palindromes detected on the Y chromosome."
fi

# If the masked X FASTA exists, process the X chromosome similarly
if [ "$PROCESS_X" = true ]; then
  # Define temporary alignment output for the X chromosome
  ALIGN_OUT_X="$PAL_OUTPUT_DIR/chrX_alignments.tsv"

  # Run LASTZ self-alignment on the X chromosome
  echo "Running LASTZ self-alignment on X chromosome..."
  lastz "$MASKED_X_CHROMOSOME" "$MASKED_X_CHROMOSOME" \
    --format=general:name1,zstart1,end1,name2,strand2,zstart2+,end2+,id%,cigarx \
    --identity=90 \
    --output="$ALIGN_OUT_X"

  if [ ! -s "$ALIGN_OUT_X" ]; then
    echo "Error: LASTZ did not produce any output for the X chromosome."
    exit 1
  fi

  # Run palindrover.py to detect palindromes on X chromosome alignments
  echo "Detecting palindromes on X chromosome using palindrover.py..."
  python "$PAL_PATH" \
    --minidentity=0.9 \
    --minlength=100000 \
    --output="$PAL_OUTPUT_DIR/chrX_palindromes.bed" < "$ALIGN_OUT_X"

  if [ -s "$PAL_OUTPUT_DIR/chrX_palindromes.bed" ]; then
    echo "Palindrome detection complete for X chromosome. Output:"
    echo "$PAL_OUTPUT_DIR/chrX_palindromes.bed"
  else
    echo "Warning: No palindromes detected on the X chromosome."
  fi
fi

echo "Palindrome detection for specified chromosomes complete."
