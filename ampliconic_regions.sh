#!/bin/bash
# --------------------------------------------------------------------
# Ampliconic Regions Identification Script
# This script processes a masked chromosome FASTA to:
#   - Filter chromosome indexes using the FAI file.
#   - Create sliding windows.
#   - Extract FASTA sequences for BLAST.
#   - Run BLASTn in parallel.
#   - Merge BLAST hits and combine them with palindrome detections,
#     subtracting PAR and satellite regions.
#
# Usage: ./ampliconic_regions.sh
# This example processes the Y chromosome.
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# STEP 1: Set up dirs
# --------------------------------------------------------------------

# Base directory is the current working directory
BASE_DIR="$(pwd)"

# Define paths for PAR file, masked FASTA, and FAI index
PAR_FILE="$BASE_DIR/PAR_annotation/PAR_region/PAR.bed"
MASKED_Y_CHROMOSOME="$BASE_DIR/PAR_annotation/combined_mask_fastas/softmasked_chrY.fa"
FAI_Y_CHROMOSOME="$MASKED_Y_CHROMOSOME.fai"
REPEATMASKER_FILE="$BASE_DIR/PAR_annotation/repeatmasker_output/chrY.fa.out"

# Define non-PAR annotation directories
non_PAR_ANNOT_DIR="$BASE_DIR/non_PAR_annotation"
AMP_OUTPUT_DIR="$non_PAR_ANNOT_DIR/amplicon_output"
SAT_OUTPUT_DIR="$non_PAR_ANNOT_DIR/satellite_output"   # Created in satellite_regions.sh
PAL_OUTPUT_DIR="$non_PAR_ANNOT_DIR/palindrome_output"    # Created in palindrome_regions.sh

# Create output directory if needed
mkdir -p "$non_PAR_ANNOT_DIR" "$AMP_OUTPUT_DIR"

# Check if required files exist
if [ ! -f "$MASKED_Y_CHROMOSOME" ]; then
  echo "Error: Masked Y chromosome FASTA not found at: $MASKED_Y_CHROMOSOME"
  exit 1
fi

if [ ! -f "$FAI_Y_CHROMOSOME" ]; then
  echo "Error: FAI index for Y chromosome not found at: $FAI_Y_CHROMOSOME"
  exit 1
fi

if [ ! -f "$PAR_FILE" ]; then
  echo "Error: PAR file not found at: $PAR_FILE"
  exit 1
fi

# Set local chromosome tag for Y chromosome (adjust if processing chrX)
LOCAL_CHROMOSOME="chrY_pri"

# Create subdirectories for BLAST queries and results within AMP_OUTPUT_DIR
mkdir -p "${AMP_OUTPUT_DIR}/blast_queries" "${AMP_OUTPUT_DIR}/blast_results"

# --------------------------------------------------------------------
# STEP 2: Hardmask Repeats in the Masked FASTA
# --------------------------------------------------------------------
# (This step will replace selected repeat regions with Ns.)
#
# Check if the RepeatMasker file exists
if [ ! -f "$REPEATMASKER_FILE" ]; then
  echo "Error: RepeatMasker output file not found at: $REPEATMASKER_FILE"
  exit 1
fi

echo "Extracting satellite/repeat features for hardmasking..."
# Option 1: Filter specific repeat classes (uncomment to use)
awk 'tolower($11) ~ /satellite\/centr|simple_repeat|low_complexity/' "$REPEATMASKER_FILE" > "${AMP_OUTPUT_DIR}/chrY_SAT_features.out"
# Option 2: Use the entire RepeatMasker file as source (here we simply copy it)
# cp "$REPEATMASKER_FILE" "${AMP_OUTPUT_DIR}/chrY_SAT_features.out"

echo "Converting repeat features to BED format (0-based)..."
awk -v oc="$LOCAL_CHROMOSOME" 'BEGIN { OFS="\t" }
     NR > 2 && $6 ~ /^[0-9]+$/ && $7 ~ /^[0-9]+$/ {
         start = $6 - 1;
         end = $7;
         print oc, start, end, $11
     }' "$REPEATMASKER_FILE" > "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_SAT_features.bed"

echo "Applying hardmasking (Ns) to the masked FASTA..."
bedtools maskfasta -fi "$MASKED_Y_CHROMOSOME" \
  -bed "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_SAT_features.bed" \
  -fo "${AMP_OUTPUT_DIR}/chrY_hardmasked.fa" \
  -mc N

# Update the FASTA used downstream to the hardmasked version:
MASKED_Y_CHROMOSOME="${AMP_OUTPUT_DIR}/chrY_hardmasked.fa"

# --------------------------------------------------------------------
# STEP 3: Chromosome Filtering
# --------------------------------------------------------------------
echo "Filtering chromosome indexes..."
grep '^>' "$MASKED_Y_CHROMOSOME" | sed 's/^>//' | cut -d' ' -f1 > "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_valid_chroms.txt"
awk 'NR==FNR {chrom[$1]=1; next} $1 in chrom' \
  "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_valid_chroms.txt" \
  "$FAI_Y_CHROMOSOME" > "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_filtered.fai"

# --------------------------------------------------------------------
# STEP 4: Window Creation
# --------------------------------------------------------------------
echo "Generating sliding windows..."
bedtools makewindows -g "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_filtered.fai" \
  -w 5000 -s 2000 > "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_windows.bed"

# --------------------------------------------------------------------
# STEP 5: Extract FASTA Sequences for BLAST
# --------------------------------------------------------------------
echo "Extracting FASTA sequences for BLAST..."
bedtools getfasta -fi "$MASKED_Y_CHROMOSOME" \
  -bed "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_windows.bed" \
  -fo "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_windows.fa"

# --------------------------------------------------------------------
# STEP 6: Split FASTA into Individual BLAST Query Files
# --------------------------------------------------------------------
echo "Splitting FASTA into individual query files..."
awk '/^>/ {
   name = substr($1, 2);
   gsub(/[:]/, "_", name);  # Replace colon with underscore in the name
   filename = "'"${AMP_OUTPUT_DIR}"'/blast_queries/" name ".fa";
   print > filename;
   next
}
{ print > filename }' "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_windows.fa"

# --------------------------------------------------------------------
# STEP 7: Run BLASTn on Query Files
# --------------------------------------------------------------------
echo "Running BLASTn with GNU Parallel..."

# Clean previous results
mkdir -p "${AMP_OUTPUT_DIR}/blast_results" "${AMP_OUTPUT_DIR}/filtered_temp"

# Export variables needed for parallel processing
export MASKED_Y_CHROMOSOME AMP_OUTPUT_DIR LOCAL_CHROMOSOME

# Process files using parallel
find "${AMP_OUTPUT_DIR}/blast_queries" -name "*.fa" | \
parallel --progress --joblog "${AMP_OUTPUT_DIR}/blast_joblog.txt" \
    --eta --halt-on-error 1 \
    --env MASKED_Y_CHROMOSOME --env AMP_OUTPUT_DIR --env LOCAL_CHROMOSOME \
    -j 10 \
'
    INPUT="{}"
    WINDOW_ID=$(basename "$INPUT" .fa)
    BLAST_OUTPUT="${AMP_OUTPUT_DIR}/blast_results/${WINDOW_ID}_blast.tsv"
    FILTERED_TEMP="${AMP_OUTPUT_DIR}/filtered_temp/${WINDOW_ID}_filtered.tsv"

    # Skip empty input files
    if [ ! -s "$INPUT" ]; then
        echo "Skipping empty file: $INPUT"
        exit 0
    fi

    # Run BLASTn
    blastn -query "$INPUT" \
        -subject "$MASKED_Y_CHROMOSOME" \
        -perc_identity 50 \
        -outfmt "6 qseqid qstart qend sseqid sstart send" \
        -out "$BLAST_OUTPUT"

    # Process BLAST results
    if [ -s "$BLAST_OUTPUT" ]; then
        awk -v qseqid="$WINDOW_ID" '\''
            BEGIN {
                split(qseqid, parts, /-/)
                qchr = parts[1]
                qstart = parts[2]
                qend = parts[3]
            }
            # Filter out self-hits on same chromosome with overlap
            !($4 == qchr && ($5 <= qend && $6 >= qstart))
            { print }
        '\'' "$BLAST_OUTPUT" > "$FILTERED_TEMP"
    else
        touch "$FILTERED_TEMP"
    fi
'

# Combine filtered results
cat "${AMP_OUTPUT_DIR}/filtered_temp/"*.tsv > "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_filtered_blast.tsv"

echo "Parallel BLAST processing complete. Results combined in: ${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_filtered_blast.tsv"

# --------------------------------------------------------------------
# STEP 8: Check for satellite region file 
# --------------------------------------------------------------------
SATELLITE_FILE="${SAT_OUTPUT_DIR}/chrY_combined_final_filtered.bed"
echo "Checking for satellite file: ${SATELLITE_FILE}"
if [ ! -f "${SATELLITE_FILE}" ]; then
  echo "WARNING: Satellite file not found for ${LOCAL_CHROMOSOME}. Creating empty file."
  SATELLITE_FILE="${AMP_OUTPUT_DIR}/empty_satellites.bed"
  touch "${SATELLITE_FILE}"
fi

# --------------------------------------------------------------------
# STEP 9: Merge BLAST Hits with Coordinate Correction
# --------------------------------------------------------------------
echo "Merging BLAST hits..."
if [ -s "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_filtered_blast.tsv" ]; then
  awk '{ if ($5 > $6) print $4, $6, $5; else print $4, $5, $6 }' OFS='\t' "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_filtered_blast.tsv" | \
    bedtools sort | \
    bedtools merge -d 100000 | \
    awk '$3 - $2 >= 90000' > "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_blast_hits.bed"
  echo "Merged BLAST hits saved to: ${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_blast_hits.bed"
else
  echo "WARNING: No BLAST hits found; creating empty blast hits file."
  touch "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_blast_hits.bed"
fi

# --------------------------------------------------------------------
# STEP 10: Generate Final Ampliconic Regions
# --------------------------------------------------------------------
echo "Generating final ampliconic regions..."
cat <(grep -v '^#' "${PAL_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_palindromes.bed" | cut -f1-3) \
    "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_blast_hits.bed" | \
  bedtools sort | \
  bedtools merge | \
  bedtools subtract -a - -b "$PAR_FILE" | \
  bedtools subtract -a - -b "${SATELLITE_FILE}" \
    > "${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_amplicons_final.bed"
echo "Final ampliconic regions saved to: ${AMP_OUTPUT_DIR}/${LOCAL_CHROMOSOME}_amplicons_final.bed"

# --------------------------------------------------------------------
# STEP 11: Cleanup Temporary Files
# --------------------------------------------------------------------
echo "Cleaning up temporary files..."
# rm -rf "${AMP_OUTPUT_DIR}/blast_queries" "${AMP_OUTPUT_DIR}/blast_results"
# rm -f "${AMP_OUTPUT_DIR}/empty_satellites.bed"

echo "Ampliconic region analysis complete."
