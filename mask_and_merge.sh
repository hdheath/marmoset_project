#!/bin/bash
# --------------------------------------------------------------------
# Masking and Merging Script
# This script runs masking tools (RepeatMasker, WindowMasker, TRF),
# converts outputs to BED files, merges masking results, applies soft-masking,
# combines the softmasked FASTA files, and indexes them.
#
# Usage: ./mask_and_merge.sh path_to_Y_chromosome_FASTA path_to_X_chromosome_FASTA
# --------------------------------------------------------------------

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 path_to_Y_chromosome_FASTA path_to_X_chromosome_FASTA"
  exit 1
fi

# User-supplied FASTA file paths
Y_FASTA="$1"
X_FASTA="$2"

# Base directory is the current working directory
BASE_DIR="$(pwd)"

# Define directories for outputs and logs
PAR_ANNOT_DIR="$BASE_DIR/PAR_annotation"
RM_OUTPUT_DIR="$PAR_ANNOT_DIR/repeatmasker_output"
WM_OUTPUT_DIR="$PAR_ANNOT_DIR/windowmasker_output"
TRF_OUTPUT_DIR="$PAR_ANNOT_DIR/trf_output"
COMBINED_DIR="$PAR_ANNOT_DIR/combined_mask_fastas"
LOG_DIR="$PAR_ANNOT_DIR/logs"

# Create all necessary directories
mkdir -p "$PAR_ANNOT_DIR" "$RM_OUTPUT_DIR" "$WM_OUTPUT_DIR" "$TRF_OUTPUT_DIR" "$COMBINED_DIR" "$LOG_DIR"

# Progress tracking
total_steps=6
current_step=0
function update_progress {
  current_step=$((current_step+1))
  percent=$(( current_step * 100 / total_steps ))
  bar_length=50
  filled_length=$(( bar_length * current_step / total_steps ))
  bar=$(printf "%${filled_length}s")
  bar=${bar// /#}
  empty_length=$(( bar_length - filled_length ))
  empty=$(printf "%${empty_length}s")
  empty=${empty// /-}
  echo "Progress: [${bar}${empty}] ${percent}% complete"
}

#####################
# RepeatMasker Step #
#####################
# Run RepeatMasker on Y chromosome in background
RepeatMasker -s -pa 16 -e ncbi "$Y_FASTA" \
  -species "Callithrix jacchus" \
  -xsmall \
  -align \
  -dir "$RM_OUTPUT_DIR" > "$LOG_DIR/repeatmasker_Y.log" 2>&1 &
pid_rm_y=$!

# Run RepeatMasker on X chromosome in background
RepeatMasker -s -pa 16 -e ncbi "$X_FASTA" \
  -species "Callithrix jacchus" \
  -xsmall \
  -align \
  -dir "$RM_OUTPUT_DIR" > "$LOG_DIR/repeatmasker_X.log" 2>&1 &
pid_rm_x=$!

# Wait for both processes to complete
wait $pid_rm_y $pid_rm_x
echo "RepeatMasker completed for both chromosomes. Logs in $LOG_DIR."
update_progress

#####################
# WindowMasker Step #
#####################
# Process Y chromosome
windowmasker -in "$Y_FASTA" -mk_counts -out "$WM_OUTPUT_DIR/Y_wm_counts.txt" > "$LOG_DIR/windowmasker_Y_counts.log" 2>&1
windowmasker -in "$Y_FASTA" -ustat "$WM_OUTPUT_DIR/Y_wm_counts.txt" \
             -out "$WM_OUTPUT_DIR/Y_masked_genome.fa" \
             -outfmt fasta > "$LOG_DIR/windowmasker_Y_fasta.log" 2>&1 &
pid_wm_y_fasta=$!
windowmasker -in "$Y_FASTA" -ustat "$WM_OUTPUT_DIR/Y_wm_counts.txt" \
             -out "$WM_OUTPUT_DIR/Y_masked_regions.bed" \
             -outfmt interval > "$LOG_DIR/windowmasker_Y_bed.log" 2>&1 &
pid_wm_y_bed=$!

# Process X chromosome
windowmasker -in "$X_FASTA" -mk_counts -out "$WM_OUTPUT_DIR/X_wm_counts.txt" > "$LOG_DIR/windowmasker_X_counts.log" 2>&1
windowmasker -in "$X_FASTA" -ustat "$WM_OUTPUT_DIR/X_wm_counts.txt" \
             -out "$WM_OUTPUT_DIR/X_masked_genome.fa" \
             -outfmt fasta > "$LOG_DIR/windowmasker_X_fasta.log" 2>&1 &
pid_wm_x_fasta=$!
windowmasker -in "$X_FASTA" -ustat "$WM_OUTPUT_DIR/X_wm_counts.txt" \
             -out "$WM_OUTPUT_DIR/X_masked_regions.bed" \
             -outfmt interval > "$LOG_DIR/windowmasker_X_bed.log" 2>&1 &
pid_wm_x_bed=$!

# Wait for all WindowMasker processes to finish
wait $pid_wm_y_fasta $pid_wm_y_bed $pid_wm_x_fasta $pid_wm_x_bed
echo "WindowMasker completed for both chromosomes. Logs in $LOG_DIR."
update_progress

###########
# TRF Step#
###########
TRF_PARAMS="2 7 7 80 10 50 500"

# Run TRF on Y chromosome
echo "Running TRF on $Y_FASTA" | tee "$LOG_DIR/trf_Y.log"
trf "$Y_FASTA" $TRF_PARAMS -f -d -m > "$LOG_DIR/trf_Y_output.log" 2>&1 &
pid_trf_y=$!
wait $pid_trf_y
mv ./"$(basename "$Y_FASTA")".2.7.7.* "$TRF_OUTPUT_DIR/" 2>> "$LOG_DIR/trf_Y.log"
echo "TRF processing completed for $Y_FASTA" | tee -a "$LOG_DIR/trf_Y.log"
update_progress

# Run TRF on X chromosome
echo "Running TRF on $X_FASTA" | tee "$LOG_DIR/trf_X.log"
trf "$X_FASTA" $TRF_PARAMS -f -d -m > "$LOG_DIR/trf_X_output.log" 2>&1 &
pid_trf_x=$!
wait $pid_trf_x
mv ./"$(basename "$X_FASTA")".2.7.7.* "$TRF_OUTPUT_DIR/" 2>> "$LOG_DIR/trf_X.log"
echo "TRF processing completed for $X_FASTA" | tee -a "$LOG_DIR/trf_X.log"
update_progress

echo "TRF processing completed for both chromosomes."

# --------------------------------------------------------------------
# Create BED Files from Raw Masking Outputs
# --------------------------------------------------------------------
echo "Creating BED files from WindowMasker, RepeatMasker, and TRF outputs..."

# --- WindowMasker Conversion ---
awk 'BEGIN{OFS="\t"} /^>/ {chrom=substr($0,2)} /^[0-9]/ {split($0, a, " - "); print chrom, a[1], a[2]}' \
  "$WM_OUTPUT_DIR/X_masked_regions.bed" > "$WM_OUTPUT_DIR/X_masked_regions_converted.bed"
awk 'BEGIN{OFS="\t"} /^>/ {chrom=substr($0,2)} /^[0-9]/ {split($0, a, " - "); print chrom, a[1], a[2]}' \
  "$WM_OUTPUT_DIR/Y_masked_regions.bed" > "$WM_OUTPUT_DIR/Y_masked_regions_converted.bed"

# --- RepeatMasker Conversion ---
awk 'BEGIN{OFS="\t"} NR>1 && !/^SW/ && !/^score/ && NF>0 {print $5, $6, $7}' \
  "$RM_OUTPUT_DIR/$(basename "$X_FASTA").out" > "$RM_OUTPUT_DIR/chrX_repeatmasker.bed"
awk 'BEGIN{OFS="\t"} NR>1 && !/^SW/ && !/^score/ && NF>0 {print $5, $6, $7}' \
  "$RM_OUTPUT_DIR/$(basename "$Y_FASTA").out" > "$RM_OUTPUT_DIR/chrY_repeatmasker.bed"

# --- TRF Conversion ---
awk 'BEGIN{OFS="\t"} /^[0-9]/ {print "chrX_pri", $1, $2}' \
  "$TRF_OUTPUT_DIR/$(basename "$X_FASTA" .fa).fa.2.7.7.80.10.50.500.dat" > "$TRF_OUTPUT_DIR/chrX_trf.bed"
awk 'BEGIN{OFS="\t"} /^[0-9]/ {print "chrY_pri", $1, $2}' \
  "$TRF_OUTPUT_DIR/$(basename "$Y_FASTA" .fa).fa.2.7.7.80.10.50.500.dat" > "$TRF_OUTPUT_DIR/chrY_trf.bed"

echo "BED file conversion completed."

# --------------------------------------------------------------------
# Merge Masking Results, Apply Soft-Masking, Combine FASTAs, and Index
# --------------------------------------------------------------------
echo "Performing merge masking results, applying soft-masking, combining FASTAs, and indexing..."

###############
# Chromosome X
###############
cat "$WM_OUTPUT_DIR/X_masked_regions_converted.bed" "$RM_OUTPUT_DIR/chrX_repeatmasker.bed" "$TRF_OUTPUT_DIR/chrX_trf.bed" | \
  sort -k1,1 -k2,2n | bedtools merge -d 100 > "$COMBINED_DIR/merged_X.bed"
bedtools maskfasta -fi "$X_FASTA" -bed "$COMBINED_DIR/merged_X.bed" \
  -fo "$COMBINED_DIR/softmasked_chrX.fa" -soft

###############
# Chromosome Y
###############
cat "$WM_OUTPUT_DIR/Y_masked_regions_converted.bed" "$RM_OUTPUT_DIR/chrY_repeatmasker.bed" "$TRF_OUTPUT_DIR/chrY_trf.bed" | \
  sort -k1,1 -k2,2n | bedtools merge -d 100 > "$COMBINED_DIR/merged_Y.bed"
bedtools maskfasta -fi "$Y_FASTA" -bed "$COMBINED_DIR/merged_Y.bed" \
  -fo "$COMBINED_DIR/softmasked_chrY.fa" -soft

###############
# Combine FASTA Files 
###############
cat "$COMBINED_DIR/softmasked_chrX.fa" "$COMBINED_DIR/softmasked_chrY.fa" > "$COMBINED_DIR/combined_softmasked.fa"

###############
# Index FASTA Files
###############
samtools faidx "$COMBINED_DIR/softmasked_chrY.fa"
samtools faidx "$COMBINED_DIR/softmasked_chrX.fa"
samtools faidx "$COMBINED_DIR/combined_softmasked.fa"

echo "Combined and indexed softmasked FASTA files are available in $COMBINED_DIR."
