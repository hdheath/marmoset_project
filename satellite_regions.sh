#!/bin/bash
# --------------------------------------------------------------------
# Satellite Track Processing (Optimized for Actual Repeat Classes)
# --------------------------------------------------------------------

BASE_DIR="$(pwd)"
PAR_FILE="$BASE_DIR/PAR_annotation/PAR_region/PAR.bed"
REPEATMASKER_OUTPUT="$BASE_DIR/PAR_annotation/repeatmasker_output"
non_PAR_ANNOT_DIR="$BASE_DIR/non_PAR_annotation"
SAT_OUTPUT_DIR="$non_PAR_ANNOT_DIR/satellite_output"

mkdir -p "$non_PAR_ANNOT_DIR" "$SAT_OUTPUT_DIR"

process_chromosome() {
  local chr=$1

  echo "Processing $chr satellite regions..."

  # Get chromosome name (column 5)
  local actual_chr=$(awk 'NR > 2 && $5 != "" {print $5; exit}' "${REPEATMASKER_OUTPUT}/${chr}.fa.out")
  [ -z "$actual_chr" ] && { echo "Error: Chromosome name not found"; return 1; }

  # Step 1: Extract all repeat features (column 11)
  awk 'NR > 2 && tolower($11) ~ /satellite\/centr|simple_repeat|low_complexity|sine|alu/' \
    "${REPEATMASKER_OUTPUT}/${chr}.fa.out" > "${SAT_OUTPUT_DIR}/${chr}_all_features.out"

  [ ! -s "${SAT_OUTPUT_DIR}/${chr}_all_features.out" ] && { echo "Warning: No features found"; return; }

  # Step 2: Split into categories
  # Alpha satellites (Satellite/centr)
  awk 'tolower($11) ~ /satellite\/centr/ {print $0}' \
    "${SAT_OUTPUT_DIR}/${chr}_all_features.out" > "${SAT_OUTPUT_DIR}/${chr}_alpha_satellites.out"

  # Non-alpha satellites (Simple_repeat)
  awk 'tolower($11) ~ /simple_repeat|low_complexity|sine|alu/ {print $0}' \
    "${SAT_OUTPUT_DIR}/${chr}_all_features.out" > "${SAT_OUTPUT_DIR}/${chr}_satellites_no_alpha.out"

  # Process categories
  for category in satellites_no_alpha alpha_satellites combined; do
    echo "Processing $category..."
    
    case $category in
      "combined") 
        input_file="${SAT_OUTPUT_DIR}/${chr}_all_features.out"
        prefix="${SAT_OUTPUT_DIR}/${chr}_combined" ;;
      *)
        input_file="${SAT_OUTPUT_DIR}/${chr}_${category}.out"
        prefix="${SAT_OUTPUT_DIR}/${chr}_${category}" ;;
    esac

    [ ! -s "$input_file" ] && { echo "Warning: No $category regions"; continue; }

    # Convert to BED (0-based)
    awk -v oc="$actual_chr" 'BEGIN {OFS="\t"} 
        { 
          start = $6 - 1  # Convert to 0-based
          end = $7
          print oc, start, end, $11
        }' "$input_file" > "${prefix}.bed"

    #################
    # Note : the original paper uses a merge of 1000, and filter size of 250,000. 
    # However, none matched the criteria here. So we increased the merge gap and reduced filter size
    #################

    # Merge and subtract PARs
    bedtools merge -i "${prefix}.bed" -d 5000 > "${prefix}_merged.bed"
    bedtools subtract -a "${prefix}_merged.bed" -b "$PAR_FILE" > "${prefix}_final_noPAR.bed"
    
    # Size filters
    if [[ "$category" == "alpha_satellites" ]]; then
      # No size filter for alpha
      cp "${prefix}_final_noPAR.bed" "${prefix}_final_filtered.bed"
    else
      # 250kb filter for non-alpha
      awk -v OFS="\t" '$3 - $2 >= 250000' "${prefix}_final_noPAR.bed" > "${prefix}_final_filtered.bed"
    fi

    [ -s "${prefix}_final_filtered.bed" ] || echo "Warning: Empty final $category"
  done
}

# Execute processing
process_chromosome "chrY"
process_chromosome "chrX"

echo "Processing complete. Check outputs in: $SAT_OUTPUT_DIR"
