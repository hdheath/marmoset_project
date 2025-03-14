# Marmoset Chromosome Annotation Project

This repository provides scripts to annotate and characterize genomic regions on the marmoset Y and X chromosomes, including repeat masking and identification of pseudoautosomal (PAR), satellite, palindromic, ampliconic, and X-transposed regions.

## 1. Repeat Masking

**Script:** `mask_and_merge.sh`  
Runs masking tools (RepeatMasker, WindowMasker, TRF), converts outputs to BED format, merges results, applies soft-masking, combines FASTA files, and indexes the final masked sequences.

**Usage:**

```bash
./mask_and_merge.sh path_to_Y_chromosome_FASTA path_to_X_chromosome_FASTA

## 2. Annotate PAR Regions
Script: sedef_and_filter.sh
Identifies candidate pseudoautosomal regions (PAR) using Sedef on combined soft-masked FASTA files and filters results to produce refined PAR regions.

Note: Manual selection of a PAR region should be saved to:
PAR_annotation/PAR_region/PAR.bed

## 3. Annotate Large Satellite Regions
Script: satellite_regions.sh
Extracts satellite regions from RepeatMasker output, merges overlapping regions (within 1 kb), and subtracts previously annotated PAR regions.

Note: Ensure the existence of the file:
PAR_annotation/PAR_region/PAR.bed

## 4. Annotate Palindromic Regions
Script: palindrome_regions.sh
Detects palindromic sequences by performing self-alignments using LASTZ and applies filtering (Y chromosome only).

Note: You'll be prompted to provide the full path to palindrover.py.

Script for Both Chromosomes: x_and_y_palindrome_regions.sh
Detects palindrome sequences for both X and Y chromosomes.

## 5. Annotate Ampliconic Regions
Script: ampliconic_regions.sh
Annotates ampliconic sequences from the Y masked chromosome FASTA by:

Creating sliding windows across chromosome sequences.
Running BLASTn searches in parallel.
Merging BLAST results with palindromic annotations.
Excluding overlaps with PAR and satellite annotations.

## 6. Annotate X-Transposed Regions
Script: x_transposed_regions.sh
Identifies X-transposed regions by aligning chrY to chrX using LASTZ. Alignments are filtered for high identity (>94%) and length (>1 kb). Results exclude overlaps with PAR and ampliconic regions using bedtools.

## Plot Creation
Moddotplot for Chr Y
bash
Copy
python -m moddotplot static -f {path to Y chromosome}/chrY.fa -o {output path}
Create Bed Files for Plotting
Script: create_bed_for_plots.py

Plot Annotated Regions
Script: plot_chrs.py

Cactus Alignment
Create Distance Tree
Script: mashtree.sh
Creates a distance tree from input FASTA files.

Binarize the Distance Tree
Script: binarize_tree.r
Binarizes the distance tree for Cactus input.

Note: Must change static paths before use.
