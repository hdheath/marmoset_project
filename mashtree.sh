#!/bin/bash
#SBATCH --job-name=mashtree
#SBATCH --partition=medium
#SBATCH --mail-user=hdheath@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=16gb  # Adjust based on input size
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16  # Use multiple CPUs for faster execution
#SBATCH --output=mashtree_%j.log
#SBATCH --time=02:00:00  # Adjust as needed

# Load Conda environment (ensure you have Mashtree installed in this environment)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mashtree_env

# Define input and output directories
INPUT_DIR="/private/groups/cgl/hdheath/assembly/CGP_y"
OUTPUT_DIR="/private/groups/cgl/hdheath/assembly/CGP_y"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Move to the directory containing the Y chromosome FASTA files
cd "$INPUT_DIR" || { echo "Error: Input directory not found"; exit 1; }

# Run Mashtree with 16 CPU cores and save output to the specified directory
echo "Running Mashtree on Y chromosome FASTA files..."
mashtree --numcpus 16 *Y*.fa > "$OUTPUT_DIR/y_tree.nwk"

# Check if the tree was created successfully
if [ -s "$OUTPUT_DIR/y_tree.nwk" ]; then
    echo "Mashtree completed successfully. Tree saved to $OUTPUT_DIR/y_tree.nwk"
else
    echo "Error: y_tree.nwk file is empty or was not created."
    exit 1
fi

echo "Mashtree job finished."
