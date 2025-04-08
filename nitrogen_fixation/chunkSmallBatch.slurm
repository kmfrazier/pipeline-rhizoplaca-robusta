#!/bin/bash
#SBATCH --job-name=blastx_combined
#SBATCH --output=/grphome/grp_lichenscapstone/blastx_analysis/blastx_%A_%a.out
#SBATCH --error=/grphome/grp_lichenscapstone/blastx_analysis/blastx_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# This combined script:
# 1. Splits large FASTA files into chunks (or copies small files) into a chunks directory.
# 2. Submits a SLURM job array (with "child" argument) to run BLASTX on each chunk file.
#
# The BLAST database is local at /grphome/grp_lichenscapstone/ncbi_db/nr

# Activate Conda environment (adjust as needed)
eval "$(conda shell.bash hook)"
conda activate blast

# Define directories and parameters
FASTA_DIR="/grphome/grp_lichenscapstone/genedata_fasta"
CHUNK_DIR="/grphome/grp_lichenscapstone/genedata_fasta/chunks"
DB_DIR="/grphome/grp_lichenscapstone/ncbi_db/nr"
RESULTS_DIR="/grphome/grp_lichenscapstone/blastx_analysis"
THRESHOLD=500000000   # Files larger than 500 MB will be split
CHUNK_SIZE=1000       # Number of sequences per chunk

# Ensure necessary directories exist
mkdir -p "$CHUNK_DIR"
mkdir -p "$RESULTS_DIR"

# If not running as a "child", perform splitting and then submit the job array.
if [ "$1" != "child" ]; then
    echo "Starting FASTA file splitting..."
    for fasta in "$FASTA_DIR"/*.fasta; do
        if [ -s "$fasta" ]; then
            filesize=$(stat -c%s "$fasta")
            if [ $filesize -gt $THRESHOLD ]; then
                echo "Splitting large file: $fasta (size: $filesize bytes)"
                # Split into chunks of CHUNK_SIZE sequences
                seqkit split -s $CHUNK_SIZE "$fasta" -O "$CHUNK_DIR"
            else
                echo "Copying small file: $fasta"
                cp "$fasta" "$CHUNK_DIR"
            fi
        else
            echo "Skipping empty file: $fasta"
        fi
    done
    echo "FASTA file splitting completed."

    # Count the number of chunk files
    num_chunks=$(ls -1 "$CHUNK_DIR"/*.fasta 2>/dev/null | wc -l)
    if [ $num_chunks -eq 0 ]; then
        echo "No chunk files found. Exiting."
        exit 1
    fi
    echo "Found $num_chunks chunk files."

    # Submit a SLURM job array to process each chunk.
    echo "Submitting job array for BLASTX on chunk files."
    sbatch --array=0-$(($num_chunks - 1)) $0 child
    exit 0
fi

# ----- Child Phase: Process a single chunk file with BLASTX -----

# Build an array of chunk files
chunk_files=("$CHUNK_DIR"/*.fasta)
# Select the file corresponding to the current SLURM_ARRAY_TASK_ID
fasta="${chunk_files[$SLURM_ARRAY_TASK_ID]}"

# Check that the file exists and is not empty
if [ ! -s "$fasta" ]; then
    echo "File $fasta does not exist or is empty. Skipping..."
    exit 0
fi

# Derive a base name for output files
base_name=$(basename "$fasta" .fasta)
output_file="$RESULTS_DIR/${base_name}_blastx.txt"
error_file="$RESULTS_DIR/${base_name}_blastx.err"

echo "Running BLASTX on chunk: $fasta"
blastx -query "$fasta" \
       -db "$DB_DIR" \
       -out "$output_file" \
       -num_threads 4 \
       -evalue 1e-5 \
       -max_target_seqs 10 \
       -outfmt 6 \
       2> "$error_file"
echo "Finished processing chunk: $fasta"

