# #!/bin/bash

# # Activate Conda environment
# eval "$(conda shell.bash hook)"
# conda activate blast

# # Ensure seqkit is available
# export PATH=$(conda info --base)/envs/blast/bin:$PATH
# which seqkit || { echo "Error: seqkit not found in the activated environment"; exit 1; }

# # Define directories and parameters
# FASTA_DIR="/grphome/grp_lichenscapstone/genedata_fasta"
 CHUNK_DIR="/grphome/grp_lichenscapstone/genedata_fasta/chunks"
 DB_DIR="/grphome/grp_lichenscapstone/combined_db"
 RESULTS_DIR="/grphome/grp_lichenscapstone/blastx3-5"
# THRESHOLD=500000000  # Files larger than 500 MB will be split
# CHUNK_SIZE=500000      # Number of sequences per chunk

# # Ensure necessary directories exist
# mkdir -p "$CHUNK_DIR"
# mkdir -p "$RESULTS_DIR"

# # If not running as a "child", perform splitting and then submit separate SLURM jobs.
# echo "Starting FASTA file splitting..."
# for fasta in "$FASTA_DIR"/*.fasta; do
#     if [ -s "$fasta" ]; then
#         filesize=$(stat -c%s "$fasta") 
#         if [ $filesize -gt $THRESHOLD ]; then
#             echo "Splitting large file: $fasta (size: $filesize bytes)"
#             seqkit split -s $CHUNK_SIZE "$fasta" -O "$CHUNK_DIR"
#         else
#             echo "Copying small file: $fasta"
#             cp "$fasta" "$CHUNK_DIR"
#         fi
#     else
#         echo "Skipping empty file: $fasta"
#     fi
# done
# echo "FASTA file splitting completed."

 num_file=$(ls -1 $CHUNK_DIR | wc -l)
 sbatch --array 1-$num_file array.sh $CHUNK_DIR $DB_DIR $RESULTS_DIR
#echo $num_file
