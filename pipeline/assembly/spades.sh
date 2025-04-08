#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH -C 'intel'   # features syntax (use quotes): -C 'a&b&c&d'
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "kylemf-spades-$OUT_FOL"   # job name
#SBATCH --mail-user=kyle.matthew.frazier@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load gcc/12

GRP_DIR="../../"
IN_PATH="$GRP_DIR/raw_read_fasta/"
#OUT_FOL="run-1"
OUT_PATH="$GRP_DIR/spades-results/" #$OUT_FOL"
mkdir $OUT_PATH
INPUT_SUFFIX=".fasta"

for read in "$IN_PATH/"*.fasta
do 
    f_name="$(basename $read)"
    p_name="$(basename $read $INPUT_SUFFIX)"
    result_folder=$OUT_PATH/$p_name
    mkdir $result_folder

    $GRP_DIR/SPAdes-4.0.0-Linux/bin/spades.py -s $IN_PATH/$f_name -o $result_folder
    #Flag explanations:

    # using spades.py instead of metaspades.py b/c
    # "Currently metaSPAdes supports only a single short-read library which has to be paired-end"
    #https://ablab.github.io/spades/running.html#

    # -s
    #https://ablab.github.io/spades/running.html#examples:~:text=%2Ds%20%3Cfile_name%3E%20File%20with%20unpaired%20reads.

    # -o
    #https://ablab.github.io/spades/running.html#examples:~:text=%2Do%20%3Coutput_dir%3E%20Specify%20the%20output%20directory.%20Required%20option.

    # "Contigs/scaffolds names in SPAdes output FASTA files have the following format:"
    #https://ablab.github.io/spades/output.html#:~:text=Contigs/scaffolds%20names%20in%20SPAdes%20output%20FASTA%20files%20have%20the%20following%20format

done
