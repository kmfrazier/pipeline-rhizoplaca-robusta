#!/bin/bash
#***** NOTE: run this using: sg grp_lichenscapstone "sbatch thefilename"

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=5120M   # memory per CPU core
#SBATCH -J "blastgenes"   # job name
#SBATCH --mail-user=ncc32@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

if [ "$(id -gn)" != "grp_lichenscapstone" ]; then
    echo '*!*!*' This job is not running as the intended group. If you want to run it as grp_lichenscapstone, run sbatch as follows:  sg grp_lichenscapstone '"'sbatch thefilename'"'
    exit 1
fi


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


module load blast-plus/2.14.1-szs4syo
tblastn -query photosystemI_fasta_file_path -db sl20079_comb.fasta_database_path -out outfile_path_and_name