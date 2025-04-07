#!/bin/bash
#***** NOTE: run this using: sg group_name "sbatch thefilename"

#SBATCH --time=6-00:00:00   # walltime
#SBATCH --ntasks=64   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=20480M   # memory per CPU core
#SBATCH -J "REPLACE_JOB_NAME"   # job name
#SBATCH --mail-user=sample_user@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

if [ "$(id -gn)" != "group_name" ]; then
    echo '*!*!*' This job is not running as the intended group. If you want to run it as group_name, run sbatch as follows:  sg group_name '"'sbatch thefilename'"'
    exit 1
fi


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


module load blast-plus/2.14.1-szs4syo
blastn -query ../../home/USER/groups/group_name/file.fasta -db ../../apps/blast/databases/nt -out ../../home/USER/groups/group_name/file.out -evalue 1e-5