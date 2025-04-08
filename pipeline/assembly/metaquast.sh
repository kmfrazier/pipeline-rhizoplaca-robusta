#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH -C 'intel'   # features syntax (use quotes): -C 'a&b&c&d'
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "kylemf-metaquast-$OUT_FOL"   # job name
#SBATCH --mail-user=kyle.matthew.frazier@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
grep -qxF '[[ -f ~/.bashrc ]] && . ~/.bashrc' ~/.bash_profile || echo '[[ -f ~/.bashrc ]] && . ~/.bashrc' >> ~/.bash_profile
source ~/.bash_profile
module load miniconda3 gcc/12
conda init
conda activate quast_env

GRP_DIR="../../"
IN_PATH="$GRP_DIR/spades-results"
NCBI_PATH="/apps/blast/databases/nr"  #-db ../../apps/blast/databases/nr
#OUT_FOL="test-6"
OUT_PATH="$GRP_DIR/metaquast-results/" #$OUT_FOL"
EXECUTABLE="/~/$USER/.conda/envs/quast_env/bin/metaquast.py"
mkdir $OUT_PATH

declare -a input_files
unset input_files
declare -a ref_files
unset ref_files

for folder in "$IN_PATH/"*
do 
    input_file="$folder/scaffolds.fasta"
    ref_file="$folder/contigs.fasta"
    if [ -f $input_file ] && [ -f $ref_file ]; then
        input_files+=($input_file)
        ref_files+=($ref_file)
    fi

done

python $EXECUTABLE "${input_files[@]}" -r "${ref_files[@]}" --blast-db $NCBI_PATH -L -o $OUT_PATH
#REFERENCE COMMAND: python metaquast.py contigs_1 contigs_2 ... -r reference_1,reference_2,reference_3,...

#Flag explanations:

# input_files and ref_files content:
#https://quast.sourceforge.net/docs/manual.html#faq_q6:~:text=user%20can%20compare%20results%20for%20real%20scaffolds%20and%20%22reconstructed%20contigs%22%20and%20find%20out%20whether%20scaffolding%20step%20was%20useful%20or%20not

# -r
#https://quast.sourceforge.net/docs/manual.html#:~:text=2.4%20Metagenomic%20assemblies

# --blast-db 
#https://quast.sourceforge.net/docs/manual.html#:~:text=reference%20search%20algorithm.-,%2D%2Dblast%2Ddb%20%3Cpath%3E,-Use%20custom%20BLAST

# --split-scaffolds
#don't use, b/c #https://quast.sourceforge.net/docs/manual.html#faq_q6:~:text=If%20you%20have%20both%20contigs.fasta%20and%20scaffolds.fasta%20it%20is%20better%20to%20specify%20both%20files%20to%20QUAST%20and%20don%27t%20set%20%2D%2Dsplit%2Dscaffolds%20option.
#otherwise, use b/c https://quast.sourceforge.net/docs/manual.html#:~:text=set%20to%204.-,Advanced%20options,-%3A

# -L
#https://quast.sourceforge.net/docs/manual.html#:~:text=2.4%2C%20IDBA%2DUD%22-,%2DL,-Take%20assembly%20names

# -o
#https://quast.sourceforge.net/docs/manual.html#:~:text=Options%3A-,%2Do,-%3Coutput_dir%3E
