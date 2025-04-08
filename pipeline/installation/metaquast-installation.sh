#!/bin/bash
# Load Miniconda module
module load miniconda3

# Initialize Conda for bash
conda init bash

# Ensure ~/.bash_profile sources ~/.bashrc
grep -qxF '[[ -f ~/.bashrc ]] && . ~/.bashrc' ~/.bash_profile || echo '[[ -f ~/.bashrc ]] && . ~/.bashrc' >> ~/.bash_profile
source ~/.bash_profile

# Activate the base Conda environment
conda activate

# Create and activate the quast_env Conda environment
conda create --name quast_env -y
conda activate quast_env

# Install QUAST and dependencies using Mamba
mamba install -c bioconda quast -y

echo "QUAST environment setup is complete. Use 'conda activate quast_env' to activate the environment."
