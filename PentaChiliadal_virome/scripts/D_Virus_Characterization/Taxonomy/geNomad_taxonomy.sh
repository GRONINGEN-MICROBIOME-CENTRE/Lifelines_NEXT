#!/bin/bash
#SBATCH --job-name=geNomad_tax
#SBATCH --output=geNomad_tax.out
#SBATCH --mem=60gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

virus_file=$1 #path to FASTA file with viral sequences
output_dir=$2 #path to output directory

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /home2/p304845/Conda_envs/geNomad_conda/; conda list 

# Run geNomad
genomad annotate \
        --cleanup \
        $virus_file \
        $output_dir \
        /home2/p304845/Conda_envs/geNomad_conda/genomad_db/

conda deactivate
