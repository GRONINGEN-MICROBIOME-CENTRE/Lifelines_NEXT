#!/bin/bash
#SBATCH --job-name=BACPHLIP_present_vOTUs
#SBATCH --output=BACPHLIP_present_vOTUs.out
#SBATCH --mem=30gb
#SBATCH --time=30:59:00
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

vOTU_seqs=$1 #FASTA file of representative sequences of vOTUs present in the samples

# Clean environment, load modules and activate conda environment
source /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/miniconda3/bin/activate \
    /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/BACPHLIP; conda list

# Run BACPHLIP
bacphlip -i $vOTU_seqs --multi_fasta
    
conda deactivate

