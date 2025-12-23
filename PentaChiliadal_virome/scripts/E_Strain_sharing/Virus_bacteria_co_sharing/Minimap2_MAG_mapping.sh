#!/bin/bash
#SBATCH --job-name=Minimap2_MAG_mapping
#SBATCH --output=Minimap2_MAG_mapping.out
#SBATCH --mem=60gb
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

MAG_file=$1   # path to FASTA file with concatenated MAGs
PHAGE_file=$2 # path to FASTA file with your predicted phage genomes

mkdir -p MINIMAP_RESULTS

# Load modules
module purge; ml minimap2/2.26-GCCcore-12.3.0; ml list  

# Create minimap2 index of MAGs
minimap2 -d MINIMAP_RESULTS/MAGs.mmi $MAG_file

# Map phages to MAGs, with output in PAF
minimap2 -x asm5 -t ${SLURM_CPUS_PER_TASK} \
  MINIMAP_RESULTS/MAGs.mmi $PHAGE_file \
  > MINIMAP_RESULTS/phages_vs_MAGs.paf

