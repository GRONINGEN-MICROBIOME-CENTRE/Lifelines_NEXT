#!/bin/bash
#SBATCH --job-name=External_gene_calls 
#SBATCH --output=External_gene_calls.out
#SBATCH --mem=32gb
#SBATCH --time=3-0
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

prodigalgv_CSV=$1 #path to CSV file with prodigal-gv output
prodigalgv_proteins=$2 #path to FASTA file with prodigal-gv output

# Clean environment, load modules 
ml purge; ml R; ml list

# Run the R script
Rscript LLNEXT_DGR_Activity_Gene_Calls.R $prodigalgv_CSV $prodigalgv_proteins
