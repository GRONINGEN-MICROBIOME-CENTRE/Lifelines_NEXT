#!/bin/bash
#SBATCH --job-name=HMMER_results_processing
#SBATCH --output=HMMER_results_processing.out
#SBATCH --mem=32gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

metadata_file=$1 # Path to metadata file of protein families from public DBs
hmmer_result_NCBI_Pfam=$2 # Path to the hmmsearch results file (processed)
hmmer_result_PHROGs=$3 # Path the hmmsearch results file (processed)
output_dir=$4 # Path to output directory

mkdir -p $output_dir

# Load necessary modules
module load R

# Run the R script 
Rscript HMMER_result_processing.R $metadata_file $hmmer_result_NCBI_Pfam $hmmer_result_PHROGs $output_dir

