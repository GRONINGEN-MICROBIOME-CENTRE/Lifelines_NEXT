#!/bin/bash
#SBATCH --job-name=HMMER_results_processing
#SBATCH --output=HMMER_results_processing.out
#SBATCH --mem=16gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

metadata_file=$1 # Path to metadata file of protein families from public DBs
hmmer_result_PHROGs=$2 # Path the hmmsearch results file (processed)
output_dir=$3 # Path to output directory

mkdir -p $output_dir

# Load necessary modules
module load R

# Run the R script 
Rscript HMMER_result_processing.R $metadata_file $hmmer_result_PHROGs $output_dir

