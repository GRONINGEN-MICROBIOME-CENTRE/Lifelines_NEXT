#!/bin/bash
#SBATCH --job-name=DB_origin_process_duplicates
#SBATCH --output=DB_process_duplicates.out
#SBATCH --mem=40gb
#SBATCH --time=00:19:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

viral_fasta_dir=$1 # directory with all the FASTA files to be used as input for dereplication

# Clean environment, load modules
module purge; ml R; module list

# Concatenate all the viral sequences to be used for dereplication
cat $viral_fasta_dir/*fna > Dereplication_input_sequences.fa
cat $viral_fasta_dir/*fa  >> Dereplication_input_sequences.fa
cat $viral_fasta_dir/*fasta >> Dereplication_input_sequences.fa

# Execute the R script
Rscript DB_origin_check_duplicated_sequences.R Dereplication_input_sequences.fa

# Set permissions
chmod 440 Dereplication_combined_sequences_nodup_renamed.fa Dereplication_input_sequences_nodup_DB_origin.txt
