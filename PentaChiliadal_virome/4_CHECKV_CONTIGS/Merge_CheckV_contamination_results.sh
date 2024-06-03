#!/bin/bash
#SBATCH --job-name=merge_CheckV_contamination_results
#SBATCH --output=merge_CheckV_contamination.out
#SBATCH --mem=4gb
#SBATCH --time=00:29:59
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

CheckV_output=$1 # CheckV output folder
mkdir -p CheckV_merged_results # directory to store results

# Concatenate viruses.fna files
cat $(find $CheckV_output/CheckV_results -type f -name "viruses.fna")  >> CheckV_merged_results/viruses.fna
# Concatenate proviruses.fna files
cat $(find $CheckV_output/CheckV_results -type f -name "proviruses.fna") >> CheckV_merged_results/proviruses.fna

# Generate a single FASTA file with all viral sequences detected by CheckV after trimming
sed 's, ,_,g' CheckV_merged_results/proviruses.fna > CheckV_merged_results/proviruses_mod.fna 
cat CheckV_merged_results/viruses.fna CheckV_merged_results/proviruses_mod.fna  > CheckV_merged_results/all_viruses.fna

# Concatenate contamination.tsv
find $CheckV_output/CheckV_results -name "contamination.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/contamination_initial.tsv
find $CheckV_output/CheckV_results -type f -name "contamination.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/contamination_initial.tsv

# Remove folder with split contigs and generate the new batch of split FASTA files
rm -r SPLIT_CONTIGS CheckV_results
sbatch Split_contigs.sh ${CheckV_output}/CheckV_merged_results/all_viruses.fna
