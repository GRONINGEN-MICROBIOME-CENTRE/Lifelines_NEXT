#!/bin/bash
#SBATCH --job-name=merge_HMMER_results
#SBATCH --output=merge_HMMER.out
#SBATCH --mem=40gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

HMMER_output=$1 # HMMER output folder

# Concatenate HMMER_result and HMMER_table files
# Table: header 3 lines and end description 10 lines
# Result: header 12 lines and end description 1 line

# Concatenate Tables (keeping only the first header and excluding description at the end)
head -n 3 $(find $HMMER_output -type f -name "*HMMER_table" | head -n 1) > $HMMER_output/HMMER_table.tsv

for file in $HMMER_output/*HMMER_table; do
    tail -n +4 "$file" | head -n -10 >> $HMMER_output/HMMER_table.tsv
done

# Concatenate Results (keeping only the first header and excluding description at the end)
head -n 12 $(find $HMMER_output -type f -name "*HMMER_result" | head -n 1) > $HMMER_output/HMMER_result.tsv
for file in $HMMER_output/*HMMER_result; do
    tail -n +13 "$file" | head -n -1 >> $HMMER_output/HMMER_result.tsv
done

# Remove individual files
rm $HMMER_output/*table $HMMER_output/*result

# Retrive the relevant information from the HMMER_table.tsv 
# “target name” and “query name”, E-value, Score and the protein description columns
# When "accession" column is available, retrieve instead if "query name" - Useful for Pfam) 
awk 'BEGIN {OFS="\t"; print "Protein_ID", "Match", "E_value", "Score", "Description"} NR>3 && $4 != "-" { print $1, $4, $5, $6, $NF } NR>3 && $4 == "-" { print $1, $3, $5, $6, $NF }' $HMMER_output/HMMER_table.tsv | grep -E -v "^#" > $HMMER_output/HMMER_processed_table.tsv
