#!/bin/bash
#SBATCH --job-name=merge_KOFamScan_results
#SBATCH --output=merge_KOFamScan.out
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

KoFamScan_output=$1 # KoFamScan output folder

# Generate empty file to store output
ANNOTATIONS_MERGED="${KoFamScan_output}/ko-annotations.tsv"
> $ANNOTATIONS_MERGED

# Add header from the first file (lines 1–2 for *_ko-annotations.tsv)
first_file=$(ls ${KoFamScan_output}/*_ko-annotations.tsv | head -n 1)
head -n 2 "$first_file" >> $ANNOTATIONS_MERGED

# Append the rest of all *_ko-annotations.tsv files (skip first 2 lines)
for file in ${KoFamScan_output}/*_ko-annotations.tsv; do
    echo "Processing $file"
    tail -n +3 $file >> $ANNOTATIONS_MERGED
done

# Remove intermediate files
rm ${KoFamScan_output}/*part*
