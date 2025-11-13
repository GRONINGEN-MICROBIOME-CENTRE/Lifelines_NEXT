#!/bin/bash
#SBATCH --job-name=AMGP_KFS_2
#SBATCH --output=./out/AMGP/AMGP_KFS_2_%A.out
#SBATCH --error=./err/AMGP/AMGP_KFS_2_%A.err
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

KoFamScan_output=$1 # KoFamScan output folder

# Generate empty files to store output
ANNOTATIONS_MERGED="${KoFamScan_output}/ko-annotations.tsv"
FILTERED_MERGED="${KoFamScan_output}/ko-annotations-filtered.tsv"
> $ANNOTATIONS_MERGED
> $FILTERED_MERGED

# Add header from the first file (lines 1–2 for *_ko-annotations.tsv)
first_file=$(ls ${KoFamScan_output}/*_ko-annotations.tsv | head -n 1)
head -n 2 "$first_file" >> $ANNOTATIONS_MERGED

# Append the rest of all *_ko-annotations.tsv files (skip first 2 lines)
for file in ${KoFamScan_output}/*_ko-annotations.tsv; do
    echo "Processing $file"
    tail -n +3 $file >> $ANNOTATIONS_MERGED
done

# Add header from the first file (line 1 for *_ko-annotations-filtered.tsv)
first_filtered=$(ls ${KoFamScan_output}/*_ko-annotations-filtered.tsv | head -n 1)
head -n 1 $first_filtered >> $FILTERED_MERGED

# Append the rest of all *_ko-annotations-filtered.tsv files (skip 1 line)
for file in ${KoFamScan_output}/*_ko-annotations-filtered.tsv; do
    echo "Processing $file"
    tail -n +2 $file >> $FILTERED_MERGED
done

# Optional: Remove parts
# rm ${KoFamScan_output}/*part*
