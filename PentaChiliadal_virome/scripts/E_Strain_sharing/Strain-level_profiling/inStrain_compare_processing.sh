#!/bin/bash
#SBATCH --job-name=inStrain_compare_processing
#SBATCH --output=inStrain_compare_processing.out
#SBATCH --mem=8gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

IS_compare=$1 # Path to directory with the inStrain compare results

# Create or empty the output file
output_file="inStrain_filtered_075_comparisonsTable.tsv"
> $output_file

# Variable to track if the header has been added
header_added=false

# Find all *_comparisonsTable.tsv files and process them
find "$IS_compare" -type f -name "*_comparisonsTable.tsv" | while read -r file; do
    if [ "$header_added" = false ]; then
        # Add header from the first file to the output
        head -n 1 "$file" >> "$output_file"
        header_added=true
    fi
    # Filter rows where the 6th column (percentage of the genome compared) is >= 0.75 and append to the output file
    awk 'NR > 1 && $6 >= 0.75' "$file" >> "$output_file"
done

echo "Filtered rows with headers have been concatenated into $output_file"


