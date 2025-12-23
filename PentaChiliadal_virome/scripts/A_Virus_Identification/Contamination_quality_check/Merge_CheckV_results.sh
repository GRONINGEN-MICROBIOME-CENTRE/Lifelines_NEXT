#!/bin/bash
#SBATCH --job-name=merge_CheckV_results
#SBATCH --output=Merge_CheckV_final.out
#SBATCH --mem=4gb
#SBATCH --time=00:29:59
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

CheckV_output=$1 # CheckV output folder
mkdir -p CheckV_merged_results # directory to store results

echo -e '\n-------------------- MERGING CHECKV OUTPUT --------------------'

# Concatenate quality_summary.tsv, contamination.tsv, completeness.tsv and complete_genomes.tsv
# Important to use 'NR>1' to avoid concatenating headers
find $CheckV_output -name "quality_summary.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/quality_summary.tsv
find $CheckV_output -type f -name "quality_summary.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/quality_summary.tsv
find $CheckV_output -name "contamination.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/contamination_final.tsv
find $CheckV_output -type f -name "contamination.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/contamination_final.tsv
find $CheckV_output -name "completeness.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/completeness.tsv
find $CheckV_output -type f -name "completeness.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/completeness.tsv
find $CheckV_output -name "complete_genomes.tsv" -type f -exec head -1 {} \; -quit > CheckV_merged_results/complete_genomes.tsv
find $CheckV_output -type f -name "complete_genomes.tsv" -exec sed -e 1d '{}' \; >> CheckV_merged_results/complete_genomes.tsv


echo -e '\n-------------------- EXTRACTING VIRAL GENOMES OF INTEREST --------------------'

# Clean environment and load modules
module purge; ml seqtk; module list

cd CheckV_merged_results

# Get the names of contigs whose viral regions will be selected (completeness > 50% and host/viral gene ratio < 1) or discarded.
awk 'NR>1' quality_summary.tsv | awk '$8 != "Low-quality" && $8 != "Not-determined"' | awk '$6 > $7' | cut -f1 | sort > selected_CheckV_contigs.txt
awk 'NR>1' quality_summary.tsv | awk '$8 != "Low-quality" && $8 != "Not-determined"' | awk '$7 >= $6' | cut -f1 | sort > filtered_CheckV_contigs.txt
awk 'NR>1' quality_summary.tsv | awk '$8 == "Low-quality" || $8 == "Not-determined"' |  cut -f1 | sort >> filtered_CheckV_contigs.txt

# Repeat the previous process discarding also all possible chimeras detected by CheckV (contig >1.5x longer than expected genome length)
awk 'NR>1' quality_summary.tsv | awk '$8 != "Low-quality" && $8 != "Not-determined"' | awk '$6 > $7' | grep -v 'contig >1.5x longer than expected genome length' | cut -f1 | sort > selected_no_chimeras_CheckV_contigs.txt
awk 'NR>1' quality_summary.tsv | awk '$8 != "Low-quality" && $8 != "Not-determined"' | awk '$7 >= $6' | cut -f1 | sort > temp_filtered_chimeras_CheckV_contigs.txt
awk 'NR>1' quality_summary.tsv | awk '$8 == "Low-quality" || $8 == "Not-determined"'  | cut -f1 | sort >> temp_filtered_chimeras_CheckV_contigs.txt
awk 'NR>1' quality_summary.tsv | grep 'contig >1.5x longer than expected genome length' | cut -f1 | sort >> temp_filtered_chimeras_CheckV_contigs.txt
sort temp_filtered_chimeras_CheckV_contigs.txt | uniq > filtered_chimeras_CheckV_contigs.txt 
rm temp_filtered_chimeras_CheckV_contigs.txt

# Extract from the all_viruses.fna file the select viruses (after excluding potential chimeras)    
seqtk subseq \
    -l 80 \
    all_viruses.fna \
    selected_no_chimeras_CheckV_contigs.txt > CheckV_sequences_no_chimeras.fna  

# Set permissions
chmod 440 CheckV_sequences_no_chimeras.fna selected_CheckV_contigs.txt selected_no_chimeras_CheckV_contigs.txt
