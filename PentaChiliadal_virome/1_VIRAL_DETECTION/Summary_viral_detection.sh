#!/bin/bash
#SBATCH --job-name=Summary_viral_detection
#SBATCH --output=Summary_viral_detection.out
#SBATCH --mem=4gb
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --export=NONE
#SBATCH --partition=regular

results_folder=$1 #Folder with Viral detection results per sample

ls $results_folder > list_samples.txt
while read i
do
	sample=$(echo $i)
	echo ${sample} | tr '\n' '\t'

# Check VS2 results
	if [ -f "${results_folder}/${sample}/virome_discovery/VirSorter2/final-viral-combined.fa" ]; then # check if VirSorter2 output file is present
		echo "yes" | tr '\n' '\t'
		awk 'NR>1 {print $1}' ${results_folder}/${sample}/virome_discovery/VirSorter2/final-viral-score.tsv | sed -E 's/\|\|[a-z0-9_]+$//' | sort | uniq | wc -l | tr '\n' '\t'
	else
		echo -e 'no\t-' | tr '\n' '\t'
	fi

# Check DVF results
	if [[ -f "${results_folder}/${sample}/virome_discovery/DeepVirFinder/${sample}_metaspades_contigs.fa_gt10000bp_dvfpred.txt" ]]; then #check if DeepVirFinder output file is present
		echo "yes" | tr '\n' '\t'
		echo $(cat ${results_folder}/${sample}/virome_discovery/DeepVirFinder/${sample}_metaspades_contigs.fa_gt10000bp_dvfpred.txt | tail -n +2 | awk '$3 >= 0.95 && $4 < 0.01'| wc -l) | tr '\n' '\t'
    	else
      		echo -e 'no\t-' | tr '\n' '\t'
    	fi

# Check geNomad result
	if [[ -f "${results_folder}/${sample}/virome_discovery/geNomad/${sample}_metaspades_contigs_virus_summary.tsv" ]]; then #check if geNomad output file is present
		echo "yes" | tr '\n' '\t'
		awk 'NR>1' ${results_folder}/${sample}/virome_discovery/geNomad/${sample}_metaspades_contigs_virus_summary.tsv | awk -F '\t' '$2 >= 10000 && $8 < 0.05' | wc -l 
    	else
      		echo -e 'no' 
    	fi

done < list_samples.txt > Viral_detection_results_summary.txt

# Remove intermediate files
rm list_samples.txt

# Add header to results file
sed -i $'1 i\\\nSample\tVS2_file\tVS2_contigs\tDVF_file\tDVF_Contigs p<0.01 & score>=0.95\tgeNomad_file\tgeNomad_Contigs FDR<0.05' Viral_detection_results_summary.txt
