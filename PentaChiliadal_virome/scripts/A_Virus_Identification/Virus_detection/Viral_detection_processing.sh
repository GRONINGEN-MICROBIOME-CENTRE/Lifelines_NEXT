#!/bin/bash
#SBATCH --job-name=Viral_detection_processing
#SBATCH --output=Viral_detection_processing.out
#SBATCH --mem=16gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

#Load modules
module purge; ml SeqKit; module list

completed_samples_dir=$1 # Directory with results of Viral detection step
MSP_contigs_dir=$2 # Directory with metaSPAdes contigs/scaffolds

ls $completed_samples_dir > list_samples_completed.txt

# For each completed sample, filter selected contigs and create a combined file
while read i
do
	sample=$i
	results_path=virome_discovery/processed_results
	mkdir -p $completed_samples_dir/$sample/$results_path
	awk 'NR>1 {print $0}' $completed_samples_dir/$sample/virome_discovery/DeepVirFinder/${sample}_metaspades_contigs.fa_gt10000bp_dvfpred.txt | awk '$3 >= 0.94' | cut -f1 | sed "s/^/${sample}_/" > $completed_samples_dir/$sample/$results_path/${sample}_DVF_contigs.txt
	awk 'NR>1 {print $1}' $completed_samples_dir/$sample/virome_discovery/VirSorter2/final-viral-score.tsv | sed -E 's/\|\|[a-z0-9_]+$//'| sort | uniq | sed "s/^/${sample}_/" > $completed_samples_dir/$sample/$results_path/${sample}_VS2_contigs.txt
	awk 'NR>1' $completed_samples_dir/$sample/virome_discovery/geNomad/${sample}_metaspades_contigs_virus_summary.tsv | awk -F '\t' '$2 >= 10000' | cut -f1 | sed "s/^/${sample}_/" > $completed_samples_dir/$sample/$results_path/${sample}_geNomad_contigs.txt
	
	cat $completed_samples_dir/$sample/$results_path/${sample}_VS2_contigs.txt $completed_samples_dir/$sample/$results_path/${sample}_DVF_contigs.txt $completed_samples_dir/$sample/$results_path/${sample}_geNomad_contigs.txt | sort | uniq > $completed_samples_dir/$sample/$results_path/${sample}_combined_contigs.txt
  
# Generate the FASTA file with the selected contigs using sekqit. Add sample name to contig_IDs
  sed "s/^>/>${sample}_/" "$MSP_contigs_dir/${sample}_metaspades_contigs.fa" > $MSP_contigs_dir/${sample}_metaspades_contigs_mod.fa
  seqkit grep -f $completed_samples_dir/$sample/$results_path/${sample}_combined_contigs.txt $MSP_contigs_dir/${sample}_metaspades_contigs_mod.fa > $completed_samples_dir/$sample/$results_path/${sample}_predicted_viral_contigs.fa

done	< list_samples_completed.txt

# Generate merged files
cat $(find $completed_samples_dir -type f -name "*_predicted_viral_contigs.fa") > all_predicted_viral_contigs.fa
cat $(find $completed_samples_dir -type f -name "*_VS2_contigs.txt") > all_VS2_viral_contigs.txt
cat $(find $completed_samples_dir -type f -name "*_DVF_contigs.txt") > all_DVF_viral_contigs.txt
cat $(find $completed_samples_dir -type f -name "*_geNomad_contigs.txt") > all_geNomad_viral_contigs.txt
cat $(find $completed_samples_dir -type f -name "*_combined_contigs.txt") > all_predicted_viral_contigs.txt

# Generate a merged file the origin (tool) of each contig in long format
sed 's/^/VS2 /' all_VS2_viral_contigs.txt > modified_VS2_viral_contigs.txt
sed 's/^/DVF /' all_DVF_viral_contigs.txt > modified_DVF_viral_contigs.txt
sed 's/^/geNomad /' all_geNomad_viral_contigs.txt > modified_geNomad_viral_contigs.txt

# Concatenate the modified files into a single merged file (long format)
cat modified_VS2_viral_contigs.txt modified_DVF_viral_contigs.txt \
modified_geNomad_viral_contigs.txt  > Viral_table_of_origin.txt

# Remove intermediate files
rm list_samples_completed.txt modified_*txt
