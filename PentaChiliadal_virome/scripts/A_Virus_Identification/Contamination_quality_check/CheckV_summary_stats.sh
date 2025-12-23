#!/bin/bash
#SBATCH --job-name=CheckV_summary_stats
#SBATCH --output=CheckV_summary.out
#SBATCH --mem=2gb
#SBATCH --time=00:19:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

CheckV_dir=$1 #directory with CheckV results
cd $CheckV_dir
mkdir -p CheckV_summary_stats_files
exec > CheckV_summary_stats_files/CheckV_summary_info.txt
echo -e "\n" 

echo -e "################################### CHECKV SUMMARY STATS ###################################\n"
# 1. Get files and summary stats of CheckV results
n_provirus=$(cat proviruses.fna | grep ">" | wc -l)
n_virus=$(cat viruses.fna | grep ">" | wc -l)
echo "The number of viral sequences detected is: $n_virus (viruses.fna file)"
echo -e "The number of proviral sequences detected is: $n_provirus (proviruses.fna file)\n"
echo "NOTE: The total number of proviruses should be extracted fromm both geNomad and CheckV outputs"
sum=$(( $n_provirus + $n_virus ))
echo -e "The total number of viral/proviral sequences identified by CheckV is: $sum\n"

# Get the summary stats of the total number of contigs to be kept (completeness > 50% and n viral genes > host genes) or discarded (with and without potential chimeras).
n_selected_contigs=$(cat selected_CheckV_contigs.txt| wc -l)
n_filtered_contigs=$(cat filtered_CheckV_contigs.txt| wc -l)
n_selected_contigs_no_chimera=$(cat selected_no_chimeras_CheckV_contigs.txt | wc -l)
n_filtered_contigs_chimera=$(cat filtered_chimeras_CheckV_contigs.txt | wc -l)
awk 'NR>1' quality_summary.tsv | grep "1.5x longer" | cut -f1 | sort > potential_chimeras.txt
n_chimeras=$(awk 'NR>1' quality_summary.tsv | grep "1.5x longer" | wc -l)
awk 'NR>1' quality_summary.tsv | grep "Not-determined" | cut -f1 | sort > not_determined_comp_CheckV_contigs.txt
n_not_det_contigs=$(cat not_determined_comp_CheckV_contigs.txt| wc -l)
echo "The total number of contigs selected (completeness >50%; viral/host gene ratio >1) is: $n_selected_contigs"
echo "The total number of contigs filtered (completeness =<50%; Not-determined; viral/host gene ratio =<1) is: $n_filtered_contigs"
echo "The total number of non-chimeric contigs selected (completeness >50%; viral/host gene ratio >1) is: $n_selected_contigs_no_chimera"
echo -e "The total number of contigs filtered (completeness =<50%; Not-determined; viral/host gene ratio =<1; potential chimeras) is: $n_filtered_contigs_chimera\n"

echo -e "The total number of contigs detected as potential chimeras is: $n_chimeras"
echo -e "The total number of contigs with not-determined completeness is: $n_not_det_contigs\n"

# 2. Get the file and summary stats of viral contigs selected but with 0 viral genes
grep -Fwf selected_CheckV_contigs.txt quality_summary.tsv | grep "no viral genes detected" | cut -f1 | sort > selected_nonviralgenes_CheckV_contigs_no_chimera.txt
n_selected_contigs_no_viral_genes=$(cat selected_nonviralgenes_CheckV_contigs_no_chimera.txt| wc -l)
echo -e "The total number of non-chimeric contigs selected but with 0 viral genes is: $n_selected_contigs_no_viral_genes\n"

echo -e "The final sequences/contigs are available in CheckV_sequences.fna and CheckV_sequences_no_chimeras.fna files\n"

# Move all the generated summary stats files to a folder
mv potential_chimeras.txt *CheckV*.txt CheckV_summary_stats_files
echo -e "########################################### END ###########################################\n"
