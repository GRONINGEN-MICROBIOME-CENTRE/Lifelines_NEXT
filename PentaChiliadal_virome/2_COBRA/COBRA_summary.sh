#!/bin/bash
#SBATCH --job-name=COBRA_summary
#SBATCH --output=COBRA_summary.out
#SBATCH --mem=8gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

COBRA_dir=$1 # directory with the COBRA results

# Generate final FASTA file after COBRA extension (with extended and non-extended contigs)
cat $(find $COBRA_dir -name "*fna" -type f) > $COBRA_dir/all_predicted_viral_sequences_extended.fa

# Create a file with the summary of the results
module purge; ml SeqKit; module list
exec > ${COBRA_dir}/COBRA_summary_info.txt
echo -e "\n" 

echo -e "################################### COBRA SUMMARY STATS ###################################\n"

# Get summary stats of COBRA results
n_extended_circular=$(cat $(find $COBRA_dir -type f -name "*extended_circular_unique.fasta.summary.txt") | grep -v SeqID | wc -l)
n_extended_partial=$(cat $(find $COBRA_dir -type f -name "*extended_partial_unique.fasta.summary.txt") | grep -v SeqID | wc -l)
n_extended_failed=$(cat $(find $COBRA_dir -type f -name "*extended_failed.fasta.summary.txt") | grep -v SeqID | wc -l)
n_extended_failed_orphan_end=$(cat $(find $COBRA_dir -type f -name "*end.fasta.summary.txt") | grep -v SeqID | wc -l)
n_no_extended_selfcircular=$(cat $(find $COBRA_dir -type f -name "*self_circular.fasta.summary.txt") | grep -v SeqID | wc -l)


echo "The number of newly extended circular genomes: $n_extended_circular"
echo "The number of newly partially extended contigs: $n_extended_partial"
echo "The number of contigs that were NOT extended (total): $(( $n_extended_failed + $n_extended_failed_orphan_end + $n_no_extended_selfcircular ))"
echo "The number of contigs that failed to be extended (COBRA rules): $n_extended_failed" 
echo "The number of contigs that failed to be extended (orphan end): $n_extended_failed_orphan_end" 
echo -e "The number of contigs that failed to be extended (already circular): $n_no_extended_selfcircular\n" 


n_viral_contigs=$(cat ${COBRA_dir}/all_predicted_viral_sequences_extended.fa | grep ">"| wc -l)
summary_viral_contigs=$(seqkit stats -a ${COBRA_dir}/all_predicted_viral_sequences_extended.fa)
echo "The total number of predicted viral contigs after COBRA re-assembly is: $n_viral_contigs"
echo "The summary stats of the predicted viral contigs after COBRA re-assembly:"
echo -e "$summary_viral_contigs\n"

echo -e "all_predicted_viral_sequences_extended.fa contains the sequences of all contigs after COBRA re-assembly. It can be used as input for CheckV.\n"
echo -e "###################################################### END ######################################################\n"

