#!/bin/bash
#SBATCH --job-name=vOTU_clustering_processing
#SBATCH --output=vOTU_clustering_proc.out
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

clusters=$1 #CheckV aniclust.py clusters.tsv output file
viral_genomes_file=$2 #file with all the viral genomes used as input for the dereplication

#Load modules
ml purge; ml Python/3.10.8-GCCcore-12.2.0 seqtk; ml list

#Extract the names of the predicted viral sequences in the negative control sample
cat $viral_genomes_file | grep 'Neg_control' | cut -c2- > neg_control_IDs.txt

#Execute the python script that outputs Final_viral_sequences_after_dereplication.txt file with only the sequences that do not cluster with negative control sequences
python vOTUs_clustering_processing.py $clusters neg_control_IDs.txt

#Extract from the FASTA file with the contigs those that are not clustered with negative control sequences (after dereplication)
seqtk subseq $viral_genomes_file Final_viral_sequences_after_dereplication.txt > Dereplication_viral_sequences_no_neg.fa 
