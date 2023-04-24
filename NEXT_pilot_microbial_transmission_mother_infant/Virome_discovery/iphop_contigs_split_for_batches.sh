#!/bin/bash

#SBATCH --job-name=iphop_contigs_split
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --output=iphop_contigs_split.out
#SBATCH --error=iphop_contigs_split.err
#SBATCH --mem=1GB

module purge; ml Anaconda3; conda activate iphop_env

iphop split \
	--input_file /scratch/umcg-nkuzub/iphop/iphop_pilot_filtered_shuffled_splitted/viral_noneg405_99_der95_decontaminated_filtered_shuffled.fasta \
	--split_dir /scratch/umcg-nkuzub/iphop/iphop_pilot_filtered_shuffled_splitted \
	--n_seq 500
