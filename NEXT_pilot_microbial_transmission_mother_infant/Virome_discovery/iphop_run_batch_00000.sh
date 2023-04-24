#!/bin/bash

#SBATCH --job-name=iphop_batch_00000
#SBATCH --time=17:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=iphop_batch_00000.out
#SBATCH --error=iphop_batch_00000.err
#SBATCH --mem=55GB

module purge; ml Anaconda3; conda activate iphop_env

mkdir /scratch/umcg-nkuzub/iphop/RESULTS/batch_00000

iphop predict \
	--fa_file /scratch/umcg-nkuzub/iphop/iphop_pilot_filtered_shuffled_splitted/batch_00000.fna \
	--db_dir /scratch/umcg-nkuzub/iphop/iphop_database/Sept_2021_pub \
	--out_dir /scratch/umcg-nkuzub/iphop/RESULTS/batch_00000 \
	--num_threads 1
