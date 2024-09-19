#!/bin/bash

#SBATCH --job-name=12.a.HP_NCP
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=./out/12a.hps/12.a.HP_NCP.out
#SBATCH --error=./err/12a.hps/12.a.HP_NCP.err
#SBATCH --mem=10GB

module purge; ml Anaconda3/2022.05; conda activate iphop133_env
cd /scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/host_prediction_w_neg_der95_NCP
mkdir ./splitted_fasta

iphop split \
	--input_file ../virus_contigs/NCP_vOTU_representatives_w_neg_der95.fasta \
	--split_dir ./splitted_fasta \
	--n_seq 500

ls splitted_fasta > batch_list.txt
sed -i 's/.*batch_//g' batch_list.txt
sed -i 's/.fna//g' batch_list.txt

