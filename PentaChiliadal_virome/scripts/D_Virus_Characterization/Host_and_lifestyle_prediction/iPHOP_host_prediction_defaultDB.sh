#!/bin/bash
#SBATCH --job-name=HP_PNCHLDL_host_prediction
#SBATCH --output=./out/HP_PNCHLDL_host_prediction_%A_%a.out
#SBATCH --error=./err/HP_PNCHLDL_host_prediction_%A_%a.err
#SBATCH --time=71:59:00
#SBATCH --mem=70GB
#SBATCH --cpus-per-task=5
#SBATCH --open-mode=truncate

BATCH_LIST=$1
echo ${BATCH_LIST}

BATCH_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${BATCH_LIST})
echo "BATCH_ID=${BATCH_ID}"

module purge; ml Anaconda3/2022.05; conda activate iphop133_env

cd /scratch/p309176/host_prediction_pchldl
mkdir -p ./results

mkdir ${TMPDIR}/${BATCH_ID}

iphop predict \
	--fa_file ./splitted_fasta/${BATCH_ID}.fna \
	--db_dir /scratch/hb-llnext/databases/iphop133_db/Aug_2023_pub_rw \
	--out_dir ${TMPDIR}/${BATCH_ID} \
	--num_threads 5

mv ${TMPDIR}/${BATCH_ID}/Date_and_version.log ./results/Date_and_version_${BATCH_ID}.log
mv ${TMPDIR}/${BATCH_ID}/Host_prediction_to_genome_m90.csv ./results/Host_prediction_to_genome_m90_${BATCH_ID}.csv
mv ${TMPDIR}/${BATCH_ID}/Host_prediction_to_genus_m90.csv ./results/Host_prediction_to_genus_m90_${BATCH_ID}.csv
mv ${TMPDIR}/${BATCH_ID}/Detailed_output_by_tool.csv ./results/Detailed_output_by_tool_${BATCH_ID}.csv

