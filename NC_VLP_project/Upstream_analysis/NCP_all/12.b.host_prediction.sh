#!/bin/bash
#SBATCH --job-name=12.b.HP_NCP
#SBATCH --output=./out/12b.hpr/12.b.HP_NCP_%A_%a.out
#SBATCH --error=./err/12b.hpr/12.b.HP_NCP_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=70GB
#SBATCH --cpus-per-task=5
#SBATCH --open-mode=truncate

BATCH_LIST=$1
echo ${BATCH_LIST}

BATCH_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${BATCH_LIST})
echo "BATCH_ID=${BATCH_ID}"

module purge; ml Anaconda3/2022.05; conda activate iphop133_env

cd /scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/host_prediction_w_neg_der95_NCP
# mkdir -p ./results

mkdir ${TMPDIR}/batch_${BATCH_ID}

iphop predict \
	--fa_file ./splitted_fasta/batch_${BATCH_ID}.fna \
	--db_dir /scratch/hb-llnext/databases/iphop133_db/Aug_2023_pub_rw \
	--out_dir ${TMPDIR}/batch_${BATCH_ID} \
	--num_threads 5

mv ${TMPDIR}/batch_${BATCH_ID}/Date_and_version.log ./results/Date_and_version_batch_${BATCH_ID}.log
mv ${TMPDIR}/batch_${BATCH_ID}/Host_prediction_to_genome_m90.csv ./results/Host_prediction_to_genome_m90_batch_${BATCH_ID}.csv
mv ${TMPDIR}/batch_${BATCH_ID}/Host_prediction_to_genus_m90.csv ./results/Host_prediction_to_genus_m90_batch_${BATCH_ID}.csv
mv ${TMPDIR}/batch_${BATCH_ID}/Detailed_output_by_tool.csv ./results/Detailed_output_by_tool_batch_${BATCH_ID}.csv
