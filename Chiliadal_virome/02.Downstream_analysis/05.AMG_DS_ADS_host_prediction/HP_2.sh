#!/bin/bash

#SBATCH --job-name=HP_2
#SBATCH --output=./out/HP/HP_2_%A_%a.out
#SBATCH --error=./err/HP/HP_2_%A_%a.err
#SBATCH --time=03:59:00
#SBATCH --mem=55GB
#SBATCH --cpus-per-task=5
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

BATCH_LIST=$1
BATCH_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${BATCH_LIST})
OUTPUT_DIR=$2

echo ${BATCH_LIST}
echo "BATCH_ID=${BATCH_ID}"
echo "Output directory: $(realpath "$OUTPUT_DIR")"

module purge
module load Anaconda3/2022.05
conda activate iphop133_env

module list
conda list

mkdir ${TMPDIR}/${BATCH_ID}

iphop predict \
	--fa_file "${OUTPUT_DIR}/splitted_fasta/${BATCH_ID}.fna" \
	--db_dir /scratch/hb-llnext/databases/iphop133_db/Aug_2023_pub_rw \
	--out_dir "${TMPDIR}/${BATCH_ID}" \
	--num_threads ${SLURM_CPUS_PER_TASK}

mv "${TMPDIR}/${BATCH_ID}/Date_and_version.log" "${OUTPUT_DIR}/host_predictions/Date_and_version_${BATCH_ID}.log"
mv "${TMPDIR}/${BATCH_ID}/Host_prediction_to_genome_m90.csv" "${OUTPUT_DIR}/host_predictions/Host_prediction_to_genome_m90_${BATCH_ID}.csv"
mv "${TMPDIR}/${BATCH_ID}/Host_prediction_to_genus_m90.csv" "${OUTPUT_DIR}/host_predictions/Host_prediction_to_genus_m90_${BATCH_ID}.csv"
mv "${TMPDIR}/${BATCH_ID}/Detailed_output_by_tool.csv" "${OUTPUT_DIR}/host_predictions/Detailed_output_by_tool_${BATCH_ID}.csv"

rm -r ${TMPDIR}/${BATCH_ID}

conda deactivate
module purge
