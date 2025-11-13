#!/bin/bash

#SBATCH --job-name=PP_PNT_2
#SBATCH --output=./out/PP/PP_PNT_2_%A_%a.out
#SBATCH --error=./err/PP/PP_PNT_2_%A_%a.err
#SBATCH --time=01:29:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

BATCH_LIST=$1
BATCH_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${BATCH_LIST})
OUTPUT_DIR=$2 

echo ${BATCH_LIST}
echo "BATCH_ID=${BATCH_ID}"
echo "Output directory: $(realpath "$OUTPUT_DIR")"

module purge
ml Anaconda3/2024.02-1
conda activate /home1/p309176/conda/envs/phanotate_env

module list
conda list

phanotate.py "${OUTPUT_DIR}/splitted_fasta/${BATCH_ID}.fna" \
	-o "${OUTPUT_DIR}/PP_PNT/${BATCH_ID}.faa" \
	-f faa

phanotate.py "${OUTPUT_DIR}/splitted_fasta/${BATCH_ID}.fna" \
        -o "${OUTPUT_DIR}/PP_PNT/${BATCH_ID}.gff" \
        -f gff

conda deactivate
module purge

# Removing special charachters from the .faa file, so it's suitable for the other tools
 sed -i 's/[#+]$//' "${OUTPUT_DIR}/PP_PNT/${BATCH_ID}.faa"
 sed -i '/tRNA/d' "${OUTPUT_DIR}/PP_PNT/${BATCH_ID}.gff"
