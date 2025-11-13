#!/bin/bash
#SBATCH --job-name=AMGP_PRG_1
#SBATCH --output=./out/AMGP/AMGP_PRG_1_%A_%a.out
#SBATCH --error=./err/AMGP/AMGP_PRG_1_%A_%a.err
#SBATCH --time=00:59:00
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

VIRAL_PROTEINS_DIR=$1 #directory with the split files of predicted viral proteins
OUTPUT_DIR=$2 #full path to output directory
BATCH_LIST=$3 #file with the list of all files in the directory
BATCH_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${BATCH_LIST})

mkdir -p ${OUTPUT_DIR}/HMMER_RESULTS

# Clean environment and load modules 
module purge; ml HMMER; module list

# Analyze the viral protein sequences using our HMM DB
hmmsearch \
    -o ${OUTPUT_DIR}/HMMER_RESULTS/${BATCH_ID}_HMMER_result \
    --tblout ${OUTPUT_DIR}/HMMER_RESULTS/${BATCH_ID}_HMMER_table \
    -E 0.001 \
    --cpu ${SLURM_CPUS_PER_TASK} \
    /scratch/hb-llnext/databases/PHROGs_HMM_db/PHROGS_new.hmm \
    ${VIRAL_PROTEINS_DIR}/${BATCH_ID}.faa
