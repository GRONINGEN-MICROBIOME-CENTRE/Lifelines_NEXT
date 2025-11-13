#!/bin/bash
#SBATCH --job-name=AMGP_KFS_1
#SBATCH --output=./out/AMGP/AMGP_KFS_1_%A_%a.out
#SBATCH --error=./err/AMGP/AMGP_KFS_1_%A_%a.err
#SBATCH --mem=10gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

HMM_DB=$1 #directory with KOFam HMM profiles
VIRAL_PROTEINS_DIR=$2 #directory with the split files of predicted viral proteins
OUTPUT_DIR=$3 #full path to output directory
BATCH_LIST=$4 #file with the list of all files in the directory
BATCH_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${BATCH_LIST})

echo "WORKING WITH '${BATCH_ID}' BATCH"

module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/KOFamScan_env; conda list 

mkdir -p ${TMPDIR}/${BATCH_ID} 
mkdir -p ${OUTPUT_DIR}

# Annotate proteins
exec_annotation -p ${HMM_DB}/profiles \
    -k ${HMM_DB}/ko_list \
    --cpu ${SLURM_CPUS_PER_TASK} \
    -f detail-tsv \
    --tmp-dir ${TMPDIR}/${BATCH_ID} \
    -o ${OUTPUT_DIR}/${BATCH_ID}_ko-annotations.tsv \
    ${VIRAL_PROTEINS_DIR}/${BATCH_ID}.faa
 
 # Filter significant hits and keep only the one KO per protein (with the highest significance)
 bit-filter-KOFamScan-results -i ${OUTPUT_DIR}/${BATCH_ID}_ko-annotations.tsv \
     -o ${OUTPUT_DIR}/${BATCH_ID}_ko-annotations-filtered.tsv

 rm -r ${TMPDIR}/${BATCH_ID}
