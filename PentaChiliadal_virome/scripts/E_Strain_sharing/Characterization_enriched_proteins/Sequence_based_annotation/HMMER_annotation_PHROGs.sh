#!/bin/bash
#SBATCH --job-name=HMMER_annotation_PHROGs
#SBATCH --output=HMMER_annotation_PHROGs_%A_%a.out
#SBATCH --mem=40gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

HMM_DB=$1 #directory with HMM profiles of PHROGs
proteins_dir=$2 #directory with the split files of predicted viral proteins
file_list=${proteins_dir}/$3 #file with the list of all files in the directory
FILE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file_list})

mkdir -p HMMER_RESULTS_PHROGs OUTPUT_files

echo -e '\n-------------------- WORKING WITH '${FILE_ID}' FILE ------------------'

# Clean environment and load modules 
module purge; ml HMMER; module list

# Analyze the viral protein sequences using our HMM DB
hmmsearch \
    -o HMMER_RESULTS_PHROGs/${FILE_ID}_HMMER_result \
    --tblout HMMER_RESULTS_PHROGs/${FILE_ID}_HMMER_table \
    -E 1e-5 \
    --cpu ${SLURM_CPUS_PER_TASK} \
    ${HMM_DB}/PHROGS_new.hmm \
    ${proteins_dir}/${FILE_ID}
