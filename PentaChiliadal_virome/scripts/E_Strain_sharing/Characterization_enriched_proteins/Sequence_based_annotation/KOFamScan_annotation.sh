#!/bin/bash
#SBATCH --job-name=KOFamScan_annotation
#SBATCH --output=KOFamScan_annotation_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

HMM_DB=$1 #directory with KOFam HMM profiles
proteins_dir=$2 #directory with the split files of predicted viral proteins
file_list=${proteins_dir}/$3 #file with the list of all files in the directory
FILE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file_list})

echo -e '\n-------------------- WORKING WITH '${FILE_ID}' FILE ------------------'
    
echo -e '\n---- RUNNING KOFamScan annotation ----'
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/KOFamScan_env; conda list 

mkdir -p ${TMPDIR}/${FILE_ID} RESULTS

# Annotate proteins
exec_annotation -p ${HMM_DB}/profiles \
    -k ${HMM_DB}/ko_list \
    --cpu ${SLURM_CPUS_PER_TASK} \
    -f detail-tsv \
    --tmp-dir ${TMPDIR}/${FILE_ID} \
    -o RESULTS/${FILE_ID}_ko-annotations.tsv \
    ${proteins_dir}/${FILE_ID}
 
 # Filter significant hits and keep only the one KO per protein (with the highest significance)
 #bit-filter-KOFamScan-results -i RESULTS/${FILE_ID}_ko-annotations.tsv \
 #    -o RESULTS/${FILE_ID}_ko-annotations-filtered.tsv
     
 rm -r ${TMPDIR}/${FILE_ID}

