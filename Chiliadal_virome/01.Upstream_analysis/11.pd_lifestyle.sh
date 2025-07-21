#!/bin/bash
#SBATCH --job-name=VIR_DB
#SBATCH --output=./out/11.lfs/PD_lifestyle_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

DB_LIST=$1

echo "DB list: ${DB_LIST}"

DB=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${DB_LIST})

echo "DB=${DB}"

# --- LOAD MODULES ---
module load Python/3.8.16-GCCcore-11.2.0
source /scratch/p282752/tools/python_envs/bacphlip/bin/activate
module load HMMER/3.3.2-gompi-2021b
module list

# --- WORKING IN TMP ---
mkdir ${TMPDIR}/${DB}_LIFESTYLE

cp /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/virus_contigs/CLEANED/fragmented/${DB}.fa ${TMPDIR}/${DB}_LIFESTYLE

# --- RUNNING BACPHLIP ---
bacphlip --multi_fasta -i ${TMPDIR}/${DB}_LIFESTYLE/${DB}.fa

# --- COPYING RESULTS ---
cp ${TMPDIR}/${DB}_LIFESTYLE/${DB}.fa.bacphlip /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/virus_contigs/CLEANED/fragmented/

deactivate

module purge
