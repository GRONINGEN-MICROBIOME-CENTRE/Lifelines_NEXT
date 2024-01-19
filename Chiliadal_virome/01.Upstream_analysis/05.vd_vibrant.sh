#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --error=./err/05.vib/VD_Chiliadal_%A_%a.err
#SBATCH --output=./out/05.vib/VD_Chiliadal_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load prodigal-gv/2.11.0-GCCcore-12.2.0 
module load Python/3.11.3-GCCcore-12.3.0

# --- PREDICTING ORFs ---
echo "> Running parallel prodigal-gv"

# https://raw.githubusercontent.com/apcamargo/prodigal-gv/master/parallel-prodigal-gv.py

python parallel-prodigal-gv.py \
	-t ${SLURM_CPUS_PER_TASK} \
	-q \
	-i ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	-a ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.AA.fasta \
	-o ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_prodigal.out

# --- CLEAN ENV --- 
module purge

# --- LOAD MODULES ---
module load Anaconda3/2022.05
conda activate /scratch/hb-llnext/conda_envs/Vibrant_env

# --- RUNNING VIBRANT ---
echo "> Running VIBRANT"

/scratch/hb-llnext/conda_envs/Vibrant_env/bin/VIBRANT_run.py \
	-i ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.AA.fasta \
	-folder ../SAMPLES/${SAMPLE_ID}/virome_discovery/ \
	-f prot \
	-t ${SLURM_CPUS_PER_TASK} \
	-l 1000 \
	-virome \
	-no_plot

if [ $(grep 'End' ../SAMPLES/${SAMPLE_ID}/virome_discovery/VIBRANT_${SAMPLE_ID}_contigs.min1kbp.AA/VIBRANT_log_run_${SAMPLE_ID}_contigs.min1kbp.AA.log | wc -l) == 1 ]; then
	echo "VIBRANT is done"
fi

# --- REMOVING BYPRODUCTS ---
echo "> Removing and renaming byproducts"
rm -r ../SAMPLES/${i}/virome_discovery/VIBRANT_CHV200013F12_contigs.min1kbp.AA/VIBRANT_HMM_tables_unformatted_CHV200013F12_contigs.min1kbp.AA

conda list
conda deactivate

module list

module purge
