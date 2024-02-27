#!/bin/bash
#SBATCH --job-name=VIB_ProPhageBench
#SBATCH --error=./err/VIB.err
#SBATCH --output=./out/VIB.out
#SBATCH --mem=16gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module purge
module load prodigal-gv/2.11.0-GCCcore-12.2.0 
module load Python/3.11.3-GCCcore-12.3.0

# --- PREDICTING ORFs ---
echo "> Running parallel prodigal-gv"

# https://raw.githubusercontent.com/apcamargo/prodigal-gv/master/parallel-prodigal-gv.py

python /scratch/p282752/ANALYSIS_CHILIADAL/scripts/parallel-prodigal-gv.py \
	-t ${SLURM_CPUS_PER_TASK} \
	-q \
	-i ../all_positive_check.fasta \
	-a ../all_positive_check.AA.fasta \
	-o ../prodigal.out

# --- CLEAN ENV --- 
module purge

# --- LOAD MODULES ---
module load Anaconda3/2022.05
conda activate /scratch/hb-llnext/conda_envs/Vibrant_env

# --- RUNNING VIBRANT ---
echo "> Running VIBRANT"

/scratch/hb-llnext/conda_envs/Vibrant_env/bin/VIBRANT_run.py \
	-i ../all_positive_check.AA.fasta \
	-folder ../ \
	-f prot \
	-t ${SLURM_CPUS_PER_TASK} \
	-l 1000 \
	-virome \
	-no_plot

conda list
conda deactivate

module list

module purge
