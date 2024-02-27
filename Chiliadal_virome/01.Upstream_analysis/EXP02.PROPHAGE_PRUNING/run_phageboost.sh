#!/bin/bash
#SBATCH --job-name=PHB_ProPhageBench
#SBATCH --error=./err/PHB.err
#SBATCH --output=./out/PHB.out
#SBATCH --mem=16gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module purge
module load Python/3.8.16-GCCcore-11.2.0
source /scratch/p282752/tools/venvs/PhageBoost_env/bin/activate

# --- RUNNING PhageBoost ---
echo "> Running PhageBoost"

PhageBoost \
	-f ../all_positive_check.fasta \
	-o ../PhageBoost \
	-j ${SLURM_CPUS_PER_TASK} \
	-cs 1000 \
	-meta 1

deactivate

module list

module purge
