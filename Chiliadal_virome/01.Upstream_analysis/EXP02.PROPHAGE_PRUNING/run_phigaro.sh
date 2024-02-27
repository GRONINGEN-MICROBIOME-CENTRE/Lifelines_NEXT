#!/bin/bash
#SBATCH --job-name=PHI_ProPhage_Bench
#SBATCH --error=./err/PHI.err
#SBATCH --output=./out/PHI.out
#SBATCH --mem=16gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module purge

# --- RUNNING PHIGARO ---
echo "> Running Phigaro"

apptainer exec -B /scratch/p282752 /scratch/p282752/tools/phigaro_latest.sif \
	/root/miniconda3/bin/phigaro \
	-f /scratch/p282752/benchmark_prophage/all_positive_check.fasta \
	-c /root/.phigaro/config.yml \
	-d \
	-e tsv \
	-o /scratch/p282752/benchmark_prophage/Phigaro \
	-t ${SLURM_CPUS_PER_TASK} \
	--save-fasta

module list

module purge
