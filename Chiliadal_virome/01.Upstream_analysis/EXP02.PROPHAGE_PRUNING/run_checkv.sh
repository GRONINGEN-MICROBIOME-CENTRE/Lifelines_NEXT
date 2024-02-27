#!/bin/bash
#SBATCH --job-name=CHV_ProPhageBench
#SBATCH --error=./err/CHV.err
#SBATCH --output=./out/CHV.out
#SBATCH --mem=48gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module purge
module load CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8
module list

# --- RUNNING CHECKV --- 
checkv end_to_end \
    ../all_positive_check.fasta \
    ../CheckV \
    -t ${SLURM_CPUS_PER_TASK} \
    -d /scratch/hb-llnext/databases/checkv-db-v1.5

module purge
