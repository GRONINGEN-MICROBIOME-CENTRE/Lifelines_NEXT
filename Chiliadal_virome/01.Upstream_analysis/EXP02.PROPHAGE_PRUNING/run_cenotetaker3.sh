#!/bin/bash
#SBATCH --job-name=ProPhageBench
#SBATCH --error=./err/CT3.err
#SBATCH --output=./out/CT3.out
#SBATCH --mem=32gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module purge
module load Anaconda3
conda activate /scratch/p282752/tools/conda_envs/ct3.2.1_env

# --- RUNNING Cenote-Taker3 ---
echo "> Running Cenote-Taker3"

cenotetaker3 \
	-c ../all_positive_check.fasta \
	-wd ../ \
	-r CenoteTaker3 \
	-p T \
	--genbank F \
	-t ${SLURM_CPUS_PER_TASK}

rm -r ../CenoteTaker3/ct_processing # intermediate
rm -r ../CenoteTaker3/sequin_and_genome_maps # we will re-create it afterwards for selected contigs

conda list
conda deactivate

module list

module purge
