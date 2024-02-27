#!/bin/bash
#SBATCH --job-name=gen_ProPhageBench
#SBATCH --error=./err/GND.err
#SBATCH --output=./out/GND.out
#SBATCH --mem=32gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module load ARAGORN/1.2.41-foss-2021b
module load Python/3.9.5-GCCcore-10.3.0
source /scratch/p282752/tools/python_envs/geNomad/bin/activate

# --- RUNNING geNomad ---
echo "> Running geNomad"

genomad \
	end-to-end \
	--enable-score-calibration \
	--cleanup \
	../all_positive_check.fasta \
	../geNomad \
	/scratch/p282752/databases/genomad_db

genomad --version

deactivate

module list

module purge
