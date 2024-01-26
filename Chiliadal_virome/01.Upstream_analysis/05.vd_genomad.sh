#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --error=./err/05.gnd/VD_Chiliadal_%A_%a.err
#SBATCH --output=./out/05.gnd/VD_Chiliadal_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

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
	../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	../SAMPLES/${SAMPLE_ID}/virome_discovery/geNomad \
	/scratch/p282752/databases/genomad_db

genomad --version

deactivate

module list

module purge
