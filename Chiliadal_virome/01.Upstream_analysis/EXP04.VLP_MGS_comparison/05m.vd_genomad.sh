#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --output=./out/05.mgnd/VD_PentaChiliadal_%A_%a.out
#SBATCH --mem=64gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

SAMPLE_DIR=$2

echo ${SAMPLE_LIST}
echo ${SAMPLE_DIR}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# PREPARING TMP
mkdir -p ${TMPDIR}/${SAMPLE_ID}/geNomad

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
	${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	${TMPDIR}/${SAMPLE_ID}/geNomad \
	/scratch/p282752/databases/genomad_db

genomad --version

deactivate

# --- Moving geNomad results to scratch ---
mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/geNomad/
cp ${TMPDIR}/${SAMPLE_ID}/geNomad/${SAMPLE_ID}_contigs.min1kbp_summary/${SAMPLE_ID}_contigs.min1kbp_virus_summary.tsv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/geNomad/

module list

module purge
