#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --output=./out/05.mdvf/VD_PentaChiliadal_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=12:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

SAMPLE_DIR=$2

echo ${SAMPLE_LIST}
echo ${SAMPLE_DIR}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- SWITCH DVF run TO TMP --- 
mkdir -p ${TMPDIR}/${SAMPLE_ID}/DVF
cd ${TMPDIR}/${SAMPLE_ID}/DVF

# --- LOAD MODULES --- 
module purge
module load Anaconda3
source activate /scratch/hb-llnext/conda_envs/DeepVirFinder_env

# --- SWITCH Theano to TMP ---
export THEANO_FLAGS="base_compiledir=${TMPDIR}/${USER}/${SAMPLE_ID}/theano"
mkdir -p /${TMPDIR}/${USER}/${SAMPLE_ID}/theano
python -c "import theano; print(theano.config.base_compiledir)"

# --- RUNNING DeepVirFinder ---
echo "> Running DeepVirFinder"

python /scratch/hb-llnext/conda_envs/DeepVirFinder_env/DeepVirFinder/dvf.py \
	-i ${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	-o ${TMPDIR}/${SAMPLE_ID}/DVF \
	-l 1000

conda list
conda deactivate

# --- MOVING TO /SCRATCH ---
mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/DeepVirFinder
cp ${TMPDIR}/${SAMPLE_ID}/DVF/${SAMPLE_ID}_contigs.min1kbp.fasta_gt1000bp_dvfpred.txt ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/DeepVirFinder/

module list

module purge
