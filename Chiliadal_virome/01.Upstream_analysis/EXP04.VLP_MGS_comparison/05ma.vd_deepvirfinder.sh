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

# since DVF stumbles on some extra large sequences:
cp ${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta ${TMPDIR}/${SAMPLE_ID}/DVF
/scratch/p282752/ANALYSIS_CHILIADAL/scripts/fastasplitn ${TMPDIR}/${SAMPLE_ID}/DVF/${SAMPLE_ID}_contigs.min1kbp.fasta 2

# creates 2 files: frag001.fa and frag002.fa, frag001.fa should contain the largest seq (here ~ 2.1 Mb), let's check it:
echo $("Largest contig that excluded from DVF run:" head -n 1 frag001.fa)

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
	-i frag002.fa \
	-o ${TMPDIR}/${SAMPLE_ID}/DVF \
	-l 1000

conda list
conda deactivate

# --- MOVING TO /SCRATCH ---
mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/DeepVirFinder
cat ${TMPDIR}/${SAMPLE_ID}/DVF/*_gt1000bp_dvfpred.txt ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/DeepVirFinder/${SAMPLE_ID}_contigs.min1kbp.fasta_gt1000bp_dvfpred.txt

module list

module purge
