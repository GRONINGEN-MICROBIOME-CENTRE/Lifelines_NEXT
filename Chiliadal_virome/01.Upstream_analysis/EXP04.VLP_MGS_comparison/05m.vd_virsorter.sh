#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --output=./out/05.mvs2/VD_PentaChiliadal_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=12:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

SAMPLE_DIR=$2

echo ${SAMPLE_LIST}
echo ${SAMPLE_DIR}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- PREPARING TMPDIR ---
mkdir -p ${TMPDIR}/${SAMPLE_ID}/VirSorter2

# --- LOAD MODULES --- 
module purge
module load Anaconda3
source activate /scratch/hb-llnext/conda_envs/VirSorter2_env

# --- RUNNING VIRSORTER2 ---
echo "> Running VirSorter2"

virsorter run \
	-w ${TMPDIR}/${SAMPLE_ID}/VirSorter2/ \
	-i ${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	--min-length 1000 \
	--keep-original-seq \
	--include-groups "dsDNAphage,RNA,NCLDV,ssDNA,lavidaviridae" \
	--db-dir /scratch/hb-llnext/conda_envs/VirSorter2_env/db \
	-j ${SLURM_CPUS_PER_TASK} \
	all

conda list
conda deactivate

# --- PREPARING scratch ---
mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/VirSorter2

cp ${TMPDIR}/${SAMPLE_ID}/VirSorter2/final-viral-boundary.tsv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/VirSorter2/${SAMPLE_ID}_final-viral-boundary.tsv

module list

module purge
