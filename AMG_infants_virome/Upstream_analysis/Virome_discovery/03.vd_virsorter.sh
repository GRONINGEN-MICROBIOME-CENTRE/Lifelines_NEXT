#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --error=./err/03.vs2/VD_AMG_%A_%a.err
#SBATCH --output=./out/03.vs2/VD_AMG_%A_%a.out
#SBATCH --mem=16gb
#SBATCH --time=12:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST})

echo "SAMPLE_ID=${SAMPLE_ID}"

mkdir -p ../SAMPLES/${SAMPLE_ID}/virome_discovery/VirSorter2

# --- LOAD MODULES --- 
module purge
module load Anaconda3
source activate /scratch/hb-llnext/conda_envs/VirSorter2_env

# --- RUNNING VIRSORTER2 ---
echo "> Running VirSorter2"

virsorter run \
	-w ../SAMPLES/${SAMPLE_ID}/virome_discovery/VirSorter2/ \
	-i ../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	--min-length 1000 \
	--keep-original-seq \
	--include-groups "dsDNAphage,RNA,NCLDV,ssDNA,lavidaviridae" \
	--db-dir /scratch/hb-llnext/conda_envs/VirSorter2_env/db \
	-j ${SLURM_CPUS_PER_TASK} \
	all

# --- REMOVING BYPRODUCTS ---
echo "> Removing byproducts"
rm -r ../SAMPLES/${SAMPLE_ID}/virome_discovery/VirSorter2/iter-0
rm -r ../SAMPLES/${SAMPLE_ID}/virome_discovery/VirSorter2/log
rm ../SAMPLES/${SAMPLE_ID}/virome_discovery/VirSorter2/config.yaml

conda list
conda deactivate

module list

module purge
