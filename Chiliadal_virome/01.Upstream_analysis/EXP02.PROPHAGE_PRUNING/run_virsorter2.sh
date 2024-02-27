#!/bin/bash
#SBATCH --job-name=VS2
#SBATCH --error=./err/VS2.err
#SBATCH --output=./out/VS2.out
#SBATCH --mem=16gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

mkdir -p ../VirSorter2

# --- LOAD MODULES --- 
module purge
module load Anaconda3
source activate /scratch/hb-llnext/conda_envs/VirSorter2_env
which python
# --- RUNNING VIRSORTER2 ---
echo "> Running VirSorter2"

virsorter run \
	-w ../VirSorter2/ \
	-i ../all_positive_check.fasta \
	--min-length 1000 \
	--keep-original-seq \
	--include-groups "dsDNAphage,RNA,NCLDV,ssDNA,lavidaviridae" \
	--db-dir /scratch/hb-llnext/conda_envs/VirSorter2_env/db \
	-j ${SLURM_CPUS_PER_TASK} \
	all

# --- REMOVING BYPRODUCTS ---
#echo "> Removing byproducts"
#rm -r ../VirSorter2/iter-0
#rm -r ../VirSorter2/log
#rm ../VirSorter2/config.yaml

conda list
conda deactivate

module list

module purge
