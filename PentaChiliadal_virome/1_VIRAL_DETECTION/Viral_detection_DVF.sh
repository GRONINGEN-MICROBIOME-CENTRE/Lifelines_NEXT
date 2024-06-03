#!/bin/bash
#SBATCH --job-name=Viral_detection_DVF_PENTA
#SBATCH --output=Viral_detection_DVF_PENTA_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

sample_dir=$1 #directory with the metaSPAdes contigs/scaffolds
sample_list=${sample_dir}/$2 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo '---- RUNNING DeepVirFinder WITH '${SAMPLE_ID}' SAMPLE ----'

mkdir -p ${sample_dir}/${SAMPLE_ID}/virome_discovery

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/DeepVirFinder_env; conda list

# Run DeepVirFinder
python3 /scratch/hb-llnext/conda_envs/DeepVirFinder_env/DeepVirFinder/dvf.py \
	-i ${sample_dir}/${SAMPLE_ID}_metaspades_contigs.fa \
	-o ${sample_dir}/${SAMPLE_ID}/virome_discovery/DeepVirFinder \
	-l 10000

conda deactivate
