#!/bin/bash
#SBATCH --job-name=Viral_detection_VS2_PENTA
#SBATCH --output=Viral_detection_VS2_PENTA_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

sample_dir=$1 #directory with the metaSPADEs scaffolds
sample_list=${sample_dir}/$2 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo '---- RUNNING VirSorter2 WITH '${SAMPLE_ID}' SAMPLE ----'

mkdir -p ${sample_dir}/${SAMPLE_ID}/virome_discovery

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/VirSorter2_env; conda list

# Run VirSorter2
virsorter run \
	-w ${sample_dir}/${SAMPLE_ID}/virome_discovery/VirSorter2 \
	-i ${sample_dir}/${SAMPLE_ID}_metaspades_contigs.fa \
	--min-length 10000 \
	--keep-original-seq \
	--include-groups "dsDNAphage,NCLDV,ssDNA,lavidaviridae" \
	--db-dir /scratch/hb-llnext/conda_envs/VirSorter2_env/db \
	-j ${SLURM_CPUS_PER_TASK} \
	all

rm -r ${sample_dir}/${SAMPLE_ID}/virome_discovery/VirSorter2/iter-0
rm -r ${sample_dir}/${SAMPLE_ID}/virome_discovery/VirSorter2/log
rm ${sample_dir}/${SAMPLE_ID}/virome_discovery/VirSorter2/config.yaml

conda deactivate
