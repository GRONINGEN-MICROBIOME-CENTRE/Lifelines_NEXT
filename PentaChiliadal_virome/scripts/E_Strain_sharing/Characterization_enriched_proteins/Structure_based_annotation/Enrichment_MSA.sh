#!/bin/bash
#SBATCH --job-name=Enrichment_MSAs
#SBATCH --output=Enrichment_MSAs_%A_%a.out
#SBATCH --mem=60gb
#SBATCH --time=00-01:00 
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular


MSA_dir=$1 #directory with MSAs of unnanotated clusters
PF_list=$2 # TXT file with the list of unnanotated clusters
DB=$3 #path to Uniref30 DB
enriched_MSAs=$4 #path to output directory

cluster_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${PF_list} | sed 's/\.faa$//')


echo -e '\n-------------------- WORKING WITH '${cluster_ID}' PROTEIN CLUSTER --------------------'

mkdir -p $enriched_MSAs $enriched_MSAs/HHR_files OUTPUT_files

hhblits \
    -d ${DB}/UniRef30_2023_02 \
    -i ${MSA_dir}/${cluster_ID}.faa \
    -o $enriched_MSAs/HHR_files/${cluster_ID}.hhr \
    -oa3m $enriched_MSAs/${cluster_ID}.a3m \
    -id 80 \
    -n 3 \
    -cpu ${SLURM_CPUS_PER_TASK} \
    -e 0.1
    
