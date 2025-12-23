#!/bin/bash
#SBATCH --job-name=inStrain_profile
#SBATCH --output=inStrain_profile_%A_%a.out
#SBATCH --mem=100gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

mapping_dir=$1 #path to directory with the mapping BAM files
file_list=${mapping_dir}/$2 #file with the list of mapping files
vOTU_genomes=$3 #path to FASTA with vOTU representative genomes used for mapping
mapping_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file_list})
SAMPLE_ID=$(basename "$mapping_file" | cut -d . -f1)

echo '-------------------- WORKING WITH '${SAMPLE_ID}' MAPPING FILE --------------------'

# Activate environment
ml purge; module load Python/3.10.8-GCCcore-12.2.0; ml list
source /scratch/hb-llnext/python_venvs/inStrain/bin/activate

mkdir -p 2_INSTRAIN_PROFILES OUTPUT_files

# Run inStrain
inStrain profile \
    ${mapping_dir}/${mapping_file} \
    $vOTU_genomes \
    -o 2_INSTRAIN_PROFILES/${SAMPLE_ID}_inStrain \
    --min_mapq 0 \
    --min_insert 160 \
    --skip_plot_generation \
    -p ${SLURM_CPUS_PER_TASK} 
