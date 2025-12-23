#!/bin/bash
#SBATCH --job-name=Anvio_profile
#SBATCH --output=Anvio_profile_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

mapping_dir=$1 #path to directory with the mapping BAM files
file_list=$2 #file with the list of all mapping files
vOTU_DB=$3 #path to the vOTU DB built in the previous step
mapping_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${file_list})
SAMPLE_ID=$(basename "$mapping_file" | cut -d . -f1)

echo '-------------------- WORKING WITH '${SAMPLE_ID}' MAPPING FILE --------------------'

# Clean environment, load modules 
ml purge; ml Anaconda3; ml list
conda activate /scratch/hb-llnext/conda_envs/anvio-8; conda list

mkdir -p 1_ANVIO_PROFILES OUTPUT_files

# Run anvi-profile (including codon frequency characterization)
anvi-profile \
    -i ${mapping_dir}/${mapping_file} \
    -c $vOTU_DB \
    --profile-SCVs \
    --min-coverage-for-variability 10 \
    --sample-name ${SAMPLE_ID} \
    -o 1_ANVIO_PROFILES/${SAMPLE_ID}_profile \
    -T ${SLURM_CPUS_PER_TASK} 
