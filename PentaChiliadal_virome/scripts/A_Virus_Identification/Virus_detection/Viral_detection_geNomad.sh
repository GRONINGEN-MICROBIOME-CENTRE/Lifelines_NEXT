#!/bin/bash
#SBATCH --job-name=Viral_detection_geNomad_PENTA
#SBATCH --output=Viral_detection_geNomad_PENTA_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

sample_dir=$1 #directory with the metaSPADEs scaffolds
sample_list=${sample_dir}/$2 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo '---- RUNNING geNomad WITH '${SAMPLE_ID}' SAMPLE ----'

mkdir -p ${sample_dir}/${SAMPLE_ID}/virome_discovery

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /home2/p304845/Conda_envs/geNomad_conda/; conda list 

# Run geNomad
genomad end-to-end \
    --enable-score-calibration \
	  --disable-find-proviruses \
    ${sample_dir}/${SAMPLE_ID}_metaspades_contigs.fa \
    ${sample_dir}/${SAMPLE_ID}/virome_discovery/geNomad_results \
    --cleanup \
	  --threads ${SLURM_CPUS_PER_TASK} \
    /scratch/hb-llnext/databases/geNomad_db/

# Remove intermediate files
rsync -av ${sample_dir}/${SAMPLE_ID}/virome_discovery/geNomad_results/${SAMPLE_ID}_metaspades_contigs_summary/*metaspades_contigs_virus.fna ${sample_dir}/${SAMPLE_ID}/virome_discovery/geNomad/
rsync -av ${sample_dir}/${SAMPLE_ID}/virome_discovery/geNomad_results/${SAMPLE_ID}_metaspades_contigs_summary/*metaspades_contigs_virus_summary.tsv ${sample_dir}/${SAMPLE_ID}/virome_discovery/geNomad/
rm -r ${sample_dir}/${SAMPLE_ID}/virome_discovery/geNomad_results

conda deactivate
