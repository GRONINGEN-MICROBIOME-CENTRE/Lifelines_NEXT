#!/bin/bash
#SBATCH --job-name=inStrain_profile
#SBATCH --output=./out/13.dct/IS_profile_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST})

echo "SAMPLE_ID=${SAMPLE_ID}"

# COPYING THE READ ALIGNMENT MAPS:
cp ../VIR_DB/mapping/VLP_N_vir_genes/bams/${SAMPLE_ID}.sorted.bam* ../VIR_DB/decontamination/bams/

# FILTERING READ ALIGNMENT MAPS TO CONTAIN ONLY NC-IDENTIFIED VOTUS PRESENT PER FILE:
ml SAMtools
samtools view \
	-b \
	-o ../VIR_DB/decontamination/bams/${SAMPLE_ID}_filtered.sorted.bam \
	../VIR_DB/decontamination/bams/${SAMPLE_ID}.sorted.bam $(cat ../VIR_DB/decontamination/keep_vOTUs_per_sample/${SAMPLE_ID}_vOTUs_to_keep)

rm ../VIR_DB/decontamination/bams/${SAMPLE_ID}.sorted.bam*

samtools index \
	../VIR_DB/decontamination/bams/${SAMPLE_ID}_filtered.sorted.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1))

module purge

# --- LOADING MODULES ---
ml Python/3.10.8-GCCcore-12.2.0
source /scratch/p282752/tools/python_envs/instrain/bin/activate
ml SAMtools

inStrain profile \
	../VIR_DB/decontamination/bams/${SAMPLE_ID}_filtered.sorted.bam \
	../VIR_DB/TMP_EXP_DREP/NEXT_vOTUr.fasta \
	-o ../VIR_DB/decontamination/ALL_PROFILE/${SAMPLE_ID}.inStrain \
	--database_mode \
	-s ../VIR_DB/decontamination/genomes.stb \
	-p ${SLURM_CPUS_PER_TASK} \
	--scaffolds_to_profile ../VIR_DB/decontamination/vOTUs_shared_by_samples_and_ncs \
	--min_cov 1 \
	--skip_plot_generation

deactivate

module purge
