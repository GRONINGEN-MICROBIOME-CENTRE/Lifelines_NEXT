#!/bin/bash
#SBATCH --job-name=rAssQuality
#SBATCH --error=./err/09.pra/VD_Chiliadal_%A_%a.err
#SBATCH --output=./out/09.pra/VD_Chiliadal_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

mkdir -p ${TMPDIR}/${SAMPLE_ID}/ALIGN

# --- LOAD MODULES --- 
module purge
module load Bowtie2
module load SAMtools
module list

#bowtie2-build \
#	../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
#	${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}

#bowtie2 \
#	--very-sensitive \
#	-x ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID} \
#	-1 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz \
#	-2 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz \
#	--no-unal \
#	--threads ${SLURM_CPUS_PER_TASK} \
#	-S ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sam \
#	&> ../SAMPLES/${SAMPLE_ID}/paired_to_1kbp.bowtie2.log

#rm ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}*

#bowtie2-build \
#        ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta \
#        ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_vir

#bowtie2 \
#        --very-sensitive \
#        -x ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_vir \
#        -1 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz \
#        -2 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz \
#        --no-unal \
#        --threads ${SLURM_CPUS_PER_TASK} \
#        -S ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sam \
#        &> ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/paired_to_all.vir.bowtie2.log

#rm ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}*

bowtie2-build \
        ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta \
	${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_vir_ext

bowtie2 \
        --very-sensitive \
        -x ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_vir_ext \
        -1 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz \
        -2 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz \
        --no-unal \
        --threads ${SLURM_CPUS_PER_TASK} \
        -S ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sam \
        &> ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/paired_to_ext.vir.bowtie2.log

rm ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}*

bowtie2-build \
        ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_pruned_viral.fasta \
        ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_vir_ext_prun

bowtie2 \
        --very-sensitive \
        -x ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_vir_ext_prun \
        -1 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz \
        -2 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz \
        --no-unal \
        --threads ${SLURM_CPUS_PER_TASK} \
        -S ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sam \
        &> ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/paired_to_ext.prun.vir.bowtie2.log

module purge
