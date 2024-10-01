#!/bin/bash
#SBATCH --job-name=rAssQuality
#SBATCH --output=./out/12.rmc/VD_N_vir_g_Chiliadal_%A_%a.out
#SBATCH --mem=128gb
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

DB_index=$2

RM_DIR=$3

mkdir -p ${TMPDIR}/${SAMPLE_ID}/ALIGN

# --- LOAD MODULES --- 
module purge
module load Bowtie2
module load SAMtools
module load BEDTools
module list

# --- READ MAPPING --- 
bowtie2 \
	--very-sensitive \
	-x ${DB_index} \
	-1 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz \
	-2 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz \
        -U ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched.fastq.gz \
	--no-unal \
	--threads ${SLURM_CPUS_PER_TASK} \
	-S ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sam #\
	#&> ${RM_DIR}/alignment_log/${SAMPLE_ID}.bowtie2.log

# --- CONVERTING SAM TO BAM ---
samtools view \
        -bS \
        ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sam \
        > ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.bam

# --- SORTING BAM ---
samtools sort \
	${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1)) \
	-o ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sorted.bam

# --- INDEXING BAM ---
samtools index \
	${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sorted.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1))

cp ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sorted.bam* ${RM_DIR}/bams/

# --- GETTING N MAPPED READS AND VIRUS CONTIG COVERAGE ---
bedtools coverage \
	-a ${DB_index}.bed \
	-b ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}.sorted.bam \
	> ${RM_DIR}/coverage/${SAMPLE_ID}.coverage.txt

module purge
