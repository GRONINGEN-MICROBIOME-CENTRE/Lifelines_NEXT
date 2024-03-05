#!/bin/bash
#SBATCH --job-name=rAssQuality
#SBATCH --error=./err/11.rmc/VD_Chiliadal_%A_%a.err
#SBATCH --output=./out/11.rmc/VD_Chiliadal_%A_%a.out
#SBATCH --mem=96gb
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
module load BEDTools
module list

# --- READ MAPPING --- 
bowtie2 \
	--very-sensitive \
	-x ../VIR_DB/virus_contigs/w_neg_der95_index/w_neg_der95 \
	-1 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz \
	-2 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz \
        -U ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched.fastq.gz \
	--no-unal \
	--threads ${SLURM_CPUS_PER_TASK} \
	-S ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_w_neg_der95.sam \
	&> ../VIR_DB/mapping/VLP_to_w_neg_der95/alignment_log/${SAMPLE_ID}.bowtie2.log

# --- CONVERTING SAM TO BAM ---
samtools view \
        -bS \
        ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_w_neg_der95.sam \
        > ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_w_neg_der95.bam

# --- SORTING BAM ---
samtools sort \
	${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_w_neg_der95.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1)) \
	-o ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_w_neg_der95.sorted.bam

# --- INDEXING BAM ---
samtools index \
	${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_w_neg_der95.sorted.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1))

# --- GETTING N MAPPED READS AND VIRUS CONTIG COVERAGE ---
bedtools coverage \
	-a ../VIR_DB/virus_contigs/w_neg_der95_index/w_neg_der95.bed \
	-b ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_w_neg_der95.sorted.bam \
	> ../VIR_DB/mapping/VLP_to_w_neg_der95/coverage/${SAMPLE_ID}.coverage.txt

module purge
