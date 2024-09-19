#!/bin/bash
#SBATCH --job-name=rAssQuality_decontam_99_NCP
#SBATCH --error=./err/09.rmc_NCP/VD_NCP_decontam_99_%A_%a.err
#SBATCH --output=./out/09.rmc_NCP/VD_NCP_decontam_99_%A_%a.out
#SBATCH --mem=256gb
#SBATCH --time=08:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST})

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
	-x /scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/virus_contigs/der95_index_NCP_decontam_99/der95_wo_nc99_NCP \
	-1 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
	-2 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
        -U ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_unmatched.fastq.gz \
	--no-unal \
	--threads ${SLURM_CPUS_PER_TASK} \
	-S ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_der95_wo_nc99_NCP.sam \
	&> /scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/mapping/VLP_to_wo_nc99_der95_NCP/alignment_log_NCP/${SAMPLE_ID}_NCP.bowtie2.log

# --- CONVERTING SAM TO BAM ---
samtools view \
        -bS \
        ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_der95_wo_nc99_NCP.sam \
        > ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_der95_wo_nc99_NCP.bam

# --- SORTING BAM ---
samtools sort \
	${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_der95_wo_nc99_NCP.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1)) \
	-o ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_der95_wo_nc99_NCP.sorted.bam

# --- INDEXING BAM ---
samtools index \
	${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_der95_wo_nc99_NCP.sorted.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1))

# --- GETTING N MAPPED READS AND VIRUS CONTIG COVERAGE ---
bedtools coverage \
	-a /scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/virus_contigs/der95_index_NCP_decontam_99/der95_wo_nc99_NCP.bed \
	-b ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_der95_wo_nc99_NCP.sorted.bam \
	> /scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/mapping/VLP_to_wo_nc99_der95_NCP/coverage_NCP/${SAMPLE_ID}_NCP.coverage.txt

echo "${SAMPLE_ID} is done."

module purge