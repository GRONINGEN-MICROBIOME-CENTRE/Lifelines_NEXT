#!/bin/bash
#SBATCH --job-name=cpn60db
#SBATCH --error=./err/03.cpn/QC_Chiliadal_%A_%a.err
#SBATCH --output=./out/03.cpn/QC_Chiliadal_%A_%a.out
#SBATCH --mem=12gb
#SBATCH --time=07:29:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

mkdir -p ${TMPDIR}/${SAMPLE_ID}/ALIGN

# --- LOAD MODULES --- 
module purge
module load Bowtie2
module list

# --- ALIGING READS TO UT OF CPN60DB ---
echo "> Running read alignment"

bowtie2 \
        --very-sensitive \
	-x /scratch/p282752/databases/cpn60db_apr_24/cpn60db \
	-1 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_1.fastq.gz \
	-2 ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_paired_2.fastq.gz \
	-U ../SAMPLES/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_dedup_unmatched.fastq.gz \
        --no-unal \
        --threads ${SLURM_CPUS_PER_TASK} \
	-S ${TMPDIR}/${SAMPLE_ID}/ALIGN/${SAMPLE_ID}_cpn60db.sam \
	&> ../SAMPLES/${SAMPLE_ID}/bowtie2.cpn60db.log

if [ -f ../SAMPLES/${SAMPLE_ID}/bowtie2.cpn60db.log ]; then
	echo "> cpn60db alignment has finished"
fi

module purge

