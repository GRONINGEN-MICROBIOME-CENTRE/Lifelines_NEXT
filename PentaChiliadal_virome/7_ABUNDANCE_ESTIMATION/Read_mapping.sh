#!/bin/bash
#SBATCH --job-name=Read_mapping
#SBATCH --output=Read_mapping_%A_%a.out
#SBATCH --mem=80gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

paired_reads=$1 #path to the batch folder with the quality trimmed paired reads
unpaired_reads=$2 #path to the batch folder with the quality trimmed unpaired reads
sample_list=${paired_reads}/$3 #file with the list of samples
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})
Bowtie_DB=$4 # Path to Bowtie index generated from the representative sequences
viral_rep_seqs=$5 # FASTA file with all the vOTUs representative viral sequences

mkdir -p RESULTS OUTPUT_files
mkdir -p RESULTS/Alignment_results RESULTS/Breadth_coverage_results

echo '-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

# Clean environment, load modules
module purge; ml Python/3.10.8-GCCcore-12.2.0 Bowtie2 SAMtools BEDTools; module list

# Map the reads
bowtie2 \
	--very-sensitive \
	-x $Bowtie_DB/Viral_rep_seqs_DB \
	-1 $paired_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
	-2 $paired_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
  -U $unpaired_reads/${SAMPLE_ID}_kneaddata_unmatched.fastq.gz \
	--no-unal \
	--threads ${SLURM_CPUS_PER_TASK} \
	-S RESULTS/${SAMPLE_ID}_all_vir_alignments.sam

samtools view \
	-S RESULTS/${SAMPLE_ID}_all_vir_alignments.sam \
	-b \
	-o RESULTS/${SAMPLE_ID}_all_vir_alignments.bam

samtools sort \
	RESULTS/${SAMPLE_ID}_all_vir_alignments.bam \
	-o RESULTS/Alignment_results/${SAMPLE_ID}.sorted.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1))

samtools index \
	-@ $((${SLURM_CPUS_PER_TASK}-1)) \
	RESULTS/Alignment_results/${SAMPLE_ID}.sorted.bam

# Get coverage final tables 
bedtools coverage \
	-a Viral_rep_seqs_DB.bed \
	-b RESULTS/Alignment_results/${SAMPLE_ID}.sorted.bam \
	> RESULTS/Breadth_coverage_results/${SAMPLE_ID}.coverage.txt

# Remove intermediate files
rm RESULTS/${SAMPLE_ID}_all_vir_alignments.sam RESULTS/${SAMPLE_ID}_all_vir_alignments.bam
