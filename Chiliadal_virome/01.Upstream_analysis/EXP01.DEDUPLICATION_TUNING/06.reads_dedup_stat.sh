#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=reads_QC.err
#SBATCH --output=reads_QC.out
#SBATCH --mem=12gb
#SBATCH --time=4:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module purge

echo -e "SID\tEXP01 Dedup Reads 1\tEXP01 Dedup Reads 2\tEXP01 Dedup Orph 1\tEXP01 Dedup Orph 2\tEXP02 Dedup Reads 1\tEXP02 Dedup Reads 2\tEXP02 Singletons\tEXP03 Dedup Reads 1\tEXP03 Dedup Reads 2\tEXP03 Singletons\tEXP06 Dedup Reads 1\tEXP06 Dedup Reads 2\tEXP06 Singletons" >> 02_reads_dedup_stat.txt

for SAMPLE_ID in $@;
	do echo `basename ${SAMPLE_ID}` >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_paired_1.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_paired_2.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched_1.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP01/${SAMPLE_ID}_01_dedup_unmatched_2.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP02/${SAMPLE_ID}_02_dedup_repaired_1.fastq|wc -l)/4|bc >> tmp
        echo $(cat ../${SAMPLE_ID}/clean_reads/EXP02/${SAMPLE_ID}_02_dedup_repaired_2.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP02/${SAMPLE_ID}_02_dedup_singletons.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_repaired_1.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_repaired_2.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP03/${SAMPLE_ID}_03_dedup_singletons.fastq|wc -l)/4|bc >> tmp
	echo $(cat ../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_repaired_1.fastq|wc -l)/4|bc >> tmp
        echo $(cat ../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_repaired_2.fastq|wc -l)/4|bc >> tmp
        echo $(cat ../${SAMPLE_ID}/clean_reads/EXP06/${SAMPLE_ID}_06_dedup_singletons.fastq|wc -l)/4|bc >> tmp
	less tmp | paste -s >> tmp2
	paste -d "\n" tmp2 >> 02_reads_dedup_stat.txt
	rm tmp*
done

module list

module purge

