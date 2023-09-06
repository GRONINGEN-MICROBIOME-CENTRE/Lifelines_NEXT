#!/bin/bash
#SBATCH --job-name=reads_QC
#SBATCH --error=./err/assembly_qual_07.err
#SBATCH --output=./out/assembly_qual_07.out
#SBATCH --mem=12gb
#SBATCH --time=4:59:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate


# --- LOAD MODULES --- 
module purge

for EXP_ID in $(cat EXP_list); do

	echo -e "SID\t${EXP_ID} insert size\t${EXP_ID} total size contigs >= 10 kb\t${EXP_ID} N contigs >= 1kb\t${EXP_ID} Largest contig\t${EXP_ID} N50" >> ${EXP_ID}_03_assembly_stat.txt
	
	for SAMPLE_ID in $(cat ../List_VLP_samples.txt); do
	
		echo `basename ${SAMPLE_ID}` >> tmp
		grep 'Insert size = ' ../${SAMPLE_ID}/assembly/${EXP_ID}/spades.log | awk 'NR==1' | awk -F ' = ' '{print $2}' | awk -F ',' '{print $1}' >> tmp
		grep 'Total length (>= 10000 bp)' ../${SAMPLE_ID}/assembly/${EXP_ID}/quast/report.tsv | awk -F '\t' '{print $2}' >> tmp
		grep 'contigs (>= 1000 bp)' ../${SAMPLE_ID}/assembly/${EXP_ID}/quast/report.tsv | awk -F '\t' '{print $2}' >> tmp
		grep 'Largest contig' ../${SAMPLE_ID}/assembly/${EXP_ID}/quast/report.tsv | awk -F '\t' '{print $2}' >> tmp
		grep 'N50' ../${SAMPLE_ID}/assembly/${EXP_ID}/quast/report.tsv | awk -F '\t' '{print $2}' >> tmp
		less tmp | paste -s >> tmp2
		paste -d "\n" tmp2 >> ${EXP_ID}_03_assembly_stat.txt
		rm tmp*
	
	done
done
module list

module purge

