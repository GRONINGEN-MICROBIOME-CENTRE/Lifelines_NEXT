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

echo -e "SID\tRaw Reads 1\tRaw Reads 2\tAdapTr Reads 1\tAdapTr Read 2\tNo Tail 1\tNo Tail 2\tQualTr 1\tQualTr 2\tQualTr Orph 1\tQualTr Orph 2\tHuman reads 1\tHuman reads 2\tHuman Orph 1\tHuman Orph 2\tClean reads 1\tClean reads 2\tClean Orph 1\tClean Orph 2" >> 01_reads_stat.txt

for SAMPLE_ID in $@;
	do echo `basename ${SAMPLE_ID}` >> tmp
	grep 'Input:' ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | head -1 | awk -F ' ' '{print $2 / 2}' >> tmp
	grep 'Input:' ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | head -1 | awk -F ' ' '{print $2 / 2}' >> tmp
	grep 'Result:' ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | head -1 | awk -F ' ' '{print $2 / 2}' >> tmp
	grep 'Result:' ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | head -1 | awk -F ' ' '{print $2 / 2}' >> tmp
	grep -m 2 'Result:' ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | awk 'NR==2 {print $2}' >> tmp
	grep -m 3 'Result:' ../${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log | awk 'NR==3 {print $2}' >> tmp
	grep 'READ COUNT: trimmed pair1' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $5}' >> tmp
	grep 'READ COUNT: trimmed pair2' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $5}' >> tmp
	grep 'INFO: READ COUNT: trimmed orphan1' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $5}' >> tmp
	grep 'INFO: READ COUNT: trimmed orphan2' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $5}' >> tmp
	grep 'paired_contam_1.fastq' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $3}' >> tmp
	grep 'paired_contam_2.fastq' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $3}' >> tmp
	grep 'unmatched_1_contam.fastq' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $3}' >> tmp
	grep 'unmatched_2_contam.fastq' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $3}' >> tmp
	grep 'INFO: READ COUNT: final pair1' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $5}' >> tmp
	grep 'INFO: READ COUNT: final pair2' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $5}' >> tmp
	grep 'INFO: READ COUNT: final orphan1 :' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $5}' >> tmp
	grep 'INFO: READ COUNT: final orphan2 :' ../${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log | awk -F ': ' '{print $5}' >> tmp
	less tmp | paste -s >> tmp2
	paste -d "\n" tmp2 >> 01_reads_stat.txt
	rm tmp*
done

module list

module purge

