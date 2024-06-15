#!/bin/bash
#SBATCH --job-name=PostDiscovery
#SBATCH --output=./out/09.dbs/DB_launch_%A_%a.out
#SBATCH --mem=8gb
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

contig_file=$1 #path to FASTA file with the predicted viral contigs

DB_DIR=$(dirname ${contig_file})

DB=$(basename ${DB_DIR})

echo "Database: ${DB}"

# --- LOAD MODULES ---
module purge
module load QUAST
module list

# --- GETTING DB PHENOS ---
quast.py \
	${contig_file} \
        -o ${DB_DIR}/quast \
        -m $((${SLURM_MEM_PER_NODE} / 1024)) \
        --threads ${SLURM_CPUS_PER_TASK}

# --- KEEPING ONLY NECESSARY FILES ---
rm -r ${DB_DIR}/quast/basic_stats
rm -r ${DB_DIR}/quast/icarus*
rm ${DB_DIR}/quast/report.html
rm ${DB_DIR}/quast/report.pdf
rm ${DB_DIR}/quast/*.tex
rm ${DB_DIR}/quast/*.txt

## --- ESTIMATING DB SIZE ---
if [ $(grep '>' ${contig_file} | wc -l) -ge 5000 ]; then
	echo "Database size exceeds 5,000 sequences"
	# It is surprising to learn 7 years into bash scripting, that there is no internal function for number rounding in shell
	SPLIT_FACTOR=$(echo "scale=0; ( $(grep '>' "${contig_file}" | wc -l) + 2500) / 5000" | bc)
else
	SPLIT_FACTOR=1
fi 

# --- RUNNING CHECKV ---
if [ ${SPLIT_FACTOR} -eq 1 ]; then 
	echo "Running DB QC as is"
	echo "${contig_file}" > /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/01.EXTERNAL_DBs/${DB}
	sbatch --array=1-1 09.db_qc.sh /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/01.EXTERNAL_DBs/${DB}
else 
	mkdir -p ${DB_DIR}/fragmented
	cd ${DB_DIR}/fragmented
	echo "Splitting the database ${DB} in ${SPLIT_FACTOR} pieces"
	/scratch/p282752/ANALYSIS_CHILIADAL/scripts/fastasplitn ${contig_file} ${SPLIT_FACTOR}
	ls *.fa > split_list
	
	for i in $(cat split_list); do 
		i=${i%.fa}
		mkdir ${DB_DIR}/${DB}_${i}
		mv ${DB_DIR}/fragmented/${i}.fa ${DB_DIR}/${DB}_${i}/${DB}_${i}.fasta
		ls ${DB_DIR}/${DB}_${i}/${DB}_${i}.fasta >> /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/01.EXTERNAL_DBs/${DB} 
	done
	
	rm -r ${DB_DIR}/fragmented
	cd /scratch/${USER}/ANALYSIS_CHILIADAL/scripts
	head -n 1 ../VIR_DB/table_of_origin/Extended_table_of_origin > ${DB_DIR}/Extended_TOF
	sbatch --array=1-${SPLIT_FACTOR} under_dev_db_prep.sh /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/01.EXTERNAL_DBs/${DB}
fi

module purge
