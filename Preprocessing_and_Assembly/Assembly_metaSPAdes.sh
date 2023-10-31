#!/bin/bash
#SBATCH --job-name=metaSPAdes
#SBATCH --output=metaSPAdes_%A_%a.out
#SBATCH --mem=60gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=24
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

sample_dir=$1 #directory with the FASTQ files with clean reads
unmatched_sample_dir=$2 #directory with the FASTQ files of unmatched reads
sample_list=${sample_dir}/$3 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo -e '\n---- Copying files to tmpdir ----'

mkdir -p ${TMPDIR}/${SAMPLE_ID}/assembly_data/
cp ${sample_dir}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz ${TMPDIR}/${SAMPLE_ID}/assembly_data/
cp ${sample_dir}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz ${TMPDIR}/${SAMPLE_ID}/assembly_data/
cp ${unmatched_sample_dir}/${SAMPLE_ID}_kneaddata_unmatched.fastq.gz ${TMPDIR}/${SAMPLE_ID}/assembly_data/

echo -e '\n---- RUNNING metaSPAdes ON '${SAMPLE_ID}' SAMPLE ----'

mkdir -p ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/

# Clean environment, load modules 
module purge; module load SPAdes; module list

metaspades.py \
    -1 ${TMPDIR}/${SAMPLE_ID}/assembly_data/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 ${TMPDIR}/${SAMPLE_ID}/assembly_data/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    -s ${TMPDIR}/${SAMPLE_ID}/assembly_data/${SAMPLE_ID}_kneaddata_unmatched.fastq.gz \
    -o ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/ \
    -m $((${SLURM_MEM_PER_NODE} / 1024)) \
    --threads ${SLURM_CPUS_PER_TASK}

# Add sample name to output files
mv ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/contigs.fasta ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/${SAMPLE_ID}_metaspades_contigs.fa 
mv ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/scaffolds.fasta ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/${SAMPLE_ID}_metaspades_scaffolds.fa
mv ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/spades.log ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/${SAMPLE_ID}_spades.log

echo -e '\n---- Moving results to SCRATCH. Generating folders with CONTIGS, SCAFFOLDS and LOG FILES. Input for Viral Identification----'

mkdir -p metaSPAdes/CONTIGS metaSPAdes/SCAFFOLDS metaSPAdes/LOG_files 
mkdir -p metaSPAdes/CONTIGS/MD5 metaSPAdes/SCAFFOLDS/MD5 metaSPAdes/SUMMARY_RESULTS

if ! grep -Eq 'Error|FAILED' ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/${SAMPLE_ID}_spades.log; then
	echo -e "metaSPAdes assembly for sample ${SAMPLE_ID} completed successfully.\n" > metaSPAdes/SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt
  rsync -av $(find ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/ -name "${SAMPLE_ID}_spades.log" -type f) metaSPAdes/LOG_files
	rsync -av $(find ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS -name "${SAMPLE_ID}_metaspades_contigs.fa" -type f) metaSPAdes/CONTIGS
	rsync -av $(find ${TMPDIR}/${SAMPLE_ID}/assembly_data/MSP_RESULTS/ -name "${SAMPLE_ID}_metaspades_scaffolds.fa" -type f) metaSPAdes/SCAFFOLDS
else
	echo "MetaSPAdes STEP failed"
fi

echo "> Removing data from tmpdir"
rm -r ${TMPDIR}/${SAMPLE_ID}/assembly_data

echo -e '\n---- RUNNING QUAST ON '${SAMPLE_ID}' SAMPLE ----' # Quality Assessment

mkdir -p QUAST/
mkdir -p QUAST/${SAMPLE_ID}

# Clean environment, load modules 
module purge; module load QUAST; module list

if [ -f metaSPAdes/CONTIGS/${SAMPLE_ID}_metaspades_contigs.fa ]; then
	echo "> Assessing the assembly quality"
	quast.py \
        	metaSPAdes/CONTIGS/${SAMPLE_ID}_metaspades_contigs.fa \
        	-o QUAST/${SAMPLE_ID} \
        	-m $((${SLURM_MEM_PER_NODE} / 1024)) \
        	--threads ${SLURM_CPUS_PER_TASK}
	rm QUAST/${SAMPLE_ID}/*tex QUAST/${SAMPLE_ID}/*tsv
fi

echo -e '\n---- Generating md5sums ----'
md5sum metaSPAdes/CONTIGS/${SAMPLE_ID}_metaspades_contigs.fa > metaSPAdes/CONTIGS/MD5/${SAMPLE_ID}_MD5.txt
md5sum metaSPAdes/SCAFFOLDS/${SAMPLE_ID}_metaspades_scaffolds.fa > metaSPAdes/SCAFFOLDS/MD5/${SAMPLE_ID}_MD5.txt

echo -e '\n---- Assembly step DONE ----'
