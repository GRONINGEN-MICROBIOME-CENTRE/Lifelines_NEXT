#!/bin/bash
#SBATCH --job-name=Binning_taxonomy_metaWRAP
#SBATCH --output=Binning_taxonomy_metaWRAP_%A_%a.out
#SBATCH --mem=90gb
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --partition=regular

FASTQ_dir=$1 #directory with the FASTQ files with clean reads
ASSEMBLY_dir=$2 #directory with the assembly contigs 
sample_list=${FASTQ_dir}/$3 #file with the list of all samples in the directory
SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${sample_list})

echo -e '\n-------------------- WORKING WITH '${SAMPLE_ID}' SAMPLE --------------------'

echo -e '\n---- Copying files to tmpdir ----'

mkdir -p ${TMPDIR}/${SAMPLE_ID}/binning_data/
rsync -av ${FASTQ_dir}/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz ${TMPDIR}/${SAMPLE_ID}/binning_data/
rsync -av ${FASTQ_dir}/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz ${TMPDIR}/${SAMPLE_ID}/binning_data/
rsync -av ${ASSEMBLY_dir}/${SAMPLE_ID}_metaspades_contigs.fa ${TMPDIR}/${SAMPLE_ID}/binning_data/

echo -e '\n---- UNZIPPING FASTQs FROM '${SAMPLE_ID}' SAMPLE ----' #(metaWRAP does not work on compressed files)

gunzip ${TMPDIR}/${SAMPLE_ID}/binning_data/*gz

echo -e '\n---- RUNNING BINNING ON '${SAMPLE_ID}' SAMPLE ----'

mkdir -p ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-tifn/condas/conda_metawrap/; conda list 

metawrap binning \
    -o ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING \
    -t ${SLURM_CPUS_PER_TASK} \
    -a ${TMPDIR}/${SAMPLE_ID}/binning_data/${SAMPLE_ID}_metaspades_contigs.fa \
    -m ${SLURM_MEM_PER_NODE} \
    --metabat2 \
    --maxbin2 \
    --concoct ${TMPDIR}/${SAMPLE_ID}/binning_data/*fastq

# Modify names of relevant output files
mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/maxbin2_out/bin.log \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/${SAMPLE_ID}_maxbin2_binning.log

cat ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/concoct_out/args.txt \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/concoct_out/log.txt > \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/${SAMPLE_ID}_concoct_binning.log

echo -e '\n---- RUNNING BINNING REFINEMENT ON '${SAMPLE_ID}' SAMPLE ----'

metawrap bin_refinement \
    -o ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT \
    -t ${SLURM_CPUS_PER_TASK} \
    -m ${SLURM_MEM_PER_NODE} \
    -A ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/metabat2_bins/ \
    -B ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/maxbin2_bins/ \
    -C ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/concoct_bins/ \
    -c 50 \ # Minimum completeness
    -x 10 # Maximum contamination

# Modify names of relevant output files
mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/metawrap_50_10_bins.contigs \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins.contigs
mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/metawrap_50_10_bins.stats \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins.stats
mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/metawrap_50_10_bins \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins 

for file in ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins/*.fa; do
    base_name=$(basename $file)
    mv $file ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins/${SAMPLE_ID}_${base_name}
done

echo -e '\n---- RUNNING QUANTIFICATION OF BINS FROM '${SAMPLE_ID}' SAMPLE ----'

metawrap quant_bins \
    -b ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins  \
    -o ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/QUANT_BINS \
    -a ${TMPDIR}/${SAMPLE_ID}/binning_data/${SAMPLE_ID}_metaspades_contigs.fa ${TMPDIR}/${SAMPLE_ID}/binning_data/*fastq

# Modify names of relevant output files
cat ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/QUANT_BINS/assembly_index/indexing.log \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/QUANT_BINS/alignment_files/${SAMPLE_ID}_kneaddata_paired.quant/logs/salmon_quant.log > \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/${SAMPLE_ID}_quant_bins.log

mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/QUANT_BINS/quant_files/${SAMPLE_ID}_kneaddata_paired.quant.counts \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/QUANT_BINS/quant_files/${SAMPLE_ID}_quant_counts.txt
mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/QUANT_BINS/bin_abundance_table.tab \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/QUANT_BINS/${SAMPLE_ID}_bin_abundance_table.tab

echo -e '\n---- RUNNING TAXONOMIC ASSIGNMENT OF BINS FROM '${SAMPLE_ID}' SAMPLE ----'

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-tifn/condas/conda_GTDB_TK/; conda list 

GTDBTK_DATA_PATH=/scratch/hb-tifn/DBs/DB_GTDB/release214

gtdbtk classify_wf \
    --genome_dir ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins \
    --out_dir ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_CLASSIFICATION \
    --extension .fa \
    --skip_ani_screen \
    --cpus ${SLURM_CPUS_PER_TASK}

# Modify names of relevant output files
mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_CLASSIFICATION/gtdbtk.log \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/${SAMPLE_ID}_gtdbtk_tax.log
mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_CLASSIFICATION/classify/gtdbtk.bac120.summary.tsv \
${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_CLASSIFICATION/classify/${SAMPLE_ID}_gtdbtk.bac120.summary.tsv

echo -e '\n---- Moving results to SCRATCH. Generating folders with BINS, QUANTIFICATION, TAXONOMY and LOG files ----'

mkdir -p metaWRAP/BINNING metaWRAP/QUANTIFICATION metaWRAP/TAXONOMY metaWRAP/LOG_files metaWRAP/LOG_files/${SAMPLE_ID}
mkdir -p metaWRAP/SUMMARY_RESULTS metaWRAP/OUTPUT_files

LOG_DIR="${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS"
BINNING_DIR="${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT"
QUANTIFICATION_DIR="${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/QUANT_BINS"
TAXONOMY_DIR="${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_CLASSIFICATION/classify"

if ! grep -qi 'Error' ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/*.log; then
	echo -e "metaWRAP binning and taxonomy assignment for sample ${SAMPLE_ID} completed successfully.\n" > metaWRAP/SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt
  mkdir -p metaWRAP/BINNING/${SAMPLE_ID} 
  rsync -av $(find ${LOG_DIR} -name "${SAMPLE_ID}_*.log" -type f) metaWRAP/LOG_files; rm -r metaWRAP/LOG_files/${SAMPLE_ID}
  rsync -av $(find ${BINNING_DIR} -name "${SAMPLE_ID}_metawrap_50_10_bins" -type d) metaWRAP/BINNING/${SAMPLE_ID}
  rsync -av $(find ${BINNING_DIR} -name "${SAMPLE_ID}_metawrap_50_10_bins.contigs" -type f) metaWRAP/BINNING/${SAMPLE_ID}
  rsync -av $(find ${BINNING_DIR} -name "${SAMPLE_ID}_metawrap_50_10_bins.stats" -type f) metaWRAP/BINNING/${SAMPLE_ID}
	rsync -av $(find ${QUANTIFICATION_DIR}/quant_files/ -name "${SAMPLE_ID}_quant_counts.txt" -type f) metaWRAP/QUANTIFICATION
	rsync -av $(find ${QUANTIFICATION_DIR} -name "${SAMPLE_ID}_bin_abundance_table.tab" -type f) metaWRAP/QUANTIFICATION
	rsync -av $(find ${TAXONOMY_DIR} -name "${SAMPLE_ID}_gtdbtk.bac120.summary.tsv" -type f) metaWRAP/TAXONOMY
else
	echo "metaWRAP STEP for sample ${SAMPLE_ID} failed. Check the LOG files for more information." > metaWRAP/SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt
  rsync -av $(find ${LOG_DIR} -name "${SAMPLE_ID}_*.log" -type f) metaWRAP/LOG_files; rm -r metaWRAP/LOG_files/${SAMPLE_ID}
fi

echo "> Removing data from tmpdir"
rm -r ${TMPDIR}/${SAMPLE_ID}/binning_data

echo -e "\n---- Generating md5sums ---"
md5sum metaWRAP/BINNING/${SAMPLE_ID}/${SAMPLE_ID}_metawrap_50_10_bins/${SAMPLE_ID}*.fa > metaWRAP/BINNING/${SAMPLE_ID}/${SAMPLE_ID}_metawrap_50_10_bins/${SAMPLE_ID}_MD5.txt

echo -e "\n---- Binning step DONE ----"
