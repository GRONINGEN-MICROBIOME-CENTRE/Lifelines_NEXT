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
    --concoct \
    --metabat2 \
    --maxbin2 ${TMPDIR}/${SAMPLE_ID}/binning_data/*fastq

# Check if maxbin2 log exists, then move it
if [ -f "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/maxbin2_out/bin.log" ]; then
    mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/maxbin2_out/bin.log \
    ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/${SAMPLE_ID}_maxbin2_binning.log
fi

# Check if concoct output files exist, then concatenate and move them
if [ -f "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/concoct_out/args.txt" ] && \
   [ -f "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/concoct_out/log.txt" ]; then
    cat ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/concoct_out/args.txt \
    ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/work_files/concoct_out/log.txt > \
    ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/${SAMPLE_ID}_concoct_binning.log
fi

echo -e '\n---- RUNNING BINNING REFINEMENT ON '${SAMPLE_ID}' SAMPLE ----'

# Initialize an empty string for options
BINNING_OPTIONS=""

# Check if metabat2 bins exist (bin*fa files present)
if [ -d "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/metabat2_bins" ] && \
   [ "$(ls -1 ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/metabat2_bins/bin*fa 2>/dev/null)" ]; then
    BINNING_OPTIONS="${BINNING_OPTIONS} -A ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/metabat2_bins/"
fi

# Check if maxbin2 bins exist
if [ -d "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/maxbin2_bins" ] && \
   [ "$(ls -1 ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/maxbin2_bins/bin*fa 2>/dev/null)" ]; then
    BINNING_OPTIONS="${BINNING_OPTIONS} -B ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/maxbin2_bins/"
fi

# Check if concoct bins exist
if [ -d "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/concoct_bins" ] && \
   [ "$(ls -1 ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/concoct_bins/bin*fa 2>/dev/null)" ]; then
    BINNING_OPTIONS="${BINNING_OPTIONS} -C ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/INITIAL_BINNING/concoct_bins/"
fi

# Only run the bin_refinement if there are any available bins
if [ -n "$BINNING_OPTIONS" ]; then
    metawrap bin_refinement \
    -o ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT \
    -t ${SLURM_CPUS_PER_TASK} \
    -m ${SLURM_MEM_PER_NODE} \
    $BINNING_OPTIONS \
    -c 50 \
    -x 10
        
    # Check if metawrap_50_10_bins directory exists
    if [ ! -d "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/metawrap_50_10_bins" ]; then
        echo "No bins with >50% completeness and <10% contamination were found."
    else
        # Modify names of relevant output files
        mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/metawrap_50_10_bins.contigs \
        ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins.contigs
        
        mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/metawrap_50_10_bins.stats \
        ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins.stats
        
        mv ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/metawrap_50_10_bins \
        ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins

        # Rename all .fa files in the bins directory
        for file in ${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins/*.fa; do
            base_name=$(basename "$file")
            mv "$file" "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/${SAMPLE_ID}_metawrap_50_10_bins/${SAMPLE_ID}_$base_name"
        done

	echo -e '\n---- SKIPPING QUANTIFICATION OF THE BIN SET FROM '${SAMPLE_ID}' SAMPLE ----'

    	# No quantification will be done with quant_bins (abundance will be etsimated later with CoverM)
     
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
 
    fi    
else
    echo "No bins were generated from any method (metaBAT, maxBin2 or CONCOCT) for sample ${SAMPLE_ID}."
fi


echo -e '\n---- Moving results to SCRATCH. Generating folders with BINS, TAXONOMY and LOG files ----'

mkdir -p metaWRAP/BINNING metaWRAP/TAXONOMY metaWRAP/LOG_files/ metaWRAP/SUMMARY_RESULTS metaWRAP/OUTPUT_files

LOG_DIR="${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS"
BINNING_DIR="${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT"
TAXONOMY_DIR="${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_CLASSIFICATION/classify"

SUMMARY_FILE="metaWRAP/SUMMARY_RESULTS/${SAMPLE_ID}_summary.txt"

# Initialize summary for the sample
echo "Sample ${SAMPLE_ID} Pipeline Summary:" > ${SUMMARY_FILE}

# Check success of binning
if [ -n "$BINNING_OPTIONS" ]; then
    echo "- Binning completed successfully." >> ${SUMMARY_FILE}
else
    echo "- Binning failed or no bins were generated." >> ${SUMMARY_FILE}
fi

# Check success of refinement
REFINEMENT_STATUS="failed"
if [ -d "${TMPDIR}/${SAMPLE_ID}/binning_data/metaWRAP_RESULTS/BIN_REFINEMENT/metawrap_50_10_bins" ]; then
    REFINEMENT_STATUS="completed"
    mkdir -p metaWRAP/BINNING/${SAMPLE_ID}
    rsync -av $(find ${BINNING_DIR} -name "${SAMPLE_ID}_metawrap_50_10_bins" -type d) metaWRAP/BINNING/${SAMPLE_ID}
    rsync -av $(find ${BINNING_DIR} -name "${SAMPLE_ID}_metawrap_50_10_bins.contigs" -type f) metaWRAP/BINNING/${SAMPLE_ID}
    rsync -av $(find ${BINNING_DIR} -name "${SAMPLE_ID}_metawrap_50_10_bins.stats" -type f) metaWRAP/BINNING/${SAMPLE_ID}

    echo -e "\n---- Generating md5sums ---"
    md5sum metaWRAP/BINNING/${SAMPLE_ID}/${SAMPLE_ID}_metawrap_50_10_bins/${SAMPLE_ID}*.fa > metaWRAP/BINNING/${SAMPLE_ID}/${SAMPLE_ID}_metawrap_50_10_bins/${SAMPLE_ID}_MD5.txt
    
    echo "- Refinement completed successfully." >> ${SUMMARY_FILE}
else
    echo "- Refinement failed or no refined bins passed thresholds." >> ${SUMMARY_FILE}
fi

# Check success of taxonomy assignment
TAXONOMY_STATUS="failed"
if [ -f "${TAXONOMY_DIR}/${SAMPLE_ID}_gtdbtk.bac120.summary.tsv" ]; then
    TAXONOMY_STATUS="completed"
    echo "- Taxonomy assignment completed successfully." >> ${SUMMARY_FILE}
    rsync -av $(find ${TAXONOMY_DIR} -name "${SAMPLE_ID}_gtdbtk.bac120.summary.tsv" -type f) metaWRAP/TAXONOMY
else
    echo "- Taxonomy assignment failed or was skipped." >> ${SUMMARY_FILE}
fi

# Transfer log files
rsync -av $(find ${LOG_DIR} -name "${SAMPLE_ID}_*.log" -type f) metaWRAP/LOG_files

# Summarize overall results
echo -e "\n---- Overall Pipeline Summary for Sample ${SAMPLE_ID} ----" >> ${SUMMARY_FILE}
if [[ "$BINNING_OPTIONS" = "" || "$REFINEMENT_STATUS" = "failed" || "$TAXONOMY_STATUS" = "failed" ]]; then
    echo -e "The pipeline encountered issues during execution. Check log files for details.\n" >> ${SUMMARY_FILE}
    echo -e "Details of failures are as follows:" >> ${SUMMARY_FILE}
    
    [[ "$BINNING_OPTIONS" = "" ]] && echo "- Note: Binning process either failed or did not produce any bins." >> ${SUMMARY_FILE}
    [[ "$REFINEMENT_STATUS" = "failed" ]] && echo "- Note: No bins passed the refinement criteria." >> ${SUMMARY_FILE}
    [[ "$TAXONOMY_STATUS" = "failed" ]] && echo "- Note: Taxonomy assignment was unsuccessful." >> ${SUMMARY_FILE}
else
    echo -e "The pipeline completed successfully with no issues detected.\n" >> ${SUMMARY_FILE}
fi

echo "> Removing data from tmpdir"
rm -r ${TMPDIR}/${SAMPLE_ID}/binning_data

echo -e "\n---- Binning step DONE ----"
