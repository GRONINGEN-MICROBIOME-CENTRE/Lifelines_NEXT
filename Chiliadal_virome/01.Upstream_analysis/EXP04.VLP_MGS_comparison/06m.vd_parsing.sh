#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --output=./out/06.mpar/VD_PentaChiliadal_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

SAMPLE_DIR=$2

echo ${SAMPLE_LIST}
echo ${SAMPLE_DIR}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load seqtk/1.3-GCC-11.3.0

mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy

# --- EXTRACTING VIRSORTER2 ---
echo "> Parsing VirSorter2 results"

awk 'NR>1 {print $1}' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/VirSorter2/${SAMPLE_ID}_final-viral-score.tsv | \
	awk -F '\|' '{print $1}' | sort | uniq | awk '{print $0 "\tVirSorter2"}' > \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/vs2_tidy

# --- EXTRACTING DEEPVIRFINDER ---
echo "> Parsing DeepVirFinder results"

awk 'NR>1' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/DeepVirFinder/${SAMPLE_ID}_contigs.min1kbp.fasta_gt1000bp_dvfpred.txt | \
	awk '$3 >= 0.94' | awk -F '\t' '{print $1}' | sort | uniq | awk '{print $0 "\tDeepVirFinder"}' > \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/dvf_tidy

# --- EXTRACTING VIBRANT ---
echo "> Parsing VIBRANT results"

awk -F '_fragment' '{print $1}' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/VIBRANT/${SAMPLE_ID}_contigs.min1kbp.AA.phages_combined.txt | \
	sort | uniq | awk '{print $0 "\tVIBRANT"}' > \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/vib_tidy

# --- EXTRACTING GENOMAD ---
echo "> Parsing geNomad results"

awk 'NR>1 {print $1}' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/geNomad/${SAMPLE_ID}_contigs.min1kbp_virus_summary.tsv | \
	 sed 's/|provirus.*//' | sort | uniq | awk '{print $0 "\tgeNomad"}' > \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/gnd_tidy

# --- EXTRACTING CENOTE-TAKER3 ---
echo "> Parsing Cenote-Taker3"

awk 'NR>1 {print $2}' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/CenoteTaker3/${SAMPLE_ID}_CenoteTaker3_virus_summary.tsv | \
	sort | uniq | awk '{print $0 "\tCenoteTaker3"}' > \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/ct3_tidy

# --- COMBINING THE PUTATIVE VIRUSES ---
echo "> Combining putative virus contigs IDs"

cat ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/*_tidy | \
	awk -F '\t' '{print $1}' | \
	sort | uniq > \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/all_predicted_viral_ids

# --- CREATING PUTATIVE VIRUS FASTA ---
echo "> Pulling virus contigs to a separate fasta"
seqtk \
        subseq \
        -l60 \
	${SAMPLE_DIR}/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/all_predicted_viral_ids \
	> ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta

echo "> Check that none of the tools changed the contigIDs"
N_CONTIGS_IDS=$(wc -l < ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/all_predicted_viral_ids)
N_VIR_CONTIGS=$(grep '>' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta | wc -l)

if [ $N_CONTIGS_IDS -eq $N_VIR_CONTIGS ]; then
    echo "Number of contig IDs and contigs is equal"
else
    echo "Number of contig IDs and contigs is unequal"
fi

# --- CREATING TABLE OF ORIGIN ---
echo "> Creating long-format table of origin"

cat ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/*_tidy > ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_table_of_origin
cp ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_table_of_origin ../VIR_DB/MGS_table_of_origin/
rm ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/*_tidy

module purge
