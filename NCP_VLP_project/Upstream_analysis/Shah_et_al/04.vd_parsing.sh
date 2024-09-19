#!/bin/bash
#SBATCH --job-name=ViromeDiscovery
#SBATCH --error=./err/04.par/VD_AMG_%A_%a.err
#SBATCH --output=./out/04.par/VD_AMG_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST})

echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES --- 
module purge
module load seqtk/1.3-GCC-11.3.0

mkdir -p ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy

# --- EXTRACTING VIRSORTER2 ---
echo "> Parsing VirSorter2 results"

awk 'NR>1 {print $1}' ../SAMPLES/${SAMPLE_ID}/virome_discovery/VirSorter2/final-viral-score.tsv | \
	awk -F '\|' '{print $1}' | sort | uniq | awk '{print $0 "\tVirSorter2"}' > \
	../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/vs2_tidy

# --- EXTRACTING DEEPVIRFINDER ---
echo "> Parsing DeepVirFinder results"

awk 'NR>1' ../SAMPLES/${SAMPLE_ID}/virome_discovery/DeepVirFinder/${SAMPLE_ID}_contigs.min1kbp.fasta_gt1000bp_dvfpred.txt | \
	awk '$3 >= 0.94' | awk -F '\t' '{print $1}' | sort | uniq | awk '{print $0 "\tDeepVirFinder"}' > \
	../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/dvf_tidy

# --- EXTRACTING VIBRANT ---
echo "> Parsing VIBRANT results"

awk -F '_fragment' '{print $1}' ../SAMPLES/${SAMPLE_ID}/virome_discovery/VIBRANT/${SAMPLE_ID}_contigs.min1kbp.AA.phages_combined.txt | \
	sort | uniq | awk '{print $0 "\tVIBRANT"}' > \
	../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/vib_tidy

# --- EXTRACTING GENOMAD ---
echo "> Parsing geNomad results"

awk 'NR>1 {print $1}' ../SAMPLES/${SAMPLE_ID}/virome_discovery/geNomad/${SAMPLE_ID}_contigs.min1kbp_virus_summary.tsv | \
	 sed 's/|provirus.*//' | sort | uniq | awk '{print $0 "\tgeNomad"}' > \
	../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/gnd_tidy

# --- COMBINING THE PUTATIVE VIRUSES ---
echo "> Combining putative virus contigs IDs"

cat ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/*_tidy | \
	awk -F '\t' '{print $1}' | \
	sort | uniq > \
	../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/all_predicted_viral_ids

# --- CREATING PUTATIVE VIRUS FASTA ---
echo "> Pulling virus contigs to a separate fasta"
seqtk \
        subseq \
        -l60 \
	../SAMPLES/${SAMPLE_ID}/01_sc_assembly/${SAMPLE_ID}_contigs.min1kbp.fasta \
	../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/all_predicted_viral_ids \
	> ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta

echo "> Check that none of the tools changed the contigIDs"
N_CONTIGS_IDS=$(wc -l < ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/all_predicted_viral_ids)
N_VIR_CONTIGS=$(grep '>' ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta | wc -l)

if [ $N_CONTIGS_IDS -eq $N_VIR_CONTIGS ]; then
    echo "Number of contig IDs and contigs is equal"
else
    echo "Number of contig IDs and contigs is unequal"
fi

# --- CREATING TABLE OF ORIGIN ---
echo "> Creating long-format table of origin"

cat ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/*_tidy > ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_table_of_origin
cp ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_table_of_origin ../VIR_DB/table_of_origin/
rm ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/*_tidy

module purge
