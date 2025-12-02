#!/bin/bash
#SBATCH --job-name=VirusDiscovery
#SBATCH --output=./out/08.mpru/PP_PentaChiliadal_%A_%a.out
#SBATCH --mem=32gb
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --open-mode=truncate

SAMPLE_LIST=$1

SAMPLE_DIR=$2

echo ${SAMPLE_LIST}
echo ${SAMPLE_DIR}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST} | cut -d "_" -f1)

echo "SAMPLE_ID=${SAMPLE_ID}"

mkdir -p ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning

# --- GETTING Post-COBRA contig lengths ---

# --- LOAD MODULES ---
module load bioawk
module list

bioawk -c fastx '{ print $name, length($seq) }' \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta \
	> ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/POST_CBR_LENGTH

module purge

# PREPARING TMP
mkdir -p ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning

# --- LOAD MODULES ---
module load ARAGORN/1.2.41-foss-2021b
module load Python/3.9.5-GCCcore-10.3.0
source /scratch/p282752/tools/python_envs/geNomad/bin/activate
module list

# --- RUNNING geNomad for prepruning ---
echo "> Running geNomad"

genomad \
        end-to-end \
        --enable-score-calibration \
        --cleanup \
        ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta \
        ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning \
	/scratch/p282752/databases/genomad_db

genomad --version

deactivate

# 1. Get putative virus contigs that have DTR or ITR:
grep -E "(DTR|ITR)" \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_viral_summary/${SAMPLE_ID}_extended_viral_virus_summary.tsv | \
	awk -F '\t' '{print $1}' \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_DTR_ITR_vIDs

# 2. Get the provirus sequences:
## We have to take it from the sample virus summary, because some proviruses might be discarded by the post-classification filtering process
grep "Provirus" \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_viral_summary/${SAMPLE_ID}_extended_viral_virus_summary.tsv | \
	awk -F '\t' '{print $1}' | \
	sort | \
	uniq \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_provirus_vIDs

awk -F '\|' '{print $1}' \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_provirus_vIDs | \
	sort | \
	uniq \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_provirus_source_vIDs

# 3. Get the rest of putative virus contigs (including those not recognized by geNomad as viral)
grep '>' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta | \
	sed 's/>//g' | \
	sort \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_viral_IDs

cat ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_DTR_ITR_vIDs \
       ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_provirus_source_vIDs | \
       sort | \
       uniq \
       > ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_DTR_ITR_provirus_vIDs

comm -23 \
	<(sort ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_viral_IDs) \
	<(sort ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_DTR_ITR_provirus_vIDs) \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_rest_vIDs

# --- LOAD MODULES ---
module purge
module load seqtk/1.3-GCC-11.3.0
module list

# 4. Create fastas
seqtk \
        subseq \
        -l60 \
	${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_DTR_ITR_vIDs \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_DTR_ITR.fasta

seqtk \
        subseq \
        -l60 \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_viral_find_proviruses/${SAMPLE_ID}_extended_viral_provirus.fna \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_provirus_vIDs \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_provirus.fasta

seqtk \
        subseq \
        -l60 \
        ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_rest_vIDs \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_rest.fasta

cat ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_provirus.fasta \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_rest.fasta \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/POST_GND.fasta

N_VIRUSES_PRE=$(grep '>' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_viral.fasta | wc -l)
N_VIRUSES_POST=$(grep '>' ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/POST_GND.fasta | sed 's/|.*//' | sort | uniq | wc -l)
N_VIRUSES_COMPLETE=$(cat ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_DTR_ITR_vIDs | wc -l)

if [ ${N_VIRUSES_PRE} -eq $((N_VIRUSES_POST + N_VIRUSES_COMPLETE)) ]; then
    echo "All putative viruses were processed"
else
    echo "N viruses in input and intermediate output is different"
fi

# --- RUNNING CheckV for pruning ---

# --- LOAD MODULES --- 
module purge
module load CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8
module list

# 1. Running CheckV
checkv end_to_end \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/POST_GND.fasta \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/CHV \
	-t ${SLURM_CPUS_PER_TASK} \
    	-d /scratch/hb-llnext/databases/checkv-db-v1.5

# --- CONCATENATING VIRUSES WITH PRUNED PROVIRUSES ---
sed -i 's/\ /_/g' ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/CHV/proviruses.fna
sed -i 's/\//_/g' ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/CHV/proviruses.fna

# --- RUNNING CheckV for quality assessment ---

# 1. Concatenating complete genomes & CheckV-pruned seqeunces & CheckV-untouched sequences
cat ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/geNomad_DTR_ITR.fasta \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/CHV/proviruses.fna \
	${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/CHV/viruses.fna \
	> ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_pruned_viral.fasta

# 2. Running CheckV
checkv end_to_end \
        ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_pruned_viral.fasta \
        ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/CHV_new \
        -t ${SLURM_CPUS_PER_TASK} \
        -d /scratch/hb-llnext/databases/checkv-db-v1.5

# --- COPYING the output ---
cp ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_pruned_viral.fasta ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning
cp ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_viral_summary/${SAMPLE_ID}_extended_viral_plasmid_summary.tsv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning
cp ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/${SAMPLE_ID}_extended_viral_summary/${SAMPLE_ID}_extended_viral_virus_summary.tsv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning
cp ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/CHV/contamination.tsv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning
cp ${TMPDIR}/${SAMPLE_ID}/Prophage_pruning/CHV_new/quality_summary.tsv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning

# --- CREATE NEW VIRUS CONTIGS METADATA AND IDs ---
module load R/4.3.2-gfbf-2023a
Rscript M_New_contigs_ID_and_metadata.R ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/ ${SAMPLE_ID}

# --- RENAME CONTIGS ---
cp ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/${SAMPLE_ID}_extended_pruned_viral.fasta ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/${SAMPLE_ID}_extended_pruned_viral_renamed.fasta
awk 'NR>1' ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/Extended_TOF | awk -F '\t' '{print $35"\t"$1}' | while IFS=$'\t' read -r old_id new_id; do
	sed "s/>$old_id\b/>$new_id/" -i ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/${SAMPLE_ID}_extended_pruned_viral_renamed.fasta
done

# --- CLEANING output ---
rm ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/${SAMPLE_ID}_extended_viral_plasmid_summary.tsv
rm ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/${SAMPLE_ID}_extended_viral_virus_summary.tsv
rm ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/contamination.tsv
rm ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/POST_CBR_LENGTH
rm ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/quality_summary.tsv
mv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/*.fasta ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/
mv ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning/Extended_TOF ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/
rm -r ${SAMPLE_DIR}/${SAMPLE_ID}/virome_discovery/tidy/Prophage_pruning

module purge
