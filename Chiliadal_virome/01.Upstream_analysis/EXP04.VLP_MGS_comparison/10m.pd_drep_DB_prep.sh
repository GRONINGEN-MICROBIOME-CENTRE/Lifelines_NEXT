#!/bin/bash
#SBATCH --job-name=NEXT_n_NC_PD_deRep
#SBATCH --output=./out/09.mdrp/PD_NEXT_deRep.out
#SBATCH --mem=64gb
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

mkdir -p ../VIR_DB/VLP_MGS_DREP

# --- CONCATENATING ALL PUTATIVE VIRUS CONTIGS & THEIR ETOF---
head -n 1 ../VIR_DB/table_of_origin/Extended_table_of_origin > ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## from all unique fecal VLP samples:
for SAMPLE_ID in $(cat ../VLP_to_compare); do
	cat ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_pruned_viral_renamed.fasta >> \
		../VIR_DB/VLP_MGS_DREP/all_virus.fasta 
	awk 'NR>1' ../SAMPLES/${SAMPLE_ID}/virome_discovery/tidy/Extended_TOF >> \
		../VIR_DB/VLP_MGS_DREP/Extended_TOF
done

## from all unique fecal MGS samples:
for SAMPLE_ID in $(cat ../MGS_to_compare); do
        cat ../MGS/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_extended_pruned_viral_renamed.fasta >> \
                ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
        awk 'NR>1' ../MGS/${SAMPLE_ID}/virome_discovery/tidy/Extended_TOF >> \
                ../VIR_DB/VLP_MGS_DREP/Extended_TOF
done

## from all selected external databases:
## Guerin et al. https://doi.org/10.1016/j.chom.2018.10.002
cat /scratch/p282752/databases/Guerin_crAss/Guerin_crAss_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/Guerin_crAss/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## Yutin et al. https://doi.org/10.1038/s41467-021-21350-w
cat /scratch/p282752/databases/Yutin_crAss/Yutin_crAss_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/Yutin_crAss/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## NCBI_virus (CrAssvirales only) https://doi.org/10.1093/database/baaa062
cat /scratch/p282752/databases/NCBIVirus_crAss/NCBIVirus_crAss_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/NCBIVirus_crAss/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## Gulyaeva et al. https://doi.org/10.1016/j.celrep.2021.110204
cat /scratch/p282752/databases/Gulyaeva_crAss/Gulyaeva_crAss_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/Gulyaeva_crAss/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## Benler et al. https://doi.org/10.1186/s40168-021-01017-w
cat /scratch/p282752/databases/Benler_et_al/Benler_et_al_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/Benler_et_al/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## Shah et al. https://doi.org/10.1038/s41564-023-01345-7
cat /scratch/p282752/databases/COPSAC_virome/COPSAC_virome_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/COPSAC_virome/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## Gregory et al. https://doi.org/10.1016/j.chom.2020.08.003
cat /scratch/p282752/databases/GVD/GVD_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/GVD/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## Camarillo-Guerrero et al. https://doi.org/10.1016/j.cell.2021.01.029
cat /scratch/p282752/databases/GPD/GPD_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/GPD/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## Nayfach et al. https://doi.org/10.1038/s41564-021-00928-6
cat /scratch/p282752/databases/MGV/MGV_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/MGV/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## IMG/VR4 https://doi.org/10.1093/nar/gkac1037
cat /scratch/p282752/databases/IMGVR_13012024/IMGVR_13012024_pruned_renamed.fasta >> ../VIR_DB/VLP_MGS_DREP/all_virus.fasta
awk 'NR>1' /scratch/p282752/databases/IMGVR_13012024/Extended_TOF >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF

## FILTERING VIRUS CONTIGS ACCORDING TO LENGTH, N VIRAL GENES, PLASMID STATUS AND K-MER FREQ

## --- LOAD MODULES ---
module load R
Rscript Filter_virus_contigs_to_derep.R /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/VLP_MGS_DREP/Extended_TOF 

## --- LOAD MODULES ---
module purge
module load seqtk/1.3-GCC-11.3.0

## Creating a fasta-file with all filtered virus seqeunces:
seqtk \
        subseq \
        -l60 \
        ../VIR_DB/VLP_MGS_DREP/all_virus.fasta \
        ../VIR_DB/VLP_MGS_DREP/QualFilt_virus_contigs_IDs \
        > ../VIR_DB/VLP_MGS_DREP/QualFilt_virus_contigs.fasta

## Viral_RefSeq (release 223)
cat /scratch/p282752/databases/viral_refseq_apr_24/Viral_Refseq_all.fasta >> ../VIR_DB/VLP_MGS_DREP/QualFilt_virus_contigs.fasta
awk 'NR>1' /scratch/p282752/databases/viral_refseq_apr_24/Extended_TOF_simulated >> ../VIR_DB/VLP_MGS_DREP/Extended_TOF_filtered

## --- LOAD MODULES ---
module purge
module load BLAST+/2.13.0-gompi-2022a
module list

# --- DEREPLICATION ACCORDING TO MIUViG GUIDELINES ---

## First, create a blast+ database:
makeblastdb \
    -in ../VIR_DB/VLP_MGS_DREP/QualFilt_virus_contigs.fasta \
    -dbtype nucl \
    -out ../VIR_DB/VLP_MGS_DREP/NEXT_VIR_DB

## SPLITTING THE DB
cd ../VIR_DB/VLP_MGS_DREP/
/scratch/p282752/ANALYSIS_CHILIADAL/scripts/fastasplitn QualFilt_virus_contigs.fasta 1000
ls *.fa | sed 's/\.fa//g'> split_list

mkdir OUTPUT

cd /scratch/p282752/ANALYSIS_CHILIADAL/scripts
## Running the all vs all blast per fragment
sbatch --array=1-1000 10m.pd_drep_run_a.sh ../VIR_DB/VLP_MGS_DREP/split_list

module purge
