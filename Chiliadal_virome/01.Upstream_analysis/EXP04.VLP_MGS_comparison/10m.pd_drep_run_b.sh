#!/bin/bash
#SBATCH --job-name=PD_deRep
#SBATCH --output=./out/09.mdrp/PD_dRep_b.out
#SBATCH --mem=64gb
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

# Concatenate per-fragment results of all vs all blast
cat ../VIR_DB/VLP_MGS_DREP/OUTPUT/NEXT_frag*_blast.tsv > ../VIR_DB/VLP_MGS_DREP/NEXT_viruses_blast.tsv 

# --- LOAD MODULES --- 
module purge
module load Python/3.10.8-GCCcore-12.2.0
module load CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8
module list

# Next, calculate pairwise ANI by combining local alignments between sequence pairs:
python /scratch/p282752/tools/checkv_scripts/anicalc.py \
	-i ../VIR_DB/VLP_MGS_DREP/NEXT_viruses_blast.tsv \
	-o ../VIR_DB/VLP_MGS_DREP/NEXT_viruses_ani.tsv

# anicalc.py is available at https://bitbucket.org/berkeleylab/checkv/src/master/scripts/anicalc.py

# Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
python /scratch/p282752/tools/checkv_scripts/aniclust.py \
    --fna ../VIR_DB/VLP_MGS_DREP/QualFilt_virus_contigs.fasta \
    --ani ../VIR_DB/VLP_MGS_DREP/NEXT_viruses_ani.tsv \
    --out ../VIR_DB/VLP_MGS_DREP/NEXT_viral_clusters_MGS_VLP.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0

# aniclust.py is available at https://bitbucket.org/berkeleylab/checkv/src/master/scripts/aniclust.py

awk -F '\t' '{print $1}' ../VIR_DB/VLP_MGS_DREP/NEXT_viral_clusters_MGS_VLP.tsv > ../VIR_DB/VLP_MGS_DREP/vOTUr_IDs

# --- LOAD MODULES ---
module purge
module load seqtk/1.3-GCC-11.3.0

# Creating a fasta-file with all vOTU representative seqeunces:
seqtk \
        subseq \
        -l60 \
        ../VIR_DB/VLP_MGS_DREP/QualFilt_virus_contigs.fasta \
	../VIR_DB/VLP_MGS_DREP/vOTUr_IDs \
        > ../VIR_DB/VLP_MGS_DREP/NEXT_vOTUr_MGS_VLP.fasta

# --- LOAD MODULES ---
module purge
module load R/4.3.2-gfbf-2023a

# Getting VC clustering info
Rscript dereplication_stat.R ../VIR_DB/VLP_MGS_DREP/NEXT_viral_clusters_MGS_VLP.tsv

module purge
