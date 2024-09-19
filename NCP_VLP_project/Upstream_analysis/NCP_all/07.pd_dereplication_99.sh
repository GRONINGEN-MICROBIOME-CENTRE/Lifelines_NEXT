#!/bin/bash
#SBATCH --job-name=PostDiscovery_deRep_NCP_99
#SBATCH --error=./err/07.drp/PD_initial_dereplication_NCP_99.err
#SBATCH --output=./out/07.drp/PD_initial_dereplication_NCP_99.out
#SBATCH --mem=64gb
#SBATCH --time=80:00:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES --- 
module purge
module load Python/3.10.8-GCCcore-12.2.0
module load CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8
module list

# Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (99% ANI + 85% AF):
python aniclust.py \
    --fna ../VIR_DB/virus_contigs/all_extended_pruned_viral_renamed_NCP.fasta \
    --ani ../VIR_DB/initial_dereplication/NCP_viruses_ani.tsv \
    --out ../VIR_DB/initial_dereplication/NCP_viral_clusters_99.tsv \
    --min_ani 99 \
    --min_tcov 85 \
    --min_qcov 0

# aniclust.py is available at https://bitbucket.org/berkeleylab/checkv/src/master/scripts/aniclust.py

# --- LOAD MODULES --- 
module purge
module load seqtk/1.3-GCC-11.3.0

# Creating a fasta-file with vOTU representatives:
awk -F '\t' '{print $1}' \
	../VIR_DB/initial_dereplication/NCP_viral_clusters_99.tsv \
	> ../VIR_DB/initial_dereplication/NCP_vOTU_representatives_99

seqtk \
        subseq \
        -l60 \
        ../VIR_DB/virus_contigs/all_extended_pruned_viral_renamed_NCP.fasta \
        ../VIR_DB/initial_dereplication/NCP_vOTU_representatives_99 \
        > ../VIR_DB/virus_contigs/NCP_vOTU_representatives_w_neg_der99.fasta

module purge
