#!/bin/bash
#SBATCH --job-name=PostDiscovery_deRep_NCP
#SBATCH --error=./err/07.drp/PD_dereplication_decontam_99_NCP.err
#SBATCH --output=./out/07.drp/PD_dereplication_decontam_99_NCP.out
#SBATCH --mem=64gb
#SBATCH --time=80:00:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- LOAD MODULES ---
module purge
module load BLAST+/2.13.0-gompi-2022a 
module list

# --- DEREPLICATION ACCORDING TO MIUViG GUIDELINES ---
mkdir -p ../VIR_DB/dereplication_decontam_99

# First, create a blast+ database:
makeblastdb \
    -in ../VIR_DB/decontamination/extended_pruned_viral_renamed_NCP_wo_nc99.fasta \
    -dbtype nucl \
    -out ../VIR_DB/dereplication_decontam_99/NCP_VIR_DB

# Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn \
    -query ../VIR_DB/decontamination/extended_pruned_viral_renamed_NCP_wo_nc99.fasta \
    -db ../VIR_DB/dereplication_decontam_99/NCP_VIR_DB \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -out ../VIR_DB/dereplication_decontam_99/NCP_viruses_blast_wo_nc99.tsv \
    -num_threads ${SLURM_CPUS_PER_TASK}

echo "all-vs-all blastn done!"

# --- LOAD MODULES --- 
module purge
module load Python/3.10.8-GCCcore-12.2.0
module load CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8
module list

# Next, calculate pairwise ANI by combining local alignments between sequence pairs:
python anicalc.py \
	-i ../VIR_DB/dereplication_decontam_99/NCP_viruses_blast_wo_nc99.tsv \
	-o ../VIR_DB/dereplication_decontam_99/NCP_viruses_ani_wo_nc99.tsv

# anicalc.py is available at https://bitbucket.org/berkeleylab/checkv/src/master/scripts/anicalc.py

# Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
python aniclust.py \
    --fna ../VIR_DB/decontamination/extended_pruned_viral_renamed_NCP_wo_nc99.fasta \
    --ani ../VIR_DB/dereplication_decontam_99/NCP_viruses_ani_wo_nc99.tsv \
    --out ../VIR_DB/dereplication_decontam_99/NCP_viral_clusters_wo_nc99.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0

# aniclust.py is available at https://bitbucket.org/berkeleylab/checkv/src/master/scripts/aniclust.py

# --- LOAD MODULES --- 
module purge
module load seqtk/1.3-GCC-11.3.0

# Creating a fasta-file with vOTU representatives:
awk -F '\t' '{print $1}' \
	../VIR_DB/dereplication_decontam_99/NCP_viral_clusters_wo_nc99.tsv \
	> ../VIR_DB/dereplication_decontam_99/NCP_vOTU_representatives_wo_nc99

seqtk \
        subseq \
        -l60 \
        ../VIR_DB/decontamination/extended_pruned_viral_renamed_NCP_wo_nc99.fasta \
        ../VIR_DB/dereplication_decontam_99/NCP_vOTU_representatives_wo_nc99 \
        > ../VIR_DB/virus_contigs/NCP_vOTU_representatives_der95_wo_nc99.fasta

module purge
