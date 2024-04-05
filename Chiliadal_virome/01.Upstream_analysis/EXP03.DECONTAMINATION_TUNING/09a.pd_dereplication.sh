#!/bin/bash
#SBATCH --job-name=PostDiscovery_deRep
#SBATCH --error=./err/09.drp/PD_w_UniVec_dereplication.err
#SBATCH --output=./out/09.drp/PD_w_UniVec_dereplication.out
#SBATCH --mem=64gb
#SBATCH --time=80:00:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# --- EXCLUDING VIRUS CONTIGS OF NGCTRLs---
grep 'NEXT' ../VIR_DB/virus_contigs/all_extended_pruned_viral_renamed.fasta | \
	grep -vE '(V1999|V2000)' | \
	awk -F '>' '{print $2}' > \
	../VIR_DB/virus_contigs/no_neg_extended_pruned_viral_renamed

# --- LOAD MODULES ---
module purge
module load seqtk/1.3-GCC-11.3.0

seqtk \
        subseq \
        -l60 \
        ../VIR_DB/virus_contigs/all_extended_pruned_viral_renamed.fasta \
        ../VIR_DB/virus_contigs/no_neg_extended_pruned_viral_renamed \
        > ../VIR_DB/virus_contigs/no_neg_extended_pruned_viral_renamed.fasta

cat ../../databases/UniVec_Core/UniVec_Core ../VIR_DB/virus_contigs/no_neg_extended_pruned_viral_renamed.fasta > \
	../VIR_DB/virus_contigs/w_UniVec_no_neg_extended_pruned_viral_renamed.fasta

# --- LOAD MODULES ---
module purge
module load BLAST+/2.13.0-gompi-2022a 
module list

# --- DEREPLICATION ACCORDING TO MIUViG GUIDELINES ---
mkdir -p ../VIR_DB/w_UniVec_dereplication

# First, create a blast+ database:
makeblastdb \
    -in ../VIR_DB/virus_contigs/w_UniVec_no_neg_extended_pruned_viral_renamed.fasta \
    -dbtype nucl \
    -out ../VIR_DB/w_UniVec_dereplication/w_UniVec_VIR_DB

# Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn \
    -query ../VIR_DB/virus_contigs/w_UniVec_no_neg_extended_pruned_viral_renamed.fasta \
    -db ../VIR_DB/w_UniVec_dereplication/w_UniVec_VIR_DB \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -out ../VIR_DB/w_UniVec_dereplication/w_UniVec_viruses_blast.tsv \
    -num_threads ${SLURM_CPUS_PER_TASK}

echo "all-vs-all blastn done!"

# --- LOAD MODULES --- 
module purge
module load Python/3.10.8-GCCcore-12.2.0
module load CheckV/1.0.1-foss-2021b-DIAMOND-2.1.8
module list

# Next, calculate pairwise ANI by combining local alignments between sequence pairs:
python /scratch/p282752/tools/checkv_scripts/anicalc.py \
	-i ../VIR_DB/w_UniVec_dereplication/w_UniVec_viruses_blast.tsv \
	-o ../VIR_DB/w_UniVec_dereplication/w_UniVec_viruses_ani.tsv

# anicalc.py is available at https://bitbucket.org/berkeleylab/checkv/src/master/scripts/anicalc.py

# Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
python /scratch/p282752/tools/checkv_scripts/aniclust.py \
    --fna ../VIR_DB/virus_contigs/w_UniVec_no_neg_extended_pruned_viral_renamed.fasta \
    --ani ../VIR_DB/w_UniVec_dereplication/w_UniVec_viruses_ani.tsv \
    --out ../VIR_DB/w_UniVec_dereplication/w_UniVec_viral_clusters.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0

# aniclust.py is available at https://bitbucket.org/berkeleylab/checkv/src/master/scripts/aniclust.py

Rscript under_dev_dereplication_stat.R /scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/w_UniVec_dereplication/w_UniVec_viral_clusters.tsv

module purge
