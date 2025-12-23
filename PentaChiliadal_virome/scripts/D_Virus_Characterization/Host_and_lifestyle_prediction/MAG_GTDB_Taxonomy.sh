#!/bin/bash
#SBATCH --job-name=iPHOP_DB_enrichment
#SBATCH --output=iPHOP_DB_enrichment.out
#SBATCH --mem=400gb
#SBATCH --time=10-0
#SBATCH --cpus-per-task=64
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

MAGs=$1 #directory with MAGs to be added to DB #(important to add TEST MAGs to our own MAGs)

# Activate conda environment
source /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/miniconda3/bin/activate \
    /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/GTDB_tk/gtdbtk-2.1.1; conda list

GTDBTK_DATA_PATH=/clusterfs/jgi/groups/science/metagen/shared_databases/gtdbtk_db/r214/release214/

# Generate GTDB-tk results
gtdbtk de_novo_wf \
    --genome_dir $MAGs \
    --bacteria \
    --outgroup_taxon p__Patescibacteria \
    --out_dir LLNEXT_GTDB-tk_results/ \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --force \
    --extension fa

gtdbtk de_novo_wf \
    --genome_dir $MAGs \
    --archaea \
    --outgroup_taxon p__Altiarchaeota \
    --out_dir LLNEXT_GTDB-tk_results/ \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --force \
    --extension fa

