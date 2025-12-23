#!/bin/bash
#SBATCH --job-name=iPHOP_DB_enrichment
#SBATCH --output=iPHOP_DB_enrichment.out
#SBATCH --mem=150gb
#SBATCH --time=5-0
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=truncate
#SBATCH --partition=himem

MAGs=$1 #directory with MAGs to be added to DB #(important to add TEST MAGs to our own MAGs)
GTDB_DIR=$2 #directory with GTDB taxonomy results from the MAGs
iPHOP_dir=$3 #directory with iPHOP DBs 

# Clean environment, load modules and activate conda environment
module purge; ml Anaconda3; module list
source activate /scratch/hb-llnext/conda_envs/iphop133_env; conda list

 # Generate new iPHoP DB (GTDB genomes + MAGs)
iphop add_to_db \
    --fna_dir $MAGs \
    --gtdb_dir $GTDB_DIR \
    --out_dir ${iPHOP_dir}/Aug_2023_pub_rw_LLNEXT \
    --db_dir ${iPHOP_dir}/Aug_2023_pub_rw/

conda deactivate
