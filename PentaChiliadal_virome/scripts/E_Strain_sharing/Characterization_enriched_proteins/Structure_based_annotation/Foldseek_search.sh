#!/bin/bash
#SBATCH --job-name=Foldseek_search
#SBATCH --output=Foldseek_search.out
#SBATCH --mem=60gb
#SBATCH --time=00-2:00 
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

PDB_dir=$1 # directory with PDB predictions from ColabFold
output_dir=$2 # output directory

echo -e '\n-------------------- LOOKING FOR HOMOLOGS with FOLDSEEK --------------------'

# Activate conda environment
source $JGI_SCRATCH/afernandezpato/Tools/miniconda3/bin/activate \
    $JGI_SCRATCH/afernandezpato/Tools/Foldseek; conda list

mkdir -p ${output_dir}

foldseek easy-search \
    $PDB_dir/ \
    /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/BFVD/bfvd \
    ${output_dir}/BFVD_results.tsv \
    ${output_dir}/BFVD_tmp \
    --alignment-type 1 \
    --format-output "query,target,pident,alnlen,evalue,bits,qlen,tlen,prob,lddt,qtmscore,ttmscore,alntmscore,qcov,tcov,qstart,qend,tstart,tend,qseq,tseq" \
    --threads ${SLURM_CPUS_PER_TASK}
 
foldseek easy-search \
    $PDB_dir/ \
    /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/PDB/pdb \
    ${output_dir}/PDB_results.tsv \
    ${output_dir}/PDB_tmp \
    --alignment-type 1 \
    --format-output "query,target,pident,alnlen,evalue,bits,qlen,tlen,prob,lddt,qtmscore,ttmscore,alntmscore,qcov,tcov,qstart,qend,tstart,tend,qseq,tseq" \
    --threads ${SLURM_CPUS_PER_TASK}
    
    foldseek easy-search \
    $PDB_dir/ \
    /clusterfs/jgi/scratch/science/metagen/afernandezpato/Tools/AFDB/afdb \
    ${output_dir}/AFDB_results.tsv \
    ${output_dir}/AFDB_tmp \
    --alignment-type 1 \
    --format-output "query,target,pident,alnlen,evalue,bits,qlen,tlen,prob,lddt,qtmscore,ttmscore,alntmscore,qcov,tcov,qstart,qend,tstart,tend,qseq,tseq" \
    --threads ${SLURM_CPUS_PER_TASK}
    
    # Remove intermediate files
    rm -r ${output_dir}/*tmp
 
