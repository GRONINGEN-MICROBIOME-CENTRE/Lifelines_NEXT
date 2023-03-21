# Strainphlan 4.0 analysis Lifelines NEXT 

Adapted from Biobakery (StrainPhlAn 4.0). 
https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4

Authors: Trishla Sinha
Description: The script shows how strain profiling was performed using Strainphlan 4 for all maternal and infant samples post QC, for all species.   
Languages: Bash and R.   

## Step 1: Reconstruct all species strains

https://github.com/GRONINGEN-MICROBIOME-CENTRE/gmc-mgs-pipeline/blob/main/GMH_pipe.py 

## Step 2: Profile the clades present in the samples (profileClades.sh)

```
#!/bin/bash

#SBATCH --mem=24gb
#SBATCH --time=0-07:59:59
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
#SBATCH --job-name=SP4pr
#SBATCH --error=__SP4_profile.err
#SBATCH --output=__SP4_profile.out

# NOTES:
# script profiles all clades in the dataset in given folder ($1)
# puts results in the current folder!
# Adding --mutation_rates will give a mutation rates table for each of the alignes markers and a summary table for the concatenated MSA
# Removing the --print_clades only will actually run it 

# PARAMS
N=1 # --marker_in_n_samples
S=10 # --sample_with_n_markers
DB=/data/umcg-tifn/rgacesa/conda_biobakery4/lib/python3.10/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl 

# purge modules
module purge
# load conda
ml Miniconda3/4.8.3
# load conda env
source activate /data/umcg-tifn/rgacesa/conda_biobakery4
# run clade profiling
strainphlan -s *.pkl --database /data/umcg-tifn/rgacesa/conda_biobakery4/lib/python3.10/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl --marker_in_n_samples ${N} --sample_with_n_markers ${S} --print_clades_only --phylophlan_mode accurate --output_dir . > strainphlan4_clades_${N}.txt

```
### Execution 

```
sbatch ./profileClades.sh 

```
Mon Mar 20 17:47:58 2023: Start StrainPhlAn 4.0.6 execution
Mon Mar 20 17:47:58 2023: Loading MetaPhlAn mpa_vJan21_CHOCOPhlAnSGB_202103 database...
Mon Mar 20 17:48:20 2023: Done.
Mon Mar 20 17:48:23 2023: Detecting clades...
Mon Mar 20 18:40:14 2023: Done.
Mon Mar 20 18:40:14 2023: Detected clades: 
Mon Mar 20 18:40:14 2023:       t__SGB10068: in 103 samples.
Mon Mar 20 18:40:14 2023:       t__SGB8007_group: in 86 samples.
Mon Mar 20 18:40:14 2023:       t__SGB6936: in 82 samples.
Mon Mar 20 18:40:14 2023:       t__SGB17248: in 77 samples.
Mon Mar 20 18:40:14 2023:       t__SGB6952: in 74 samples.
Mon Mar 20 18:40:14 2023:       t__SGB6939: in 73 samples.
Mon Mar 20 18:40:14 2023:       t__SGB9712_group: in 60 samples.
Mon Mar 20 18:40:14 2023:       t__SGB10120: in 58 samples.
Mon Mar 20 18:40:14 2023:       t__SGB10115: in 58 samples.
Mon Mar 20 18:40:14 2023:       t__SGB17247: in 56 samples.
Mon Mar 20 18:40:14 2023:       t__SGB7962: in 55 samples.


We next process this output file to select only the clade names

```
cat strainphlan4_clades_1.txt | grep t__ | cut -f 2 | cut -f 1 -d ':' > LLNEXT_sp_clades_names.txt

```

This will give us the names of each species found: 
t__SGB10068
t__SGB8007_group
t__SGB6936
t__SGB17248
t__SGB6952
t__SGB6939
t__SGB9712_group
t__SGB10120
t__SGB10115
t__SGB17247
t__S6_group
t__SGB2303


## Step 3: Build the multiple sequence alignment (doMarkerComparison.sh)


```
#!/bin/bash

#SBATCH --mem=50gb
#SBATCH --time=0-09:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

# NOTES:
# > $1 is clade name

# PARAMS
N=4 # --marker_in_n_samples
S=10 # --sample_with_n_markers
DB=/data/umcg-tifn/rgacesa/conda_biobakery4/lib/python3.10/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl 

# purge modules
module purge
# load conda
ml Miniconda3/4.8.3
# load conda env
source activate /data/umcg-tifn/rgacesa/conda_biobakery4

mkdir ${1}
strainphlan -s *.pkl  --database /data/umcg-tifn/rgacesa/conda_biobakery4/lib/python3.10/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl --marker_in_n_samples ${N} --sample_with_n_markers ${S} --phylophlan_mode accurate --output_dir ./${1} --clade ${1} --nprocs 8


```

### Execution

```
for i in $(cat LLNEXT_sp_clades_names.txt); do sbatch doMarkerComparisonLLNext.sh $i; done 
```
This will perform MSA and create .tre files and .aln files for each of the species you feed it in 

*** Important *** Metaphlan and Strainphlan are not designed to work with phages and other viruses so even if these fit the cut-off do not work with them! 




## Step 4: Make distance matrix from MSA file

Example: 
bash ./makeDistMat.sh ./s__Alistipes_shahii/s__Alistipes_shahii.StrainPhlAn3_concatenated.aln

```
srun --pty -c5 --mem=8g --time=0-14:00 bash
for i in $(find . -type f -name *.aln); do bash makeDistMat.sh $i; done 
```
```
#!/bin/bash
echo 'maker of distance matrix from multiple alignment'
echo ' feed it with .aln file (multiple alignment)'

module purge
ml Anaconda3/5.3.0
# load conda env
source activate /groups/umcg-dag3/tmp01/rgacesa_tools/conda/envs/dag3pipe_v3_conda


# Creates a distance matrix from a multiple sequence alignment using the EMBOSS package (https://www.bioinformatics.nl/cgi-bin/emboss/help/distmat) 
# distmat calculates the evolutionary distance between every pair of sequences in a multiple sequence alignment.
# Uses Kimura Two-Parameter distance (distances expressed in terms of the number of substitutions per 100 b.p or amino acids) 
distmat -sequence ${1} -nucmethod 2 -outfile ${1/.aln/.dmat}

```
## Step 5: Cleaning the distance matrix from MSA file 

First, load RPlus
```
ml RPlus 
```

#Example: Rscript parseDMat_LLNext.R s__Bifidobacterium_bifidum.dmat
```
for i in $(find . -type f -name *.dmat); do Rscript parseDMat_LLNext.R $i; done 
```
#parseDMat_LLNext.R consists of: 
```
library(optparse)
# CL PARSER
help_description <- ""
args <- OptionParser(usage = "%prog clade.distmat.txt metadata.txt ordination.png",
                      add_help_option = TRUE, prog='strainphlan_ordination.R',
                      description=help_description )
args_list <- parse_args(args, positional_arguments=TRUE)

# read in the file, skipping the first 8 rows and filling in empty columns, using the tab as sep, and stripping extra white space
inFile <- args_list$args[1]
data <- read.table( inFile, skip = 8, fill = TRUE, sep="\t", strip.white = T)

# remove the first column of the data as it is blank
data[1] <- NULL

# get the header as the last column of the data as a character vector
header <- lapply(data[,ncol(data)], as.character)

# remove the last column from the data as it has been stored as a header
data[ncol(data)] <- NULL

# remove the current last column from the data as it is blank
data[ncol(data)] <- NULL

# split header by space and digit to get the sample names
samples <- unlist(header)
# fix metaphlan just by taking only first 3 thingies after _ split
ss <- c()
sNR <- 0
for (s in samples) {
   sNR <- sNR + 1
   ssplit = strsplit(s,' ')
   ss <- c(ss,ssplit[[1]][1])
}

# add the sample names to the columns and rows of the data matrix
rownames(data) <- ss
colnames(data) <- ss

# make symmetric, add lower triangle to upper triangle
data[lower.tri(data)] <- t(data)[lower.tri(data)]

# save it
write.table(data,paste0(gsub('\\.dmat','',inFile),'_dmat_Rready.csv'),sep=',',row.names=T)

```
# Ready csv's for export 
```
for i in $(find . -type f -name *dmat_Rready.csv); do cp $i /groups/umcg-llnext/tmp01/umcg-tsinha/strainphlan_19_08_2022/ready_csv_export_23_08/; done

```
