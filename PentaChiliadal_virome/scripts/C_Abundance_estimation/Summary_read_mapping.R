# Get names of viral genomes and samples from command line arguments
args <- commandArgs(trailingOnly = TRUE)
viral_names <- readLines(args[1]) # vector with the names of representative viral genomes
sample_names <- readLines(args[2]) # vector of sample names
BED_output <- args[3] # path to BED coverage output
reads_table <- read.delim(args[4], header=T) # file with number of reads per sample

#Set sample names as row names in read_table (set colnames[2] to "clean_reads" if not done before)
rownames(reads_table) <- reads_table[, 1]

# Create an empty matrix to store the results
results_table <- matrix(0, nrow = length(viral_names), ncol = length(sample_names), dimnames = list(viral_names, sample_names))

# Loop through each sample in sample_names
for (sample in sample_names) {
    # Set the location of the result files from Bedtools coverage output
    file_path <- paste0(BED_output, sample, '.coverage.txt')
    # Read the table 
    cov_table <- read.table(file_path, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)
    # Subset the table to only include rows corresponding to the genomes in viral_names
    cov_table <- cov_table[viral_names, ]
    # Find the indices of viral genomes where the coverage is greater than or equal to 75%
    idx <- which(cov_table[, 6] >= 0.75)
    # Calculate the abundance values for results_table (column 3 = number of reads mapped; column 5 = viral sequence lenght)
    results_table[idx, sample] <- cov_table[idx, 3] / cov_table[idx, 5] / reads_table[sample, 'clean_reads'] * 10^9
}

# Write abundance table
write.table(results_table, file="LLNEXT_Abundance_Table_RPKM.txt", sep = "\t")
