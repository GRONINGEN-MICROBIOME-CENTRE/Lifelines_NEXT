##########################################
# Combining bedtools coverage output
# to the coverage and count tables
##########################################

##############################
# Loading libraries
##############################
library(tidyverse)

##############################
# Functions
##############################

##############################
# Input data
##############################

##############################
# ANALYSIS
##############################
# Selection of columns according to the bedtools coverage standard ouput:
# V1: Contig names; V4: number of mapped reads; V7: breadth of coverage  
columns_to_keep <- c("V1", "V4", "V7")

temp = list.files(path="/scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/mapping/VLP_to_w_neg_der95/coverage",pattern="*.coverage.txt",full.names=T)

#myfiles = lapply(temp, 
#		FUN = function(x){
#		read.table(x, sep='\t', header=F) %>%
#    		select(columns_to_keep)
#		})

myfiles = lapply(temp,
                 FUN = function(x){
                   print(gsub('\\..*', '', basename(x)))
                   df <- read.table(x, sep='\t', header=F) 
                   df <- df[,columns_to_keep]
                   return(df) 
                 })

# Extraction of sample names from the file names
file_names = sapply(temp, function(x) gsub('\\..*', '', basename(x)))

# Naming the data frames according to the file names
myfiles = setNames(myfiles, file_names)

# Select and rename columns related to counts
for_count <- lapply(names(myfiles), function(list_name) {
  
  # Extract the dataframe
  df <- myfiles[[list_name]]
  
  # Select V1 and V4 & rename V4 according to the file name
  df_selected <- df %>% 
    select(V1, V4) %>%
    rename(!!list_name := V4)
  
  return(df_selected)
})

# Merge all selected and renamed dataframes into a single count table
count_table <- for_count %>% reduce(full_join, by='V1')

# Select and rename columns related to breadth of coverage
for_coverage <- lapply(names(myfiles), function(list_name) {
  
  # Extract the dataframe
  df <- myfiles[[list_name]]
  
  # Select V1 and V7 & rename V7 according to the file name
  df_selected <- df %>% 
    select(V1, V7) %>%
    rename(!!list_name := V7) 
  
  return(df_selected)
})

# Merge all selected and renamed dataframes into a single breadth of coverage talbe
coverage_table <- for_coverage %>% reduce(full_join, by='V1')

write.table(coverage_table, '/scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/mapping/VLP_to_w_neg_der95/coverage_table.txt', sep='\t', row.names=F, quote=F)
write.table(count_table, '/scratch/p282752/ANALYSIS_CHILIADAL/VIR_DB/mapping/VLP_to_w_neg_der95/count_table.txt', sep='\t', row.names=F, quote=F)
