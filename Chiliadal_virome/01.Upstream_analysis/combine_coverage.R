##########################################
# Combining bedtools coverage output
# to the coverage and count tables
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################
library(tidyverse)
library(readr)
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
# X1: Contig names; X4: number of mapped reads; X7: breadth of coverage  
columns_to_keep <- c("X1", "X4", "X7")

temp = list.files(path=args[1], pattern=args[2], full.names=T)

#myfiles = lapply(temp, 
#		FUN = function(x){
#		read.table(x, sep='\t', header=F) %>%
#    		select(columns_to_keep)
#		})

myfiles = lapply(temp,
                 FUN = function(x){
                   print(gsub('\\..*', '', basename(x)))
                   tbl <- read_delim(x, col_names = F)
		   df <- as.data.frame(tbl)
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
  
  # Select X1 and X4 & rename X4 according to the file name
  df_selected <- df %>% 
    select(X1, X4) %>%
    rename(!!list_name := X4)
  
  return(df_selected)
})

# Merge all selected and renamed dataframes into a single count table
count_table <- for_count %>% reduce(full_join, by='X1')

# Select and rename columns related to breadth of coverage
for_coverage <- lapply(names(myfiles), function(list_name) {
  
  # Extract the dataframe
  df <- myfiles[[list_name]]
  
  # Select X1 and X7 & rename X7 according to the file name
  df_selected <- df %>% 
    select(X1, X7) %>%
    rename(!!list_name := X7) 
  
  return(df_selected)
})

# Merge all selected and renamed dataframes into a single breadth of coverage talbe
coverage_table <- for_coverage %>% reduce(full_join, by='X1')

write.table(coverage_table, paste0(args[1], "coverage_table.txt"),  sep='\t', row.names=F, quote=F)
write.table(count_table, paste0(args[1], "count_table.txt"), sep='\t', row.names=F, quote=F)
