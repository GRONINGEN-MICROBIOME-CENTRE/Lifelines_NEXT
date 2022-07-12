### Calculate summary statistics of metadata file ###

# Creator: Paula Sureda | Arnau Vich
# Modified: Trishla Sinha, to suit needs for the Lifelnes NEXT cohort 
# Last update: 12 July 2022

## STEPS##

# Copy or import function to R #

# Input files: 
# File containing categorical or numerical phenotypes

#Output files
#Global summary table (for whole table input)

#Example metadata_input

#SID        factor1   factor2    factor3  factor4
#Sample1     23.5  	   0.9   	   yes        A
#Sample2     10.9 	   0.01  	   no         B
#Sample3     50    	   0.3    	 no         O

#Example output

#factors      Type      Categories/Median   Counts/Mean     %/SD    Number_non-zeros(n)  Number_NA
#factor1    numerical         4.5               4.32        1.22            34              6
#factor2    numerical         10                 9.34       0.38            38              2
#factor4   categorical     A, AB, B, O     244, 13, 52, 269    42.21, 2.25, 9, 46.54       578            136


summary_statistics_metadata <- function (metadata_input, category_table) {
  # Packages needed
  library (psych)  #describe r function
  # Create other functions to calculate the different parameters
  ## Categorical values - create function to calculate the counts and the percentage for categorical variables
  tblFun <- function(x) {
    # Create a table
    tbl <- table(x)
    # Combine columnes/rows to get the counts and percentage (creates new table -> res)
    res <- cbind(tbl,round(prop.table(tbl)*100,2))
    # Give names to the columns
    colnames(res) <- c('Count','Percentage')
    res
  }
  ## NA sum function - counts the number of NA
  nzsum <- function(x) {
    sum (is.na(x))
  }
  if (missing(category_table)) {
    ## Calculate table1 with the whole data:
    my_results = matrix(ncol = 9, nrow = ncol(metadata_input))
    for (k in 1:ncol(metadata_input)){
      if (is.numeric(metadata_input[,k])) {
        # Keep in "x" the result from describe function (done in the columns) - for each factor
        x = describe(metadata_input[,k])
        z = nzsum(metadata_input[,k])
        # In the new table ("x"): keep different values in the different columns
        my_results[k,1] = "numerical"
        my_results[k,2] = x$median
        my_results[k,3] = x$mean
        my_results[k,4] = x$sd
        my_results[k,5] = x$n
        my_results[k,6] = z
        my_results[k,7] = x$min
        my_results[k,8] = x$max
        my_results[k,9] = x$range
      }
      # Condition: if the column values are categorical
      else {
        # Keep in "x" the result from tblFun function (done in the columns) - for each factor
        x = tblFun(metadata_input[,k])
        z = nzsum(metadata_input[,k])
        # In the new table ("x"): keep different values in the different columns
        my_results[k,1]="categorical"
        # toString to keep the possible different values/categories in the same vector/column
        my_results[k,2]=toString(rownames(x))
        # First column table x = 'Count'
        my_results[k,3]=toString(x[,1])
        # Second column table x = 'Percentage'
        my_results[k,4]=toString(x[,2])
        # Sum of the values on column1 ("x")
        my_results[k,5]=sum(x[,1])
        my_results[k,6]= z
        my_results[k,7] = NA
        my_results[k,8] = NA
        my_results[k,9] = NA
      }
    }
    # The column names from the original table = row names from the new table
    rownames(my_results) = colnames(metadata_input)
    # Give names to the columns of the new table
    colnames(my_results) = c("Type", "Categories/Median", "Counts/Mean", "%/SD", "Number_non_zeros", "Number_NA",
                             "Min", "Max", "Range")
    # Export the new table
    write.table (my_results, file = "./meta_data_summary_stats.txt" , quote = F, sep = "\t")
  }
}


# USAGE 
summary_statistics_metadata (metadata_input)

