extract_list_entry <- function(phenotype, factor) {
  # loop oover list
  for (i in seq_along(model_results)) {
    # entry where model_comparison is stored
    df <- model_results[[i]][[8]]
    
    # make sure I'm checking right columns
    if (is.data.frame(df) && "phenotype" %in% colnames(df) && "factor" %in% colnames(df)) {
      # Find matching indices in the current dataframe
      matching_indices <- which(df$phenotype == phenotype & df$factor == factor)
      
      # return any matching indicies. 
      if (length(matching_indices) > 0) {
        return(model_results[[i]])
      }
    }
  }
  # Return NULL if no match is found
  return(NULL)
}




