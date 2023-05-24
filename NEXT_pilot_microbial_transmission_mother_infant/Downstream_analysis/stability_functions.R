# Function to calculate stability metrics from initial timepoints
# Arguments:
# - metadata: Data frame containing metadata information (rows: samples, columns: phenotypes)
# - abundance_table: Abundance table of microbial species (rows: taxa, columns: samples)
# - timepoints: Vector of timepoints to consider 
# - sampleID: Column name in metadata containing sample IDs
# - calc_mode: Calculation mode ('abundance' or 'richness')
stability_initial <- function(metadata, abundance_table, timepoints, sampleID, calc_mode) {
  
  # Create an empty data frame to store the stability metrics
  stability_from_initial <- as.data.frame(matrix(NA, nrow=length( unique(metadata$Individual_ID) ), ncol=length(timepoints)*2 ))
  row.names(stability_from_initial) <- unique(metadata$Individual_ID)
  
  # Set column names based on the calculation mode
  if (calc_mode=='abundance') {
    colnames(stability_from_initial) <- c( timepoints, paste0('Total_space_', timepoints) )
    } else {
      colnames(stability_from_initial) <- c( timepoints, paste0('Richness_', timepoints) )
    }
  
  # Iterate over individuals
  for (i in row.names(stability_from_initial)) {
    
    # Iterate over timepoints
    for (timepoint in timepoints) {
      
      # Check if the individual and timepoint have samples
      if ( length(metadata[metadata$Individual_ID == i & metadata$Timepoint==timepoint,sampleID]) != 0 ) {
        
        # Calculate stability metric based on the calculation mode
        if (calc_mode=='abundance') {
          stability_from_initial[i,paste0('Total_space_',timepoint)] <- 100
          } else {
            stability_from_initial[i,paste0('Richness_',timepoint)] <- sum(abundance_table[, metadata[metadata$Individual_ID == i & metadata$Timepoint==timepoint, sampleID] ] != 0)
           }
      }
      
      # Get the sample IDs for the initial and current timepoints
      A <- metadata[metadata$Individual_ID == i & metadata$Timepoint == timepoints[1], sampleID]
      B <- metadata[metadata$Individual_ID == i & metadata$Timepoint == timepoint, sampleID]
      
      # Check if both initial and current timepoint have samples
      if (length(A) != 0 && length(B) != 0) {
        
        # Get the columns to consider in the abundance table
        columns <- c(A, B)
        
        # Calculate presence and copresence metrics
        presence <- rowSums(abundance_table[, columns] != 0) == 2
        copresent <- abundance_table[rowSums(abundance_table[, columns] != 0) == 2,columns]
        
        # Calculate stability metric based on the calculation mode
        if (calc_mode=='abundance') {
          stability_from_initial[i, timepoint] <- round(sum(copresent[,B])/sum(abundance_table[,B]) * 100, 2)
          } else {
          stability_from_initial[i, timepoint] <- sum(presence)
          }
        
      }
      
    }
    
  }
  
  # Calculate additional metrics for visualization
  for (i in timepoints) {
    
    if (calc_mode=='abundance') {
      stability_from_initial[,paste0('Not_retained_', i)] <- stability_from_initial[,paste0('Total_space_', i)] - stability_from_initial[,i]
    } else {
      stability_from_initial[,paste0('Not_retained_', i)] <- stability_from_initial[,paste0('Richness_', i)] - stability_from_initial[,i]
    }
 }
  
  # Return the stability metrics data frame
  return(stability_from_initial)
  
}


# Function to derive mean and sd values per column of tabels created with
# stability_initial() function
# This function also will calculate mean and 95% confidence interval quantiles
# using bootstrap
# Arguments:
# - df: Data frame containing metrics to summarize (rows: samples/individuals, columns: timepoints/metrics)
# - timepoints: Vector of timepoints to consider 
summary_stat_bootstrap <- function(df, timepoints) {
  
  # creating the data frame with means per column and sd per column
  basic_stat <- data.frame( colMeans(df, na.rm = T), apply(df, 2, sd, na.rm=T) )
  
  # assigning colnames to the data frame
  colnames(basic_stat) <- c("mean_value", "sd_value")
  
  # preparing the data frame for further visualization 
  basic_stat$Timepoint <- gsub(".*_", "", row.names(basic_stat))
  basic_stat$Timepoint <- factor(basic_stat$Timepoint, levels=timepoints, ordered = T)
  basic_stat$Condition <- gsub("_M.*|_P.*|_B$", "", row.names(basic_stat))
  basic_stat[basic_stat$Condition %in% timepoints,]$Condition <- 'Retained'
  
  # Set the number of bootstrap samples
  n_bootstrap <- 1000
  
  # calculate 95% confidence interval using bootstrap
  for (i in colnames(df)) {
    
    # Perform bootstrapping
    boot_samples <- replicate(n_bootstrap, sample(df[,i], replace = TRUE))
    
    # Calculate the mean for each bootstrap sample
    bootstrap_means <- apply(boot_samples, 2, mean, na.rm=T)
    
    # Calculate the 95% confidence interval and derive quantiles
    basic_stat[i,"q1_bootstrap"] <- unname(quantile(bootstrap_means, c(0.025), na.rm = T))
    basic_stat[i,"q2_bootstrap"] <- unname(quantile(bootstrap_means, c(0.975), na.rm = T))
    
    # Calculate bootstrapped mean
    basic_stat[i,'mean_bootstrap'] <- mean(bootstrap_means, na.rm=T)
  }
  
  # return the data frame
  return(basic_stat)
}

# Function to calculate prevalence of taxa that are persistent from the initial timepoint
# NOTE: function only considers samples from individuals that were taxa+ 
# at the previous timepoints
# Arguments:
# - metadata: Data frame containing metadata information (rows: samples, columns: phenotypes)
# - abundance_table: Abundance table of microbial species (rows: taxa, columns: samples)
# - timepoints: Vector of timepoints to consider 
# - sampleID: Column name in metadata containing sample IDs

persistent_taxa <- function(metadata, abundance_table, timepoints, sampleID) {
  
  # Create an empty data frame to store the taxa prevalence over time
  presence_from_initial <- as.data.frame(matrix(NA, nrow=nrow(abundance_table), ncol=length(timepoints) ) )
  
  # setting row names based on taxa names
  row.names(presence_from_initial) <- row.names(abundance_table)
  
  # setting colnames based on timepoints
  colnames(presence_from_initial) <- timepoints
  
  # creating empty list to store IDs of individuals that are taxa+ at the given timepoint
  positive_at_initial <- list()
  
  # iterate over taxa
  for (i in row.names(presence_from_initial)) {
    
    # assign the first (initial) timepoint
    j <- colnames(presence_from_initial)[1]
    
    # assign the list with sample IDs at the initial timepoint
    positive_at_initial[[i]] <- metadata[metadata$Timepoint==timepoints[1], sampleID]
    
    # loop to iterate over timepoints before the break with artificial timepoint
    while ( j != "M14" ) {
      
      # isolate the table of abundance for the given taxa and sample
      B <- abundance_table[i, metadata[ metadata$Timepoint == j, sampleID] ]
      
      # intersect IDs of individuals whose sample is taxa+ at the given timepoint with 
      # IDs of individuals that were taxa+ at the previous timepoint and keep them
      IDs_positive <- intersect(metadata[metadata[,sampleID] %in% colnames( B[i, B[i,]!=0 ]  ), ]$Individual_ID , 
                                metadata[metadata[,sampleID] %in% positive_at_initial[[i]], ]$Individual_ID  )
      
      # get the N of taxa+ individuals at the given timepoint
      presence_from_initial[i,j] <- length(IDs_positive)
      
      # re-assign the list with sample IDs of individuals that are taxa+ at the given timepoint
      positive_at_initial[[i]] <- metadata[ (metadata$Individual_ID %in% IDs_positive) & metadata$Timepoint==j, sampleID]
      
      # increment the timepoint for the loop to go to the next timepoint unless the last timepoint has been evaluated
      j <- ifelse( j !=timepoints[length(timepoints)] , colnames( presence_from_initial[ which(colnames(presence_from_initial) == j) + 1 ] ), j <- 'M14')
      
    }
    
  }
  
  # return the table of persistence
  return(presence_from_initial)
  
  }


# Function to fraction the personal microbiomes based in three categories
# PPB: Personal Persistent Bug, that is detected in >= 75% of individual's samples
# TDB: Transiently Detected Bug, that is deteceted in >= 2 of individual's samples
# Singletons: detected only in 1 of individual's samples
# The function will calculate the N of species in each category per individual and timepoint
# along with the abundance of species per category per timepoint
# Arguments:
# - metadata: Data frame containing metadata information (rows: samples, columns: phenotypes)
# - abundance_table: Abundance table of microbial species (rows: taxa, columns: samples)
# - sampleID: Column name in metadata containing sample IDs
personal_biome_fraction <- function(metadata, abundance_table, sampleID) {
  
  p_biome_frac <- as.data.frame(matrix(NA, nrow = length(unique(metadata$Individual_ID)), ncol=5))
  row.names(p_biome_frac) <- unique(metadata$Individual_ID)
  colnames(p_biome_frac) <- c('PPB', 'TDB', 'Singletons', 'IPB','N_timepoints')
  
  p_biome_frac$N_timepoints <- metadata$N_timepoints[match(row.names(p_biome_frac), metadata$Individual_ID)]
  
  PPB_bacteria <- list()
  TDB_bacteria <- list()
  Singletons_bacteria <- list()
  
  for (i in row.names(p_biome_frac) ) {
    
    A <- abundance_table[, metadata[metadata$Individual_ID==i,sampleID]  ]
    
    # N of bacterial species that are present in more than 75% of samples of the infant
    p_biome_frac[i,'PPB'] <- sum(rowSums(A!=0)>=ncol(A)*0.75)
    
    # personal persistent bacteria per infant:
    PPB_bacteria[[i]] <- row.names(A[rowSums(A!=0)>=ncol(A)*0.75,])
    
    # N of bacterial species that are present in less than 75% of samples of the infant but present in at least 2 samples
    p_biome_frac[i,'TDB'] <- sum(rowSums(A!=0) < ncol(A)*0.75 & rowSums(A!=0) >= 2)
    
    # personal transient bacteria per infant:
    TDB_bacteria[[i]] <- row.names(A[rowSums(A!=0) < ncol(A)*0.75 & rowSums(A!=0) >= 2,])
    
    # N of bacterial species that are present only in one sample of the individual
    p_biome_frac[i,'Singletons'] <- sum(rowSums(A!=0) == 1)
    
    # personal singletons:
    Singletons_bacteria[[i]] <- row.names(A[rowSums(A!=0) == 1,])
    
    # N of all bacterial species detected in the infant
    p_biome_frac[i, 'IPB'] <- sum(rowSums(A!=0) >= 1)
    
  }
  
  
  p_biome_frac_ab <- as.data.frame(matrix(NA, nrow = length(metadata[,sampleID]), ncol=6))
  row.names(p_biome_frac_ab) <- metadata[,sampleID]
  colnames(p_biome_frac_ab) <- c('ab_PPB', 'N_PPB','ab_TDB', 'N_TDB', 'ab_Singletons', 'N_Singletons')
  
  for (i in row.names(p_biome_frac_ab)) {
    A <- abundance_table[,i, drop=FALSE]
    
    p_biome_frac_ab[i,'ab_PPB'] <- sum(A[PPB_bacteria[[metadata[metadata[,sampleID]==i,]$Individual_ID]],])
    p_biome_frac_ab [i,'N_PPB'] <- sum(A[PPB_bacteria[[metadata[metadata[,sampleID]==i,]$Individual_ID]],]>0)
    
    
    p_biome_frac_ab[i,'ab_TDB'] <- sum(A[TDB_bacteria[[metadata[metadata[,sampleID]==i,]$Individual_ID]],])
    p_biome_frac_ab[i,'N_TDB'] <- sum(A[TDB_bacteria[[metadata[metadata[,sampleID]==i,]$Individual_ID]],]>0)
    
    p_biome_frac_ab[i,'ab_Singletons'] <- sum(A[Singletons_bacteria[[metadata[metadata[,sampleID]==i,]$Individual_ID]],])
    p_biome_frac_ab[i,'N_Singletons'] <- sum(A[Singletons_bacteria[[metadata[metadata[,sampleID]==i,]$Individual_ID]],]>0)
  }
  
  result <- list()
  result[["personal_biome_fractions"]] <- p_biome_frac
  result[["PPB_bacteria"]] <- PPB_bacteria
  result[["TDB_bacteria"]] <- TDB_bacteria
  result[["Singletons_bacteria"]] <- Singletons_bacteria
  result[["personal_biome_per_sample"]] <- p_biome_frac_ab
  return(result)
  
}





