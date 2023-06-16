setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore co-transmission of viruses and their 
# bacterial hosts
#############################################################

##############################
# Functions
##############################
flatten_correlation_matrix <- function(cor_matrix, pvalue_matrix) {
  # Get the row names and column names
  row_names <- row.names(cor_matrix)
  col_names <- colnames(cor_matrix)
  
  # Initialize an empty data frame to store the flattened matrix
  flattened <- data.frame(row = character(),
                          col = character(),
                          correlation = numeric(),
                          pvalue = numeric(),
                          stringsAsFactors = FALSE)
  
  # Iterate over each cell in the correlation matrix
  for (i in seq_along(row_names)) {
    for (j in seq_along(col_names)) {
      # Get the row name, column name, correlation value, and p-value
      row_name <- row_names[i]
      col_name <- col_names[j]
      correlation <- cor_matrix[i, j]
      pvalue <- pvalue_matrix[i, j]
      
      # Check if the correlation value and p-value are not NA
      if (!is.na(correlation) && !is.na(pvalue)) {
        # Create a new row in the flattened data frame
        new_row <- data.frame(row = row_name, col = col_name,
                              correlation = correlation, pvalue = pvalue)
        
        # Append the new row to the flattened data frame
        flattened <- rbind(flattened, new_row)
      }
    }
  }
  
  # Return the flattened data frame
  return(flattened)
}
##############################
# Loading libraries
##############################
#library(ggplot2)
library(vegan)
library(corrplot)

#library(ape)
#library(treeio)
#library(reshape2)
#library(ggtree)
#library(scales)
#library(ggforestplot)
#library(ggforce)

##############################
# Input data
##############################

selected_viruses <- read.table("02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt", sep='\t', header=T)

# selecting viruses and their hosts for the experiment:
check_cotransmission <- selected_viruses
# filter 0: select only those viruses that were used for the strain sharing analysis: 
check_cotransmission <- check_cotransmission[!is.na(check_cotransmission$Distances_Related_lower) & check_cotransmission$Distances_Related_lower=="YES",]
# filter 1: do not consider those viruses for which only 2 transmission/strain sharing events were detected
check_cotransmission <- check_cotransmission[check_cotransmission$Transmitted_in_N_related > 2,]
# filter 2: select only those hosts that were used for the strain sharing analysis: 
check_cotransmission <- check_cotransmission[!is.na(check_cotransmission$Host_Distances_Related_lower) & check_cotransmission$Host_Distances_Related_lower=="YES", ]

# cutpoint: 
check_cotransmission$cutpoint_virus <- ifelse(check_cotransmission$Youden_index >= check_cotransmission$FDR_ipv_Youden, check_cotransmission$FDR_ipv_Youden, check_cotransmission$Youden_index)
check_cotransmission$cutpoint_bacterium <- ifelse(check_cotransmission$Host_Youden_index >= check_cotransmission$Host_FDR_ipv_Youden, check_cotransmission$Host_FDR_ipv_Youden, check_cotransmission$Host_Youden_index)

# alternative easy name for host:
check_cotransmission$Host_easy <- paste0(check_cotransmission$Host_species, '_', gsub('SGB', '', check_cotransmission$Host_SGB) )

check_cotransmission <- check_cotransmission[order(check_cotransmission$Host_easy),]
#metadata <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/metadata_combined_for_exp.txt', sep='\t', header=T)
#metadata$Alex_ID <- paste0(metadata$FAM_ID, '_', metadata$Type, '_', substr(metadata$Short_sample_ID, 1,1), '_', metadata$source, '_', metadata$Timepoint)
#row.names(metadata) <- metadata$Alex_ID
#metadata$label <- metadata$Alex_ID

#RPKM_combined_0.95 <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/RPKM_counts_combined_0.95_UPD_final_for_exp.txt', sep='\t', header=T)

folder = "02.CLEAN_DATA/dist.matrices/"  
file_list = list.files(path=folder, pattern="*.txt")  
virus=lapply(paste0(folder, file_list), function(x) read.table(x, sep='\t', header=T))
names(virus) <- gsub(".dist.txt", '', file_list)
virus <- virus[names(virus) %in% check_cotransmission$Virus]

#### normalizing distance matrices by the median, as thresholds are calculated using the normalized distances:
virus <- lapply(virus, function(x) {
  x <- x/median( x[upper.tri(x)] )
})

#### exchanging distances values with 0 and 1 depending on the threshold for strain sharing event:
for (virusName in names(virus)) {
  
  virusN <- virus[[virusName]]
  
  virusN[virusN <= check_cotransmission[check_cotransmission$Virus==virusName,]$cutpoint_virus] <- 0
  
  virusN[virusN > check_cotransmission[check_cotransmission$Virus==virusName,]$cutpoint_virus] <- 1
  
  virus[[virusName]] <- virusN
}

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$source <- "MGS"
MGS_metadata$Alex_ID <- paste0(MGS_metadata$FAM_ID, '_', MGS_metadata$Type, '_', substr(MGS_metadata$Short_sample_ID_bact, 1,1), '_', MGS_metadata$source, '_', MGS_metadata$Timepoint)

folder = "01.RAW_DATA/strainphlan_4_distance_matrix/"  
file_list = list.files(path=folder, pattern="*.csv")  
bacterium=lapply(paste0(folder, file_list), function(x) read.csv(x, header=T))
names(bacterium) <- gsub("_dmat_.csv", '', file_list)
bacterium <- bacterium[names(bacterium) %in% check_cotransmission$Host_SGB]
bacterium <- lapply(bacterium, function(x) {
  x <- x/median( x[upper.tri(x)], na.rm=T )
})

## since my parsing scripts have certain ID format: 
bacterium <- lapply(bacterium, function(x) {
  row.names(x) <- substr(row.names(x), 0, 12) # first remove data processing artifacts
  colnames(x) <- substr(colnames(x), 0, 12) 
  
  row.names(x) <- MGS_metadata$Alex_ID[match(row.names(x), MGS_metadata$NG_ID)]
  colnames(x) <- MGS_metadata$Alex_ID[match(colnames(x), MGS_metadata$NG_ID)] # then change the sequencing IDs to IDs for parsing
  x
})

#### exchanging distances values with 0 and 1 depending on the threshold for strain sharing event:
for (bacteriumName in names(bacterium)) {
  
  bacteriumN <- bacterium[[bacteriumName]]
  
  bacteriumN[bacteriumN <= check_cotransmission[check_cotransmission$Host_SGB==bacteriumName,]$cutpoint_bacterium] <- 0
  
  bacteriumN[bacteriumN > check_cotransmission[check_cotransmission$Host_SGB==bacteriumName,]$cutpoint_bacterium] <- 1
  
  bacterium[[bacteriumName]] <- bacteriumN
}

#metaphlan <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)
#metaphlan_species <- metaphlan[grep('s__', row.names(metaphlan)),]

#virus_host_transmission <- read.table("03a.RESULTS/Bacterial_hosts_transmission_for_transmitted_viruses.txt", sep='\t', header=T)
#virus_host_transmission_reconstructed <- virus_host_transmission[!is.na(virus_host_transmission$easy_name),]

#virus_metadata <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/Viruses_shared_min_5_families_UPD_final.txt', sep='\t', header=T)
##############################
# ANALYSIS
##############################

co_transmission_cor <- data.frame(matrix(NA, ncol=length(unique(check_cotransmission$Virus)),
                                     nrow=length(unique(check_cotransmission$Host_SGB)))  )

colnames(co_transmission_cor) <- unique(check_cotransmission$Virus)
row.names(co_transmission_cor) <- unique(check_cotransmission$Host_SGB)

co_transmission_pval <- co_transmission_cor

mantel_partial_results <- list()

# calculating the correlation between matrices of strain sharing reconstructed for viruses and their predicted hosts:

# Iterate over each virus in the dataset
for (virusName in check_cotransmission$Virus) {
  
  # Obtain the virus genetic distance matrix
  virusN <- virus[[virusName]]
  
  # Parse sample names to handle cases where virus strain was reconstructed in both MGS and VLP
  sample_names <- data.frame(colnames(virusN))
  colnames(sample_names) <- 'sample'
  sample_names$Var1 <- gsub("VLP_|MGS_", "", sample_names$sample)
  sample_names$source <- NA
  sample_names[grep('VLP', sample_names$sample),]$source <- 'VLP'
  sample_names[grep('MGS', sample_names$sample),]$source <- 'MGS'
  repeated <- data.frame(table(gsub("VLP_|MGS_", "", colnames(virusN))))
  sample_names <- merge(sample_names, repeated, by='Var1')
  sample_names <- sample_names[sample_names$Freq==1 | sample_names$source=='VLP',]
  sample_names$bacterial <- gsub('VLP', 'MGS', sample_names$sample)
  
  # Iterate over each bacterial host strain predicted to be infected by the virus:
  for (bacteriumName in check_cotransmission$Host_SGB) {
    
    # Obtain the bacterial distance matrix
    bacteriumN <- bacterium[[bacteriumName]]
    
    # Merge names of samples with reconstructed virus and bacterial strains
    sample_names_bacterial <- data.frame(colnames(bacteriumN))
    colnames(sample_names_bacterial) <- 'bacterial'
    sample_names_concurrent <- merge(sample_names, sample_names_bacterial, by='bacterial')
    
    # Extract overlapping samples where both virus and bacterial strains were reconstructed
    virusN_concurrent <- virusN[sample_names_concurrent$sample, sample_names_concurrent$sample]
    bacteriumN <- bacteriumN[sample_names_concurrent$bacterial, sample_names_concurrent$bacterial]
    
    if (length(sample_names_concurrent$bacterial)>=7) {
      
      dist_matrix_contol <- virusN_concurrent
      
      for (i in colnames(dist_matrix_contol)) {
        
        for (j in row.names(dist_matrix_contol)) {
          
          if ( substr(i, 1, 16)==substr(j, 1, 16) ) {
            dist_matrix_contol[i,j] <- 0
          } else {
            dist_matrix_contol[i,j] <- 1
          }
          
        }
        
      }
      
      
      mantel_partial_results[[paste0(virusName, '_', bacteriumName)]] <- mantel.partial(virusN_concurrent, bacteriumN, dist_matrix_contol, method = "pearson", permutations = 999, na.rm = T)
      
      # Store the correlation coefficient and p-value
      co_transmission_cor[bacteriumName,virusName] <- mantel_partial_results[[paste0(virusName, '_', bacteriumName)]]$statistic
      co_transmission_pval[bacteriumName,virusName] <- mantel_partial_results[[paste0(virusName, '_', bacteriumName)]]$signif
      
    } else {
      mantel_partial_results[[paste0(virusName, '_', bacteriumName)]] <- 'No concurrent samples'
    }

  }
  
}

##############################
# INTERMEDIATE OUTPUT
##############################
# since running the mantel.partial() takes substantial time, raw output should be saved:
write.table(co_transmission_cor, "03a.RESULTS/Co_transmission_logical_pearson_viruses_vs_bacteria_corr_raw.txt", sep='\t', quote=F)
write.table(co_transmission_pval, "03a.RESULTS/Co_transmission_logical_pearson_viruses_vs_bacteria_pval_raw.txt", sep='\t', quote=F)
##############################
# INTERMEDIATE OUTPUT
##############################

colnames(co_transmission_cor) <- check_cotransmission$ContigID_easy[match(colnames(co_transmission_cor), check_cotransmission$Virus)]
row.names(co_transmission_cor) <- check_cotransmission$Host_easy[match(row.names(co_transmission_cor), check_cotransmission$Host_SGB)]

colnames(co_transmission_pval) <- check_cotransmission$ContigID_easy[match(colnames(co_transmission_pval), check_cotransmission$Virus)]
row.names(co_transmission_pval) <- check_cotransmission$Host_easy[match(row.names(co_transmission_pval), check_cotransmission$Host_SGB)]

co_transmission_cor <- co_transmission_cor[,colSums(is.na(co_transmission_cor)) != nrow(co_transmission_cor)]
co_transmission_cor <- co_transmission_cor[rowSums(is.na(co_transmission_cor))!=ncol(co_transmission_cor),]

co_transmission_cor <- co_transmission_cor[order(row.names(co_transmission_cor)),]
co_transmission_pval <- co_transmission_pval[order(row.names(co_transmission_pval)),]

co_transmission_cor <- as.matrix(co_transmission_cor)
co_transmission_pval <- as.matrix(co_transmission_pval)

# multiple test correction:
cotransmission_flat <- flatten_correlation_matrix(co_transmission_cor, co_transmission_pval)
cotransmission_flat$FDR <- p.adjust(cotransmission_flat$pvalue, method = "BH")
co_transmission_FDR <- co_transmission_pval

for (i in 1:nrow(cotransmission_flat)) {
  
  virusName <- cotransmission_flat[ i,  'col']
  bacteriumName <- cotransmission_flat[ i,  'row']
  
  co_transmission_FDR[bacteriumName,virusName] <- cotransmission_flat[i, 'FDR']
  
}

co_transmission_cor[is.nan(co_transmission_cor) | !is.finite(co_transmission_cor)] <- NA 
co_transmission_FDR[is.nan(co_transmission_FDR) | !is.finite(co_transmission_FDR)] <- NA
co_transmission_FDR[is.na(co_transmission_FDR)] <- NA

bgcolors <- matrix(0, nrow = nrow(co_transmission_cor), 
                   ncol = ncol(co_transmission_cor), 
                   dimnames = list(row.names(co_transmission_cor), colnames(co_transmission_cor)))

for (virusName in check_cotransmission$ContigID_easy) {
  
  bacteriumName <- check_cotransmission[check_cotransmission$ContigID_easy==virusName,]$Host_easy
  
  bgcolors[bacteriumName, virusName] <- NA
}

# some mantel exact stat artifact I guess:
co_transmission_cor[co_transmission_cor>1] <- 1

pdf('./04.PLOTS/Virus_bacteria_co_transmission_corrplot_FDR_0.1.pdf', width=25/2.54, height=22/2.54)
co_transmission_cor[co_transmission_FDR > 0.10] <- 0
corrplot(bgcolors, na.label = "square", na.label.col = "#E8AA42", tl.col = "white")
corrplot(co_transmission_cor,
         p.mat=co_transmission_FDR, 
         insig = "blank",
         bg=NA,
         na.label = "X", na.label.col = "grey", tl.col='black', add = T)
dev.off()

# test if co-transmission is enriched in predicted virus-host pairs vs random virus-bacteria pairs:

cotransmission_flat$virus_host_pair <- "NO"
for (i in 1:nrow(cotransmission_flat)) {
  
  bacteriumName <- cotransmission_flat[i,"row"]
  virusName <- cotransmission_flat[i,"col"]
  
  if (bacteriumName %in% check_cotransmission[check_cotransmission$ContigID_easy==virusName,"Host_easy"] ) {
    cotransmission_flat[i, "virus_host_pair"] <- "YES"
  }
  
}

dim(cotransmission_flat[cotransmission_flat$FDR <= 0.05 & cotransmission_flat$virus_host_pair=="YES",])[1]
dim(cotransmission_flat[cotransmission_flat$FDR > 0.05 & cotransmission_flat$virus_host_pair=="YES",])[1]
dim(cotransmission_flat[cotransmission_flat$FDR <= 0.05 & cotransmission_flat$virus_host_pair=="NO",])[1]
dim(cotransmission_flat[cotransmission_flat$FDR > 0.05 & cotransmission_flat$virus_host_pair=="NO",])[1]

cotransmission_stat <- matrix( c(19, # pair & cotransmitted
                                 29, # pair & not cotransmitted
                                 81, # not a pair & cotransmitted
                                 269), # not a pair & not cotransmittd
                                 nrow=2,
                                 dimnames = list(c("Cotransmitted", "Not-cotransmitted"),
                                             c("Paired", "Unpaired")))

fisher.test(cotransmission_stat, "greater")


#### calculating the background: are bacterial strains co-transmitted with each other? 

 cotransmission_BB_cor <- matrix(NA, ncol=length(unique(check_cotransmission$Host_SGB)),
                                              nrow=length(unique(check_cotransmission$Host_SGB)),
                                 dimnames = list(unique(check_cotransmission$Host_SGB), unique(check_cotransmission$Host_SGB)))  

 cotransmission_BB_pval <-  cotransmission_BB_cor
 
 mantel_partial_BB_results <- list()
 
 for (bacteriumName in colnames(cotransmission_BB_cor)) {
   
   # Obtain the virus genetic distance matrix
   bacteriumN <- bacterium[[bacteriumName]]
   
   # Parse sample names to handle cases where virus strain was reconstructed in both MGS and VLP
   sample_names <- colnames(bacteriumN)
   
   # Iterate over each bacterial host strain predicted to be infected by the virus:
   for (bacteriumName2 in colnames(cotransmission_BB_cor)) {
     
     # Obtain the bacterial distance matrix
     bacteriumN2 <- bacterium[[bacteriumName2]]
     
     # Merge names of samples with reconstructed virus and bacterial strains
     sample_names_bacterial <- colnames(bacteriumN2)
     
     sample_names_concurrent <- intersect(sample_names, sample_names_bacterial)
     
     # Extract overlapping samples where both virus and bacterial strains were reconstructed
     bacteriumN_concurrent <- bacteriumN[sample_names_concurrent, sample_names_concurrent]
     bacteriumN2 <- bacteriumN2[sample_names_concurrent, sample_names_concurrent]
     
     if (length(sample_names_concurrent)>=7) {
       
       dist_matrix_contol <- bacteriumN_concurrent
       
       for (i in colnames(dist_matrix_contol)) {
         
         for (j in row.names(dist_matrix_contol)) {
           
           if ( substr(i, 1, 16)==substr(j, 1, 16) ) {
             dist_matrix_contol[i,j] <- 0
           } else {
             dist_matrix_contol[i,j] <- 1
           }
           
         }
         
       }
       
       mantel_partial_BB_results[[paste0(bacteriumName, '_', bacteriumName2)]] <- mantel.partial(bacteriumN_concurrent, bacteriumN2, dist_matrix_contol, method = "pearson", permutations = 999, na.rm = T)
       
       # Store the correlation coefficient and p-value
       cotransmission_BB_cor[bacteriumName2,bacteriumName] <- mantel_partial_BB_results[[paste0(bacteriumName, '_', bacteriumName2)]]$statistic
       cotransmission_BB_pval[bacteriumName2,bacteriumName] <- mantel_partial_BB_results[[paste0(bacteriumName, '_', bacteriumName2)]]$signif
       
     } else {
       mantel_partial_BB_results[[paste0(bacteriumName, '_', bacteriumName2)]] <- 'No concurrent samples'
     }
     
   }
   
 }
 
 ##############################
 # INTERMEDIATE OUTPUT
 ##############################
 # since running the mantel.partial() takes substantial time, raw output should be saved:
 write.table(cotransmission_BB_cor, "03a.RESULTS/Co_transmission_logical_pearson_bacteria_vs_bacteria_corr_raw.txt", sep='\t', quote=F)
 write.table(cotransmission_BB_pval, "03a.RESULTS/Co_transmission_logical_pearson_bacteria_vs_bacteria_pval_raw.txt", sep='\t', quote=F)
 ##############################
 # INTERMEDIATE OUTPUT
 ##############################
 
 colnames(cotransmission_BB_cor) <- check_cotransmission$Host_easy[match(colnames(cotransmission_BB_cor), check_cotransmission$Host_SGB)]
 row.names(cotransmission_BB_cor) <- check_cotransmission$Host_easy[match(row.names(cotransmission_BB_cor), check_cotransmission$Host_SGB)]
  
 colnames(cotransmission_BB_pval) <- check_cotransmission$Host_easy[match(colnames(cotransmission_BB_pval), check_cotransmission$Host_SGB)]
 row.names(cotransmission_BB_pval) <- check_cotransmission$Host_easy[match(row.names(cotransmission_BB_pval), check_cotransmission$Host_SGB)]
  
 cotransmission_BB_cor <- cotransmission_BB_cor[,colSums(is.na(cotransmission_BB_cor)) != nrow(cotransmission_BB_cor)]
 cotransmission_BB_cor <- cotransmission_BB_cor[rowSums(is.na(cotransmission_BB_cor))!=ncol(cotransmission_BB_cor),]

  # since the matrix of correlations is symmetric:   
cotransmission_BB_cor[upper.tri( cotransmission_BB_cor, diag = T)] <- NA
cotransmission_BB_pval[upper.tri( cotransmission_BB_pval,diag = T)] <- NA

# multiple test correction:
cotransmission_BB_flat <- flatten_correlation_matrix(cotransmission_BB_cor, cotransmission_BB_pval)
cotransmission_BB_flat$FDR <- p.adjust(cotransmission_BB_flat$pvalue, method = "BH")

cotransmission_BB_FDR <- cotransmission_BB_pval

 for (i in 1:nrow(cotransmission_BB_flat)) {

   bacteriumName1 <- cotransmission_BB_flat[ i,  'col']
   bacteriumName2 <- cotransmission_BB_flat[ i,  'row']

   cotransmission_BB_FDR[bacteriumName2,bacteriumName1] <- cotransmission_BB_flat[i, 'FDR']

 }
  
cotransmission_BB_cor[is.nan(cotransmission_BB_cor) | !is.finite(cotransmission_BB_cor)] <- NA 
cotransmission_BB_cor[is.nan(cotransmission_BB_cor) | !is.finite(cotransmission_BB_cor)] <- NA
cotransmission_BB_FDR[is.na(cotransmission_BB_cor)] <- NA

indices <- which(cotransmission_BB_cor < -1 | cotransmission_BB_cor > 1, arr.ind = TRUE)

# Print the indices and corresponding values
for (i in 1:nrow(indices)) {
  row <- indices[i, 1]
  col <- indices[i, 2]
  value <- cotransmission_BB_cor[row, col]
  cat("Element at [", row, ",", col, "] is", value, "\n")
}

  # some mantel exact stat artifact I guess:
cotransmission_BB_cor[cotransmission_BB_cor>1] <- 1
cotransmission_BB_cor[cotransmission_BB_cor< (-1)] <- -1
  
pdf('./04.PLOTS/Bacteria_bacteria_co_transmission_corrplot_FDR_0.1_raw.pdf', width=25/2.54, height=22/2.54)
cotransmission_BB_cor[cotransmission_BB_FDR > 0.05] <- 0

 corrplot(cotransmission_BB_cor,
          type = 'lower',
          bg=NA,
          na.label = "X", na.label.col = "grey", tl.col='black')
 
 
 
  dev.off()
 # 
 # # test if co-transmission is enriched in predicted virus-host pairs vs random virus-bacteria pairs:
 # 
 # cotransmission_flat$virus_host_pair <- "NO"
 # for (i in 1:nrow(cotransmission_flat)) {
 #   
 #   bacteriumName <- cotransmission_flat[i,"row"]
 #   virusName <- cotransmission_flat[i,"col"]
 #   
 #   if (bacteriumName %in% check_cotransmission[check_cotransmission$ContigID_easy==virusName,"Host_easy"] ) {
 #     cotransmission_flat[i, "virus_host_pair"] <- "YES"
 #   }
 #   
 # }
 # 
 # dim(cotransmission_flat[cotransmission_flat$FDR <= 0.05 & cotransmission_flat$virus_host_pair=="YES",])[1]
 # dim(cotransmission_flat[cotransmission_flat$FDR > 0.05 & cotransmission_flat$virus_host_pair=="YES",])[1]
 # dim(cotransmission_flat[cotransmission_flat$FDR <= 0.05 & cotransmission_flat$virus_host_pair=="NO",])[1]
 # dim(cotransmission_flat[cotransmission_flat$FDR > 0.05 & cotransmission_flat$virus_host_pair=="NO",])[1]
 # 
 # cotransmission_stat <- matrix( c(19, # pair & cotransmitted
 #                                  29, # pair & not cotransmitted
 #                                  81, # not a pair & cotransmitted
 #                                  269), # not a pair & not cotransmittd
 #                                nrow=2,
 #                                dimnames = list(c("Cotransmitted", "Not-cotransmitted"),
 #                                                c("Paired", "Unpaired")))
 # 
 # fisher.test(cotransmission_stat, "greater")
 # colnames(co_transmission_cor) <- unique(virus_host_transmission_reconstructed$Virus)
# row.names(co_transmission_cor) <- unique(virus_host_transmission_reconstructed$Bacterium)
# 
# co_transmission_pval <- co_transmission_cor
# 
# mantel_results <- list()
# for (n in virus_host_transmission_reconstructed$Virus) {
#   
#   # virus distance matrix:
#   virusN <- virus[[n]]
#   
#   # normalizing the distances:
#   virusN <- virusN/median(unname(unlist(virusN)))
#   
#   # sample names parsing (because in some cases virus strain was reconstructed both in MGS and VLP simultaneously):
#   sample_names <- data.frame(colnames(virusN))
#   colnames(sample_names) <- 'sample'
#   sample_names$Var1 <- gsub("VLP_|MGS_", "", sample_names$sample)
#   sample_names$source <- NA
#   sample_names[grep('VLP', sample_names$sample),]$source <- 'VLP'
#   sample_names[grep('MGS', sample_names$sample),]$source <- 'MGS'
#   repeated <- data.frame(table(gsub("VLP_|MGS_", "", colnames(virusN))))
#   sample_names <- merge(sample_names, repeated, by='Var1')
#   sample_names <- sample_names[sample_names$Freq==1 | sample_names$source=='VLP',]
#   sample_names$bacterial <- gsub('VLP', 'MGS', sample_names$sample)
#   
#   
#   
#   # in some cases, multiple hosts are predicted for 1 virus: 
#   for (b in virus_host_transmission_reconstructed$Bacterium) {
#     
#     bacteriumN <- bacterium_to_test[[b]]
#     # normalizing the distance:
#     bacteriumN <- bacteriumN/median(unname(unlist(bacteriumN)), na.rm=T)
#     
#     sample_names_bacterial <- data.frame(colnames(bacteriumN))
#     colnames(sample_names_bacterial) <- 'bacterial'
#     
#     sample_names_concurrent <- merge(sample_names, sample_names_bacterial, by='bacterial')
#     
#     
#     # overalpping samples where both virus and bacterial strains were reconstructed:
#     virusN_concurrent <- virusN[sample_names_concurrent$sample, sample_names_concurrent$sample, drop=F]
#     bacteriumN <- bacteriumN[sample_names_concurrent$bacterial, sample_names_concurrent$bacterial, drop=F]
#     
#     if (length(sample_names_concurrent$bacterial)>1 & length(sample_names_concurrent$sample)>1) {
#       mantel_results[[paste0(n, '_', b)]] <- mantel(virusN_concurrent, bacteriumN, method = "spearman", permutations = 999, na.rm = T)
#       
#       co_transmission_cor[b,n] <- mantel_results[[paste0(n, '_', b)]]$statistic
#       co_transmission_pval[b,n] <- mantel_results[[paste0(n, '_', b)]]$signif
#       
#     } else {
#       mantel_results[[paste0(n, '_', b)]] <- 'No concurrent samples'
#     }
#     
#     
#     
#   }
#   
# }



# colnames(co_transmission_cor) <- gsub('.*length', 'L', colnames(co_transmission_cor))
# keep <- co_transmission_cor
# 
# co_transmission_cor <- keep
# #co_transmission_cor[co_transmission_pval>0.05] <- NA
# co_transmission_cor <- co_transmission_cor[,colSums(is.na(co_transmission_cor)) != nrow(co_transmission_pval)]
# co_transmission_cor <- co_transmission_cor[rowSums(is.na(co_transmission_cor))!=ncol(co_transmission_pval),]
# 
# virus_metadata$easy_name <- gsub('.*length', 'L', virus_metadata$V1)
# colnames(co_transmission_cor) <- virus_metadata$ContigID_easy[match(colnames(co_transmission_cor), virus_metadata$easy_name)]
# row.names(co_transmission_cor) <- virus_host_transmission_reconstructed$easy_name[match(row.names(co_transmission_cor), virus_host_transmission_reconstructed$Bacterium)]


# bgcolors <- matrix("white", nrow = length(unique(virus_host_transmission_reconstructed$Bacterium)), 
#                             ncol = length(unique(virus_host_transmission_reconstructed$Virus)))
# colnames(bgcolors) <- unique(virus_host_transmission_reconstructed$Virus)
# row.names(bgcolors) <- unique(virus_host_transmission_reconstructed$Bacterium)
# 
# for (virus_name in virus_host_transmission_reconstructed$Virus) {
#   
#   bacterium_name <- virus_host_transmission_reconstructed[virus_host_transmission_reconstructed$Virus==virus_name,]$Bacterium
#   
#   bgcolors[bacterium_name, virus_name] <- NA
# }


# colnames(co_transmission_pval) <- virus_metadata$ContigID_easy[match(colnames(co_transmission_pval), virus_metadata$V1)]
# row.names(co_transmission_pval) <- virus_host_transmission_reconstructed$easy_name[match(row.names(co_transmission_pval), virus_host_transmission_reconstructed$Bacterium)]
# 
# colnames(bgcolors) <- virus_metadata$ContigID_easy[match(colnames(bgcolors), virus_metadata$V1)]
# row.names(bgcolors) <- virus_host_transmission_reconstructed$easy_name[match(row.names(bgcolors), virus_host_transmission_reconstructed$Bacterium)]
# bgcolors <- bgcolors[row.names(co_transmission_cor),colnames(co_transmission_cor)]


co_transmission_dummy <- matrix(0, ncol=ncol(co_transmission_cor), nrow=nrow(co_transmission_cor), dimnames = dimnames(co_transmission_cor))
co_transmission_dummy[is.na(bgcolors)] <- NA

keep2 <- co_transmission_pval


# tmp <- flatten_correlation_matrix(co_transmission_cor, co_transmission_pval)
# tmp$FDR <- p.adjust(tmp$pvalue, method = "BH")
# tmp$unique_pair <- paste0(tmp$row, '_', tmp$col)
# 
# co_transmission_FDR <- co_transmission_pval
# 
# for (i in tmp$unique_pair) {
#   
#   virusN <- tmp[ tmp$unique_pair==i,  'col']
#   bacteriumN <- tmp[ tmp$unique_pair==i,  'row']
#   
#   co_transmission_FDR[bacteriumN,virusN] <- tmp[tmp$unique_pair==i, 'FDR']
#   
# }




significance_level <- 0.05
#create the color matrix from the p-value matrix, select only necessary data
pdf('./04.PLOTS/Virus_bacteria_co_transmission_corrplot_no_correction.pdf', width=25/2.54, height=22/2.54)
mycol <-ifelse(c(co_transmission_FDR > significance_level), "white", "black")
co_transmission_cor[co_transmission_pval > 0.05] <- 0
corrplot(bgcolors, na.label = "square", na.label.col = "#E8AA42", tl.col = "white")
corrplot(co_transmission_cor,
         p.mat=co_transmission_FDR, 
         insig = "blank", sig.level = 0.05,
         bg=NA,
         na.label = "X", na.label.col = "grey", tl.col='black', add = T)
dev.off()




corrplot(as.matrix(co_transmission_cor), 
         p.mat=as.matrix(co_transmission_FDR), 
         insig = "blank",sig.level = 0.05, na.label = "X", na.label.col = "grey", bg=NA, add = T)




#### calculating the co-transmission of bacteria:
co_transmission_bacteria <- data.frame(matrix(NA, ncol=length(unique(virus_host_transmission_reconstructed$Bacterium)),
                                         nrow=length(unique(virus_host_transmission_reconstructed$Bacterium)))  )
colnames(co_transmission_bacteria) <- unique(virus_host_transmission_reconstructed$Bacterium)
row.names(co_transmission_bacteria) <- unique(virus_host_transmission_reconstructed$Bacterium)

co_transmission_bacpval <- co_transmission_bacteria

mantel_results_bacteria <- list()
for (n in virus_host_transmission_reconstructed$Bacterium) {
  
  # virus distance matrix:
  bacteriumN <- bacterium_to_test[[n]]
  
  # normalizing the distances:
  bacteriumN <- bacteriumN/median(unname(unlist(bacteriumN)), na.rm=T)
  
  # in some cases, multiple hosts are predicted for 1 virus: 
  for (b in virus_host_transmission_reconstructed$Bacterium) {
    
    bacteriumK <- bacterium_to_test[[b]]
    # normalizing the distance:
    bacteriumK <- bacteriumK/median(unname(unlist(bacteriumK)), na.rm=T)
    
    concurrent_samples <- intersect(colnames(bacteriumN), colnames(bacteriumK))
    
    bacteriumN_concurrent <- bacteriumN[concurrent_samples, concurrent_samples]
    bacteriumK <- bacteriumK[concurrent_samples, concurrent_samples]
    
    if (length(concurrent_samples)>1) {
      
      mantel_results_bacteria[[paste0(n, '_', b)]] <- mantel(bacteriumN_concurrent, bacteriumK, method = "spearman", permutations = 999, na.rm = T)
      
      co_transmission_bacteria[b,n] <- mantel_results_bacteria[[paste0(n, '_', b)]]$statistic
      co_transmission_bacpval[b,n] <- mantel_results_bacteria[[paste0(n, '_', b)]]$signif
      
    } else {
      mantel_results_bacteria[[paste0(n, '_', b)]] <- 'No concurrent samples'
    }
    
    
    
  }
  
}

keep_co_transmission_bacteria <- co_transmission_bacteria
keep_co_transmission_bacpval <- co_transmission_bacpval

pdf('./04.PLOTS/Bacteria_bacteria_co_transmission_corrplot_no_correction.pdf', width=30/2.54, height=30/2.54)
corrplot(as.matrix(co_transmission_bacteria), 
         p.mat=as.matrix(co_transmission_bacpval), 
         insig = "blank",sig.level = 0.05, na.label = "X", na.label.col = "grey")
dev.off()
######## 

colnames(co_transmission_bacteria) <- virus_host_transmission_reconstructed$easy_name[match(colnames(co_transmission_bacteria), virus_host_transmission_reconstructed$Bacterium)]
row.names(co_transmission_bacteria) <- virus_host_transmission_reconstructed$easy_name[match(row.names(co_transmission_bacteria), virus_host_transmission_reconstructed$Bacterium)]
corrplot(as.matrix(co_transmission_bacteria), p.mat=as.matrix(co_transmission_bacpval), insig='blank', na.label = "X")


########  partial mantel test

virusN <- virus[["LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609"]]

# normalizing the distances:
virusN <- virusN/median(unname(unlist(virusN)))

# sample names parsing (because in some cases virus strain was reconstructed both in MGS and VLP simultaneously):
sample_names <- data.frame(colnames(virusN))
colnames(sample_names) <- 'sample'
sample_names$Var1 <- gsub("VLP_|MGS_", "", sample_names$sample)
sample_names$source <- NA
sample_names[grep('VLP', sample_names$sample),]$source <- 'VLP'
sample_names[grep('MGS', sample_names$sample),]$source <- 'MGS'
repeated <- data.frame(table(gsub("VLP_|MGS_", "", colnames(virusN))))
sample_names <- merge(sample_names, repeated, by='Var1')
sample_names <- sample_names[sample_names$Freq==1 | sample_names$source=='VLP',]
sample_names$bacterial <- gsub('VLP', 'MGS', sample_names$sample)


  bacteriumN <- bacterium_to_test[["SGB1836_group"]]
  # normalizing the distance:
  bacteriumN <- bacteriumN/median(unname(unlist(bacteriumN)), na.rm=T)
  
  sample_names_bacterial <- data.frame(colnames(bacteriumN))
  colnames(sample_names_bacterial) <- 'bacterial'
  
  sample_names_concurrent <- merge(sample_names, sample_names_bacterial, by='bacterial')
  
  
  # overalpping samples where both virus and bacterial strains were reconstructed:
  virusN_concurrent <- virusN[sample_names_concurrent$sample, sample_names_concurrent$sample, drop=F]
  bacteriumN <- bacteriumN[sample_names_concurrent$bacterial, sample_names_concurrent$bacterial, drop=F]
  
  dist_matrix_contol <- virusN_concurrent
  for (i in colnames(dist_matrix_contol)) {
    
    for (j in row.names(dist_matrix_contol)) {
      
      if ( substr(i, 1, 16)==substr(j, 1, 16) ) {
        dist_matrix_contol[i,j] <- 1
      } else {
        dist_matrix_contol[i,j] <- 0
      }
      
    }
    
  }
 

 mantel.partial(virusN_concurrent, bacteriumN, dist_matrix_contol, method = "spearman", permutations = 1000, na.rm = T)
 mantel(virusN_concurrent, bacteriumN, method = "pearson", permutations = 1000, na.rm = T)

 plot(unname(unlist(virusN_concurrent)), unname(unlist(bacteriumN)))


summary(lm(unname(unlist(virusN_concurrent))~ unname(unlist(bacteriumN))))



#### calculating correlation between transmitted viruses and their predicted hosts:
co_transmission_VB <- data.frame(matrix(NA, ncol=length(unique(virus_host_transmission_reconstructed$Virus)),
                                         nrow=length(unique(virus_host_transmission_reconstructed$Bacterium)))  )
colnames(co_transmission_VB) <- unique(virus_host_transmission_reconstructed$Virus)
row.names(co_transmission_VB) <- unique(virus_host_transmission_reconstructed$Bacterium)

co_transmission_pvalVB <- co_transmission_VB

mantel_partial_results <- list()
for (n in virus_host_transmission_reconstructed$Virus) {
  
  # virus distance matrix:
  virusN <- virus[[n]]
  
  # normalizing the distances:
  virusN <- virusN/median(unname(unlist(virusN)))
  
  # sample names parsing (because in some cases virus strain was reconstructed both in MGS and VLP simultaneously):
  sample_names <- data.frame(colnames(virusN))
  colnames(sample_names) <- 'sample'
  sample_names$Var1 <- gsub("VLP_|MGS_", "", sample_names$sample)
  sample_names$source <- NA
  sample_names[grep('VLP', sample_names$sample),]$source <- 'VLP'
  sample_names[grep('MGS', sample_names$sample),]$source <- 'MGS'
  repeated <- data.frame(table(gsub("VLP_|MGS_", "", colnames(virusN))))
  sample_names <- merge(sample_names, repeated, by='Var1')
  sample_names <- sample_names[sample_names$Freq==1 | sample_names$source=='VLP',]
  sample_names$bacterial <- gsub('VLP', 'MGS', sample_names$sample)
  
  # in some cases, multiple hosts are predicted for 1 virus: 
  for (b in virus_host_transmission_reconstructed$Bacterium) {
    
    bacteriumN <- bacterium_to_test[[b]]
    # normalizing the distance:
    bacteriumN <- bacteriumN/median(unname(unlist(bacteriumN)), na.rm=T)
    
    sample_names_bacterial <- data.frame(colnames(bacteriumN))
    colnames(sample_names_bacterial) <- 'bacterial'
    
    sample_names_concurrent <- merge(sample_names, sample_names_bacterial, by='bacterial')
    
    # overalpping samples where both virus and bacterial strains were reconstructed:
    virusN_concurrent <- virusN[sample_names_concurrent$sample, sample_names_concurrent$sample, drop=F]
    bacteriumN <- bacteriumN[sample_names_concurrent$bacterial, sample_names_concurrent$bacterial, drop=F]
    
    # since 1000 permutations only make sense when there are enough samples to permute, the cut-off is 7, as 7! > 1000
    if (length(sample_names_concurrent$bacterial)>=7 & length(sample_names_concurrent$sample)>=7) {
      
      dist_matrix_contol <- virusN_concurrent
      for (i in colnames(dist_matrix_contol)) {
        
        for (j in row.names(dist_matrix_contol)) {
          
          if ( substr(i, 1, 16)==substr(j, 1, 16) ) {
            dist_matrix_contol[i,j] <- 1
          } else {
            dist_matrix_contol[i,j] <- 0
          }
          
        }
        
      }
      
      
      mantel_partial_results[[paste0(n, '_', b)]] <- mantel.partial(virusN_concurrent, bacteriumN, dist_matrix_contol, method = "pearson", permutations = 999, na.rm = T)
      
      co_transmission_VB[b,n] <- mantel_partial_results[[paste0(n, '_', b)]]$statistic
      co_transmission_pvalVB[b,n] <- mantel_partial_results[[paste0(n, '_', b)]]$signif
      
    } else {
      mantel_partial_results[[paste0(n, '_', b)]] <- 'No concurrent samples'
    }
    
    
    
  }
  
}

keep_co_transmission_VB <- co_transmission_VB
keep_co_transmission_pvalVB <- co_transmission_pvalVB

co_transmission_VB <- as.matrix(co_transmission_VB)
co_transmission_VB[is.nan(co_transmission_VB)] <- NA
co_transmission_pvalVB <- as.matrix(co_transmission_pvalVB)

co_transmission_VB[co_transmission_pvalVB>0.05] <- 0

corrplot(co_transmission_VB)

##############################
# OUTPUT
##############################
write.table(co_transmission_cor, '03a.RESULTS/Co_transmission.txt', sep='\t', quote=F)
write.table(bgcolors, '03a.RESULTS/Co_transmission_colors.txt', sep='\t', quote=F)

