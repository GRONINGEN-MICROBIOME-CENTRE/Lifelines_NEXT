setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore inter-individual variation of bacterial
# strains and then derive cut-offs for sharing events
# between maternal and infant samples
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(vegan)
library(ggplot2)

library(ggtree)
library(scales)
library(ggforestplot)
library(ggforce)

library(cutpointr)
library(patchwork)
library(reshape2)

##############################
# Input data
##############################

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$source <- "MGS"
MGS_metadata$Alex_ID <- paste0(MGS_metadata$FAM_ID, '_', MGS_metadata$Type, '_', substr(MGS_metadata$Short_sample_ID_bact, 1,1), '_', MGS_metadata$source, '_', MGS_metadata$Timepoint)

folder = "01.RAW_DATA/strainphlan_4_distance_matrix/"  
file_list = list.files(path=folder, pattern="*.csv")  
bacterium=lapply(paste0(folder, file_list), function(x) read.csv(x, header=T))
names(bacterium) <- gsub("_dmat_.csv", '', file_list)

## since my parsing scripts have certain ID format: 

bacterium <- lapply(bacterium, function(x) {
  row.names(x) <- substr(row.names(x), 0, 12) # first remove data processing artifacts
  colnames(x) <- substr(colnames(x), 0, 12) 
  
  row.names(x) <- MGS_metadata$Alex_ID[match(row.names(x), MGS_metadata$NG_ID)]
  colnames(x) <- MGS_metadata$Alex_ID[match(colnames(x), MGS_metadata$NG_ID)] # then change the sequencing IDs to IDs for parsing
  x
})

bacterium_to_test <- list()
for (n in 1:NROW(bacterium)) {
  
  bacteriumN <- bacterium[[n]]
  bacteriumName <- names(bacterium[n])
  
  if ( length(grep('Infant', colnames(bacteriumN)))!=0 & length(grep('Mother', colnames(bacteriumN)))!=0) {
    bacterium_to_test[[bacteriumName]] <- bacterium[[bacteriumName]]
  }
  
}

species_names <- read.table('02.CLEAN_DATA/List_bacterial_strains_reconstructed_and_tested.txt', sep='\t', header=T)

list_transmitted <- read.table('02.CLEAN_DATA/List_bacteria_results_checking_transmission.txt', sep='\t', header=T)
list_transmitted <- list_transmitted[list_transmitted$FDR<= 0.05,]

bacterium_to_test <- bacterium_to_test[names(bacterium_to_test) %in% list_transmitted$Bacterium]

host_assignment <- read.table('02.CLEAN_DATA/Table_predicted_hosts_shared_strains_maximized_Youden_combined_wilcox_less.txt', sep='\t', header=T)
colnames(host_assignment)[grep("species",colnames(host_assignment))] <- "Species"
##############################
# ANALYSIS
##############################

# defining the inter-personal bacterium genetic distance variation

within_individual_distances_bacterium_infants <- list()
between_individual_distances_bacterium_infants <- list()

within_individual_distances_bacterium_mothers <- list()
between_individual_distances_bacterium_mothers <- list()

for (n in 1:NROW(bacterium_to_test)) {
  print(paste0(' > prepping within-infant distances for strain ', names(bacterium_to_test[n])))  
  
  # get the virus name and distance matrix
  bacteriumN <- bacterium_to_test[[n]]
  bacteriumN_infants <- bacteriumN[grep('Infant', row.names(bacteriumN)), grep('Infant', colnames(bacteriumN))]
  bacteriumN_mothers <- bacteriumN[grep('Mother', row.names(bacteriumN)), grep('Mother', colnames(bacteriumN))]
  
  bacteriumName <- names(bacterium_to_test[n])
  
  within_individual_distances <- list()
  between_individual_distances <- list()
  
  for (i in unique( gsub("(FAM0[0-9]{3}_Infant_[A-Za-z]).*", "\\1", colnames(bacteriumN_infants)[grep('Infant', colnames(bacteriumN_infants))] ) ) ) {
    
    if (length(grep(i, colnames(bacteriumN_infants)))>1) {
      
      A  <- bacteriumN_infants[grep(i, colnames(bacteriumN_infants)),grep(i, colnames(bacteriumN_infants))]
      within_individual_distances[[i]] <- A[upper.tri(A)]
      between_individual_distances[[i]] <- unlist(bacteriumN_infants[grep(i, colnames(bacteriumN_infants)),grep(i, colnames(bacteriumN_infants), invert=T)])
    } 
    
  }
  
  within_individual_distances_bacterium_infants[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(within_individual_distances) )))
  between_individual_distances_bacterium_infants[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(between_individual_distances) )))
  
  
  within_individual_distances <- list()
  between_individual_distances <- list()
  
  for (i in unique( gsub("(FAM0[0-9]{3}_Mother_[A-Za-z]).*", "\\1", colnames(bacteriumN_mothers)[grep('Mother', colnames(bacteriumN_mothers))] ) ) ) {
    
    if (length(grep(i, colnames(bacteriumN_mothers)))>1) {
      
      A  <- bacteriumN_mothers[grep(i, colnames(bacteriumN_mothers)),grep(i, colnames(bacteriumN_mothers))]
      within_individual_distances[[i]] <- A[upper.tri(A)]
      between_individual_distances[[i]] <- unlist(bacteriumN_mothers[grep(i, colnames(bacteriumN_mothers)),grep(i, colnames(bacteriumN_mothers), invert=T)])
    } 
    
  }
  
  within_individual_distances_bacterium_mothers[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(within_individual_distances) )))
  between_individual_distances_bacterium_mothers[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(between_individual_distances) )))
  
  
  
}

### Threshold calculation in infants 
within_individual_distances_bacterium_infants[["SGB15323"]] <- NULL
within_individual_distances_bacterium_infants[["SGB15342"]] <- NULL
within_individual_distances_bacterium_infants[["SGB2290"]] <- NULL

between_individual_distances_bacterium_infants <- between_individual_distances_bacterium_infants[names(between_individual_distances_bacterium_infants) %in% names(within_individual_distances_bacterium_infants)]

cutpoints_infants <- data.frame( matrix(NA, nrow = length(within_individual_distances_bacterium_infants), ncol=3 )  )
row.names(cutpoints_infants) <- names(within_individual_distances_bacterium_infants)
colnames(cutpoints_infants) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

infants_bacterium_hists_data <- list()

for (n in 1:NROW(within_individual_distances_bacterium_infants)) {
  
  bacteriumName <- names(within_individual_distances_bacterium_infants[n])
  
  infants_bacterium <- data.frame( c(within_individual_distances_bacterium_infants[[n]],
                                 between_individual_distances_bacterium_infants[[n]]),
                               c(rep("Within", length(within_individual_distances_bacterium_infants[[n]])),
                                 rep("Between", length(between_individual_distances_bacterium_infants[[n]])) ) ) 
  
  colnames(infants_bacterium) <- c('Distance', 'Variable')
  
  infants_bacterium$Distance <- infants_bacterium$Distance/median(unname(unlist(bacterium_to_test[[bacteriumName]][grep('Infant', colnames(bacterium_to_test[[bacteriumName]])), grep('Infant', row.names(bacterium_to_test[[bacteriumName]])) ])), na.rm=T)
  
  Youden <- cutpointr(infants_bacterium, Distance, Variable, 
                      pos_class = "Within",
                      neg_class = "Between", 
                      method = maximize_metric, 
                      metric = youden, 
                      na.rm=T)
  cutpoints_infants[n,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_infants[n,"FDR_value"] <- quantile(infants_bacterium[infants_bacterium$Variable=="Between",]$Distance, probs = 0.05)
  
  cutpoints_infants[n,"N_within_comparisons"] <- length(within_individual_distances_bacterium_infants[[n]])
  
  infants_bacterium_hists_data[[bacteriumName]] <- infants_bacterium
}

distance_histograms <- list()

for (n in 1:NROW(infants_bacterium_hists_data)) {
  
  distance_histograms[[names(infants_bacterium_hists_data[n])]] <- ggplot(infants_bacterium_hists_data[[n]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(infants_bacterium_hists_data[[n]]$Distance) ), alpha=0.7) + 
    geom_density(aes(y = (after_stat(count)/sum(after_stat(count)))*1000), alpha=0.2   ) + 
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=cutpoints_infants[ names(infants_bacterium_hists_data[n])  ,"Youden_index"]), linetype="dashed") + 
    geom_vline(aes(xintercept=cutpoints_infants[ names(infants_bacterium_hists_data[n])  ,"FDR_value"]), linetype="dashed", color="red") +
    ggtitle( gsub("\\...\\K\\d+", "", gsub(".*_length", "L"  , names(infants_bacterium_hists_data[n]) ), perl=T) ) +
    theme_bw() + 
    theme(title = element_text(size=8))
  
}


combined_plot <- distance_histograms[[1]]

for (h in 2:NROW(infants_bacterium_hists_data)) {
  combined_plot <- combined_plot + distance_histograms[[h]]
}


combined_plot <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "")
pdf('./04.PLOTS/Cutpoints_within_infants_bacterium.pdf', width=27/2.54, height=24/2.54)
combined_plot
dev.off()

#### Threshold in mothers
within_individual_distances_bacterium_mothers[["SGB17247"]] <- NULL

between_individual_distances_bacterium_mothers <- between_individual_distances_bacterium_mothers[names(between_individual_distances_bacterium_mothers) %in% names(within_individual_distances_bacterium_mothers)]

cutpoints_mothers <- data.frame( matrix(NA, nrow = length(within_individual_distances_bacterium_mothers), ncol=3 )  )
row.names(cutpoints_mothers) <- names(within_individual_distances_bacterium_mothers)
colnames(cutpoints_mothers) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

mothers_bacterium_hists_data <- list()
for (n in 1:NROW(within_individual_distances_bacterium_mothers)) {
  
  bacteriumName <- names(within_individual_distances_bacterium_mothers[n])
  
  mothers_bacterium <- data.frame( c(within_individual_distances_bacterium_mothers[[n]],
                                 between_individual_distances_bacterium_mothers[[n]]),
                               c(rep("Within", length(within_individual_distances_bacterium_mothers[[n]])),
                                 rep("Between", length(between_individual_distances_bacterium_mothers[[n]])) ) ) 
  
  colnames(mothers_bacterium) <- c('Distance', 'Variable')
  
  

    mothers_bacterium$Distance <- mothers_bacterium$Distance/median(unname(unlist(bacterium_to_test[[bacteriumName]][grep('Mother', colnames(bacterium_to_test[[bacteriumName]])), grep('Mother', row.names(bacterium_to_test[[bacteriumName]])) ])), na.rm=T)

  
  
  Youden <- cutpointr(mothers_bacterium, Distance, Variable, 
                      pos_class = "Within",
                      neg_class = "Between", 
                      method = maximize_metric, 
                      metric = youden, 
                      na.rm=T)
  cutpoints_mothers[n,"Youden_index"] <- Youden$optimal_cutpoint
  
  if (round(length(mothers_bacterium[mothers_bacterium$Variable=="Between",]$Distance)*0.05)==0) {
    cutpoints_mothers[n,"FDR_value"] <- sort(mothers_bacterium[mothers_bacterium$Variable=="Between",]$Distance)[1]
  } else {
    cutpoints_mothers[n,"FDR_value"] <- sort(mothers_bacterium[mothers_bacterium$Variable=="Between",]$Distance)[ round(length(mothers_bacterium[mothers_bacterium$Variable=="Between",]$Distance)*0.05) ]
  }
  
  cutpoints_mothers[n,"N_within_comparisons"] <- length(within_individual_distances_bacterium_mothers[[n]])
  
  mothers_bacterium_hists_data[[names(within_individual_distances_bacterium_mothers[n])]] <- mothers_bacterium
}

distance_histograms_mothers <- list()

for (n in 1:NROW(mothers_bacterium_hists_data)) {
  
  distance_histograms_mothers[[names(mothers_bacterium_hists_data[n])]] <- ggplot(mothers_bacterium_hists_data[[n]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(mothers_bacterium_hists_data[[n]]$Distance) ), alpha=0.7) + 
    geom_density(aes(y = (after_stat(count)/sum(after_stat(count)))*2000), alpha=0.2   ) + 
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=cutpoints_mothers[ names(mothers_bacterium_hists_data[n])  ,"Youden_index"]), linetype="dashed") + 
    geom_vline(aes(xintercept=cutpoints_mothers[ names(mothers_bacterium_hists_data[n])  ,"FDR_value"]), linetype="dashed", color="red") +
    ggtitle( gsub("\\...\\K\\d+", "", gsub(".*_length", "L"  , names(mothers_bacterium_hists_data[n]) ), perl=T) ) +
    theme_bw() + 
    theme(title=element_text(size=8))
  
}


combined_plot_mothers <- distance_histograms_mothers[[1]]

for (h in 2:NROW(mothers_bacterium_hists_data)) {
  combined_plot_mothers <- combined_plot_mothers + distance_histograms_mothers[[h]]
}

combined_plot_mothers <- combined_plot_mothers +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "")

pdf('./04.PLOTS/Cutpoints_within_mothers_bacteria.pdf', width=27/2.54, height=24/2.54)
combined_plot_mothers
dev.off()

#### Threshold merged: 

cutpoints_all <- data.frame( matrix(NA, nrow = length(bacterium_to_test), ncol=3 )  )
row.names(cutpoints_all) <- names(bacterium_to_test)
colnames(cutpoints_all) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

bacterium_hists_data <- list()
for (n in 1:NROW(bacterium_to_test)) {
  
  bacteriumN <- data.frame( c(within_individual_distances_bacterium_infants[[n]],
                          between_individual_distances_bacterium_infants[[n]],
                          within_individual_distances_bacterium_mothers[[n]],
                          between_individual_distances_bacterium_mothers[[n]]),
                        c(rep("Within", length(within_individual_distances_bacterium_infants[[n]])),
                          rep("Between", length(between_individual_distances_bacterium_infants[[n]])),
                          rep("Within", length(within_individual_distances_bacterium_mothers[[n]])),
                          rep("Between", length(between_individual_distances_bacterium_mothers[[n]])) ) ) 
  
  colnames(bacteriumN) <- c('Distance', 'Variable')
  
  bacteriumN <- bacteriumN[bacteriumN$Distance!=Inf,]
  
  bacteriumN$Distance <- bacteriumN$Distance/median( unname(unlist(bacterium_to_test[[n]])), na.rm = T )
  
  Youden <- cutpointr(bacteriumN, Distance, Variable, 
                      pos_class = "Within",
                      neg_class = "Between", 
                      method = oc_youden_kernel, 
                      metric = youden, 
                      na.rm=T)
  cutpoints_all[n,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_all[n,"FDR_value"] <- quantile(bacteriumN[bacteriumN$Variable=="Between",]$Distance, probs = 0.05)
  
  cutpoints_all[n,"N_within_comparisons"] <- length(c(within_individual_distances_bacterium_infants[[n]], within_individual_distances_bacterium_mothers[[n]]))
  
  bacterium_hists_data[[names(bacterium[n])]] <- bacteriumN
}

distance_histograms_all <- list()

for (n in 1:NROW(bacterium_hists_data)) {
  
  distance_histograms_all[[names(bacterium_hists_data[n])]] <- ggplot(bacterium_hists_data[[n]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(bacterium_hists_data[[n]]$Distance) ), alpha=0.7) + 
    geom_density(aes(y = (after_stat(count)/sum(after_stat(count)))*2000), alpha=0.2   ) + 
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=cutpoints_all[ names(bacterium_hists_data[n])  ,"Youden_index"]), linetype="dashed") + 
    geom_vline(aes(xintercept=cutpoints_all[ names(bacterium_hists_data[n])  ,"FDR_value"]), linetype="dashed", color="red") +
    ggtitle( gsub("\\...\\K\\d+", "", gsub(".*_length", "L"  , names(bacterium_hists_data[n]) ), perl=T) ) +
    theme_bw() + 
    theme(title=element_text(size=8))
  
}


combined_plot <- distance_histograms_all[[1]]

for (h in 2:NROW(bacterium_hists_data)) {
  combined_plot <- combined_plot + distance_histograms_all[[h]]
}

combined_plot <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "")

pdf('./04.PLOTS/Cutpoints_within_mothers_and_infants_combined_bacterium.pdf', width=27/2.54, height=24/2.54)
combined_plot
dev.off()

list_transmitted$cutpoint <- NA

for (i in list_transmitted$Bacterium) {
  list_transmitted[list_transmitted$Bacterium==i,]$cutpoint <- ifelse( cutpoints_all[i,"Youden_index"] >= cutpoints_all[i,"FDR_value"],  cutpoints_all[i,"FDR_value"],  cutpoints_all[i,"Youden_index"])
}

# based on the inter-individual strain variation, how often is the bacterium transmitted between related and unrelated individuals?

transmission_bacterium <- bacterium_to_test[names(bacterium_to_test) %in% list_transmitted$Bacterium]
mother_infant_distances_bacterium <- list()
unrelated_distances_bacterium <- list()

list_transmitted$N_related_pairs_transmitted <- NA
list_transmitted$N_unrelated_pairs_transmitted <- NA
list_transmitted$Perc_related_pairs_transmitted <- NA
list_transmitted$Perc_unrelated_pairs_transmitted <- NA
list_transmitted$N_related_pairs <- NA
list_transmitted$N_unrelated_pairs <- NA

for (n in 1:NROW(transmission_bacterium)) {
  # get the virus name and distance matrix
  
  bacteriumN <- transmission_bacterium[[n]]
  bacteriumN <- bacteriumN/median( unname(unlist(transmission_bacterium[[n]])), na.rm = T )
  bacteriumName <- names(transmission_bacterium[n])
  
  bacteriumN[bacteriumN <= list_transmitted[list_transmitted$Bacterium==bacteriumName,]$cutpoint ] <- -1
  bacteriumN <- bacteriumN[grep('Mother', row.names(bacteriumN)), grep('Infant', colnames(bacteriumN))]
  
  # two lists to store pair-wise distances for the current virus[[n]] per family
  mother_infant_distances <- list()
  mother_unrelated_distances <- list()
  
  # loop over all families
  for (i in unique( substr(row.names( bacteriumN ), 1, 7)) ) { # getting the FAM_IDs of all virus-positive families from the distance matrix 
    
    # Check if within the selected family both mother and infant have virus-positive samples
    if (any(grep(paste0(i, '_Mother'),  row.names(bacteriumN) ))==T & any(grep(paste0(i, '_Infant'),  colnames(bacteriumN) ))==T) {
      
      # Parsing pair-wise distances between maternal and infant samples from the selected family
      mother_infant_distances[[i]] <- unlist( bacteriumN[grep( paste0(i, '_Mother') , row.names(bacteriumN) ), grep( paste0(i, '_Infant'), colnames( bacteriumN ) ) ] )
      
      # Parsing pair_wise distances between maternal samples and samples of unrelated individuals:
      # Assigning all pair-wise distances between samples of related individuals as NA for the selected family
      
      mother_unrelated_distances[[i]] <- unlist( bacteriumN[grep( paste0(i, '_Mother') , row.names(bacteriumN) ), 
                                                        grep( paste0(i, '_Infant'), colnames( bacteriumN ), invert = T ) ] )
      
    } else {
      
      mother_infant_distances[[i]] <- NA
      mother_unrelated_distances[[i]] <- NA
      
    }
    
  }
  
  # collecting the distances per virus:
  mother_infant_distances_bacterium[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(mother_infant_distances) )))
  unrelated_distances_bacterium[[bacteriumName]] <- as.numeric(na.omit( unname( unlist(mother_unrelated_distances) ) ) )
  
  list_transmitted[list_transmitted$Bacterium==bacteriumName,]$N_related_pairs_transmitted <- sum(mother_infant_distances_bacterium[[bacteriumName]]==-1)
  list_transmitted[list_transmitted$Bacterium==bacteriumName,]$N_unrelated_pairs_transmitted <- sum(unrelated_distances_bacterium[[bacteriumName]]==-1)
  list_transmitted[list_transmitted$Bacterium==bacteriumName,]$Perc_related_pairs_transmitted <- sum(mother_infant_distances_bacterium[[bacteriumName]]==-1)/length(mother_infant_distances_bacterium[[bacteriumName]])
  list_transmitted[list_transmitted$Bacterium==bacteriumName,]$Perc_unrelated_pairs_transmitted <- sum(unrelated_distances_bacterium[[bacteriumName]]==-1)/length(unrelated_distances_bacterium[[bacteriumName]])
  list_transmitted[list_transmitted$Bacterium==bacteriumName,]$N_related_pairs <- length(mother_infant_distances_bacterium[[bacteriumName]])
  list_transmitted[list_transmitted$Bacterium==bacteriumName,]$N_unrelated_pairs <- length(unrelated_distances_bacterium[[bacteriumName]])
}

row.names(list_transmitted) <- list_transmitted$Bacterium
list_transmitted$p_value_kinship_transmission <- NA

testing_enrichment_transmission <- list()
for (i in list_transmitted$Bacterium) {
  
  testing_enrichment_transmission[[i]] <- matrix( c(list_transmitted[i,"N_related_pairs_transmitted"], 
                                                    (list_transmitted[i, "N_related_pairs"] - list_transmitted[i,"N_related_pairs_transmitted"]),
                                                    list_transmitted[i,"N_unrelated_pairs_transmitted"],
                                                    (list_transmitted[i, "N_unrelated_pairs"] - list_transmitted[i,"N_unrelated_pairs_transmitted"])),
                                                  nrow=2,
                                                  dimnames = list(Transmitted=c('Related', 'Unrelated'),
                                                                  Not_transmitted=c('Related', 'Unrelated'))) 
  
  check_association_with_kinship <- fisher.test(testing_enrichment_transmission[[i]])
  
  list_transmitted[i,]$p_value_kinship_transmission <- check_association_with_kinship$p.value
  
}

list_transmitted$p_value_kinship_transmission_adj <- p.adjust(list_transmitted$p_value_kinship_transmission, method = "BH")

list_transmitted$significance_level_transmission <- NA
list_transmitted[list_transmitted$p_value_kinship_transmission_adj > 0.05,]$significance_level_transmission <- 'ns'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.05,]$significance_level_transmission <- '*'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.01,]$significance_level_transmission <- '**'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.001,]$significance_level_transmission <- '***'
list_transmitted$easy_name <- paste0(list_transmitted$Bacterium, '_', list_transmitted$Species)

for_plot <- melt(list_transmitted[,c(17,10:11,16)])

pdf('./04.PLOTS/Perc_transmitted_in_pairs_maximized_Youden_bacteria.pdf', width=20.7/2.54, height=24/2.54)
ggplot(for_plot, aes(value, easy_name, fill=variable)) + 
  geom_col(position = 'dodge', alpha=0.8) +
  labs (y="Bacteria", x="% transmitted in mother-infant samples pairs") + 
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(fill="Kinship") + 
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  geom_text(aes(label = significance_level_transmission, x = 1.05, y = easy_name), size = 3, angle=270) 
dev.off()


virus_host_transmission <- merge(host_assignment, list_transmitted, by='Species', all=T)

##############################
# OUTPUT
##############################
write.table(virus_host_transmission, "03a.RESULTS/Bacterial_hosts_transmission_for_transmitted_viruses.txt", sep='\t', row.names = F, quote=F)
