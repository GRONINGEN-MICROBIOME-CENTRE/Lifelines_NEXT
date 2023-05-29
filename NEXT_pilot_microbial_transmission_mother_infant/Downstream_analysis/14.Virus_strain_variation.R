setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore inter-individual variation of virus 'strains' 
# in both VLP and MGS samples
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(vegan)
library(ggplot2)

library(ape)
library(treeio)

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

metadata <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/metadata_combined_for_exp.txt', sep='\t', header=T)
metadata$Alex_ID <- paste0(metadata$FAM_ID, '_', metadata$Type, '_', substr(metadata$Short_sample_ID, 1,1), '_', metadata$source, '_', metadata$Timepoint)
row.names(metadata) <- metadata$Alex_ID
metadata$label <- metadata$Alex_ID
  
RPKM_combined_0.95 <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/RPKM_counts_combined_0.95_UPD_final_for_exp.txt', sep='\t', header=T)

selected_viruses <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/Viruses_shared_min_5_families_UPD_final.txt', sep='\t', header=T)
colnames(selected_viruses)[1] <- 'Virus'

list_transmitted <- read.table('02.CLEAN_DATA/List_viruses_results_checking_transmission.txt', sep='\t', header=T)
list_transmitted <- list_transmitted[list_transmitted$FDR<= 0.05,]
list_transmitted <- merge(list_transmitted, selected_viruses, by='Virus')

folder = "02.CLEAN_DATA/dist.matrices/"  
file_list = list.files(path=folder, pattern="*.txt")  
virus=lapply(paste0(folder, file_list), function(x) read.table(x, sep='\t', header=T))
names(virus) <- gsub(".dist.txt", '', file_list)

virus <- virus[names(virus) %in% list_transmitted$Virus]

host_assignment <- read.csv('01.RAW_DATA/Host_assignment/MERGED_Host_prediction_to_genome_m90.csv')

hosts <- read.csv('01.RAW_DATA/Host_assignment/MERGED_Host_prediction_to_genome_m90.csv')

metaphlan <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)
metaphlan_species <- metaphlan[grep('s__', row.names(metaphlan)),]


##############################
# ANALYSIS
##############################

# defining the inter-personal virus genetic distance variation

within_individual_distances_virus_infants <- list()
between_individual_distances_virus_infants <- list()

within_individual_distances_virus_mothers <- list()
between_individual_distances_virus_mothers <- list()


for (n in 1:NROW(virus)) {
  print(paste0(' > prepping within-infant distances for strain ', names(virus[n])))  
  
  # get the virus name and distance matrix
  virusN <- virus[[n]]
  virusN_infants <- virusN[grep('Infant', row.names(virusN)), grep('Infant', colnames(virusN))]
  virusN_mothers <- virusN[grep('Mother', row.names(virusN)), grep('Mother', colnames(virusN))]
  
  virusName <- names(virus[n])
  
  within_individual_distances <- list()
  between_individual_distances <- list()
  
  for (i in unique( gsub("(FAM0[0-9]{3}_Infant_[A-Za-z]).*", "\\1", colnames(virusN_infants)[grep('Infant', colnames(virusN_infants))] ) ) ) {
    
    if (length(grep(i, colnames(virusN_infants)))>1) {
      
      A  <- virusN_infants[grep(i, colnames(virusN_infants)),grep(i, colnames(virusN_infants))]
      within_individual_distances[[i]] <- A[upper.tri(A)]
      between_individual_distances[[i]] <- unlist(virusN_infants[grep(i, colnames(virusN_infants)),grep(i, colnames(virusN_infants), invert=T)])
    } 
    
  }
  
  within_individual_distances_virus_infants[[virusName]] <- as.numeric(na.omit(unname ( unlist(within_individual_distances) )))
  between_individual_distances_virus_infants[[virusName]] <- as.numeric(na.omit(unname ( unlist(between_individual_distances) )))
  
  
  within_individual_distances <- list()
  between_individual_distances <- list()
  
  for (i in unique( gsub("(FAM0[0-9]{3}_Mother_[A-Za-z]).*", "\\1", colnames(virusN_mothers)[grep('Mother', colnames(virusN_mothers))] ) ) ) {
    
    if (length(grep(i, colnames(virusN_mothers)))>1) {
      
      A  <- virusN_mothers[grep(i, colnames(virusN_mothers)),grep(i, colnames(virusN_mothers))]
      within_individual_distances[[i]] <- A[upper.tri(A)]
      between_individual_distances[[i]] <- unlist(virusN_mothers[grep(i, colnames(virusN_mothers)),grep(i, colnames(virusN_mothers), invert=T)])
    } 
    
  }
  
  within_individual_distances_virus_mothers[[virusName]] <- as.numeric(na.omit(unname ( unlist(within_individual_distances) )))
  between_individual_distances_virus_mothers[[virusName]] <- as.numeric(na.omit(unname ( unlist(between_individual_distances) )))
  
  
  
}

### Threshold calculation

cutpoints_infants <- data.frame( matrix(NA, nrow = length(virus), ncol=3 )  )
row.names(cutpoints_infants) <- names(virus)
colnames(cutpoints_infants) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

infants_virus_hists_data <- list()
for (n in 1:NROW(virus)) {
  
  infants_virus <- data.frame( c(within_individual_distances_virus_infants[[n]],
                                 between_individual_distances_virus_infants[[n]]),
                               c(rep("Within", length(within_individual_distances_virus_infants[[n]])),
                                 rep("Between", length(between_individual_distances_virus_infants[[n]])) ) ) 
  
  colnames(infants_virus) <- c('Distance', 'Variable')
  
  infants_virus$Distance <- infants_virus$Distance/median( unname(unlist(virus[[n]])) )
  
  Youden <- cutpointr(infants_virus, Distance, Variable, 
                            pos_class = "Within",
                            neg_class = "Between", 
                            method = maximize_metric, 
                            metric = youden)
  cutpoints_infants[n,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_infants[n,"FDR_value"] <- quantile(infants_virus[infants_virus$Variable=="Between",]$Distance, probs = 0.05)
  
  cutpoints_infants[n,"N_within_comparisons"] <- length(within_individual_distances_virus_infants[[n]])
  
  infants_virus_hists_data[[names(virus[n])]] <- infants_virus
}

distance_histograms <- list()

for (n in 1:NROW(infants_virus_hists_data)) {
  
  distance_histograms[[names(infants_virus_hists_data[n])]] <- ggplot(infants_virus_hists_data[[n]], aes(x=Distance, fill=Variable)) + 
                                  geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(infants_virus_hists_data[[n]]$Distance) ), alpha=0.7) + 
                                  geom_density(aes(y = (after_stat(count)/sum(after_stat(count)))*2000), alpha=0.2   ) + 
                                  labs(x="Normalized Distance", y="proportion (%)") +
                                  geom_vline(aes(xintercept=cutpoints_infants[ names(infants_virus_hists_data[n])  ,"Youden_index"]), linetype="dashed") + 
                                  geom_vline(aes(xintercept=cutpoints_infants[ names(infants_virus_hists_data[n])  ,"FDR_value"]), linetype="dashed", color="red") +
                                  ggtitle( gsub("\\...\\K\\d+", "", gsub(".*_length", "L"  , names(infants_virus_hists_data[n]) ), perl=T) ) +
                                  theme_bw() + 
                                  theme(title = element_text(size=8))
  
}


combined_plot <- distance_histograms[[1]]

for (h in 2:NROW(infants_virus_hists_data)) {
  combined_plot <- combined_plot + distance_histograms[[h]]
}


combined_plot <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "")
pdf('./04.PLOTS/Cutpoints_within_infants_viruses.pdf', width=27/2.54, height=24/2.54)
combined_plot
dev.off()


#### Threshold in mothers

# removing those that do not have between-individual distances
between_individual_distances_virus_mothers[["LN_4E02_VL_255_NODE_333_length_8047_cov_10.923048"]] <- NULL
between_individual_distances_virus_mothers[["LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438"]] <- NULL

within_individual_distances_virus_mothers[["LN_4E02_VL_255_NODE_333_length_8047_cov_10.923048"]] <- NULL
within_individual_distances_virus_mothers[["LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438"]] <- NULL

cutpoints_mothers <- data.frame( matrix(NA, nrow = length(within_individual_distances_virus_mothers), ncol=3 )  )
row.names(cutpoints_mothers) <- names(within_individual_distances_virus_mothers)
colnames(cutpoints_mothers) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

mothers_virus_hists_data <- list()
for (n in 1:NROW(within_individual_distances_virus_mothers)) {
  
  mothers_virus <- data.frame( c(within_individual_distances_virus_mothers[[n]],
                                 between_individual_distances_virus_mothers[[n]]),
                               c(rep("Within", length(within_individual_distances_virus_mothers[[n]])),
                                 rep("Between", length(between_individual_distances_virus_mothers[[n]])) ) ) 
  
  colnames(mothers_virus) <- c('Distance', 'Variable')
  
  
  if (median( unname(unlist(within_individual_distances_virus_mothers[[n]])) ) ==0) {
    mothers_virus$Distance <- mothers_virus$Distance/sort(unname(unlist(within_individual_distances_virus_mothers[[n]])) [unname(unlist(within_individual_distances_virus_mothers[[n]]))!=0])[1]
  } else {
    mothers_virus$Distance <- mothers_virus$Distance/median( unname(unlist(within_individual_distances_virus_mothers[[n]])) )
  }
 
  
  Youden <- cutpointr(mothers_virus, Distance, Variable, 
                      pos_class = "Within",
                      neg_class = "Between", 
                      method = maximize_metric, 
                      metric = youden)
  cutpoints_mothers[n,"Youden_index"] <- Youden$optimal_cutpoint
  
  if (round(length(mothers_virus[mothers_virus$Variable=="Between",]$Distance)*0.05)==0) {
    cutpoints_mothers[n,"FDR_value"] <- sort(mothers_virus[mothers_virus$Variable=="Between",]$Distance)[1]
  } else {
    cutpoints_mothers[n,"FDR_value"] <- sort(mothers_virus[mothers_virus$Variable=="Between",]$Distance)[ round(length(mothers_virus[mothers_virus$Variable=="Between",]$Distance)*0.05) ]
  }
  
  cutpoints_mothers[n,"N_within_comparisons"] <- length(within_individual_distances_virus_mothers[[n]])
  
  mothers_virus_hists_data[[names(within_individual_distances_virus_mothers[n])]] <- mothers_virus
}

distance_histograms_mothers <- list()

for (n in 1:NROW(mothers_virus_hists_data)) {
  
  distance_histograms_mothers[[names(mothers_virus_hists_data[n])]] <- ggplot(mothers_virus_hists_data[[n]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(mothers_virus_hists_data[[n]]$Distance) ), alpha=0.7) + 
    geom_density(aes(y = (after_stat(count)/sum(after_stat(count)))*2000), alpha=0.2   ) + 
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=cutpoints_mothers[ names(mothers_virus_hists_data[n])  ,"Youden_index"]), linetype="dashed") + 
    geom_vline(aes(xintercept=cutpoints_mothers[ names(mothers_virus_hists_data[n])  ,"FDR_value"]), linetype="dashed", color="red") +
    ggtitle( gsub("\\...\\K\\d+", "", gsub(".*_length", "L"  , names(mothers_virus_hists_data[n]) ), perl=T) ) +
    theme_bw() + 
    theme(title=element_text(size=8))
  
}


combined_plot_mothers <- distance_histograms_mothers[[1]]

for (h in 2:NROW(mothers_virus_hists_data)) {
  combined_plot_mothers <- combined_plot_mothers + distance_histograms_mothers[[h]]
}

combined_plot_mothers <- combined_plot_mothers +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "")

pdf('./04.PLOTS/Cutpoints_within_mothers_viruses.pdf', width=27/2.54, height=24/2.54)
combined_plot_mothers
dev.off()

#### Threshold merged: 

cutpoints_all <- data.frame( matrix(NA, nrow = length(virus), ncol=3 )  )
row.names(cutpoints_all) <- names(virus)
colnames(cutpoints_all) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

virus_hists_data <- list()
for (n in 1:NROW(virus)) {
  
  virusN <- data.frame( c(within_individual_distances_virus_infants[[n]],
                          between_individual_distances_virus_infants[[n]],
                          within_individual_distances_virus_mothers[[n]],
                          between_individual_distances_virus_mothers[[n]]),
                               c(rep("Within", length(within_individual_distances_virus_infants[[n]])),
                                 rep("Between", length(between_individual_distances_virus_infants[[n]])),
                                 rep("Within", length(within_individual_distances_virus_mothers[[n]])),
                                 rep("Between", length(between_individual_distances_virus_mothers[[n]])) ) ) 
  
  colnames(virusN) <- c('Distance', 'Variable')
  
  virusN$Distance <- virusN$Distance/median( unname(unlist(virus[[n]])) )
  
  Youden <- cutpointr(virusN, Distance, Variable, 
                      pos_class = "Within",
                      neg_class = "Between", 
                      method = oc_youden_kernel, 
                      metric = youden)
  cutpoints_all[n,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_all[n,"FDR_value"] <- quantile(virusN[virusN$Variable=="Between",]$Distance, probs = 0.05)
  
  cutpoints_all[n,"N_within_comparisons"] <- length(c(within_individual_distances_virus_infants[[n]], within_individual_distances_virus_mothers[[n]]))
  
  virus_hists_data[[names(virus[n])]] <- virusN
}

distance_histograms_all <- list()

for (n in 1:NROW(virus_hists_data)) {
  
  distance_histograms_all[[names(virus_hists_data[n])]] <- ggplot(virus_hists_data[[n]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(virus_hists_data[[n]]$Distance) ), alpha=0.7) + 
    geom_density(aes(y = (after_stat(count)/sum(after_stat(count)))*2000), alpha=0.2   ) + 
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=cutpoints_infants[ names(virus_hists_data[n])  ,"Youden_index"]), linetype="dashed") + 
    geom_vline(aes(xintercept=cutpoints_infants[ names(virus_hists_data[n])  ,"FDR_value"]), linetype="dashed", color="red") +
    ggtitle( gsub("\\...\\K\\d+", "", gsub(".*_length", "L"  , names(virus_hists_data[n]) ), perl=T) ) +
    theme_bw() + 
    theme(title=element_text(size=8))
  
}


combined_plot <- distance_histograms_all[[1]]

for (h in 2:NROW(virus_hists_data)) {
  combined_plot <- combined_plot + distance_histograms_all[[h]]
}

combined_plot <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "")

pdf('./04.PLOTS/Cutpoints_within_mothers_and_infants_combined_viruses.pdf', width=27/2.54, height=24/2.54)
combined_plot
dev.off()


list_transmitted$cutpoint <- NA

for (i in list_transmitted$Virus) {
  list_transmitted[list_transmitted$Virus==i,]$cutpoint <- ifelse( cutpoints_all[i,"Youden_index"] >= cutpoints_all[i,"FDR_value"],  cutpoints_all[i,"FDR_value"],  cutpoints_all[i,"Youden_index"])
}

# based on the inter-individual strain variation, how often is the virus transmitted between related and unrelated individuals?

transmission_virus <- virus[names(virus) %in% list_transmitted$Virus]
mother_infant_distances_virus <- list()
unrelated_distances_virus <- list()

list_transmitted$N_related_pairs_transmitted <- NA
list_transmitted$N_unrelated_pairs_transmitted <- NA
list_transmitted$Perc_related_pairs_transmitted <- NA
list_transmitted$Perc_unrelated_pairs_transmitted <- NA
list_transmitted$N_related_pairs <- NA
list_transmitted$N_unrelated_pairs <- NA

for (n in 1:NROW(transmission_virus)) {
  # get the virus name and distance matrix
  
  virusN <- transmission_virus[[n]]
  virusN <- virusN/median( unname(unlist(transmission_virus[[n]])) )
  virusName <- names(transmission_virus[n])
  
  virusN[virusN <= list_transmitted[list_transmitted$Virus==virusName,]$cutpoint ] <- -1
  virusN <- virusN[grep('Mother', row.names(virusN)), grep('Infant', colnames(virusN))]
  
  # two lists to store pair-wise distances for the current virus[[n]] per family
  mother_infant_distances <- list()
  mother_unrelated_distances <- list()
  
  # loop over all families
  for (i in unique( substr(row.names( virusN ), 1, 7)) ) { # getting the FAM_IDs of all virus-positive families from the distance matrix 
    
    # Check if within the selected family both mother and infant have virus-positive samples
    if (any(grep(paste0(i, '_Mother'),  row.names(virusN) ))==T & any(grep(paste0(i, '_Infant'),  colnames(virusN) ))==T) {
      
      # Parsing pair-wise distances between maternal and infant samples from the selected family
      mother_infant_distances[[i]] <- unlist( virusN[grep( paste0(i, '_Mother') , row.names(virusN) ), grep( paste0(i, '_Infant'), colnames( virusN ) ) ] )
      
      # Parsing pair_wise distances between maternal samples and samples of unrelated individuals:
      # Assigning all pair-wise distances between samples of related individuals as NA for the selected family
      
      mother_unrelated_distances[[i]] <- unlist( virusN[grep( paste0(i, '_Mother') , row.names(virusN) ), 
                                                        grep( paste0(i, '_Infant'), colnames( virusN ), invert = T ) ] )
      
    } else {
      
      mother_infant_distances[[i]] <- NA
      mother_unrelated_distances[[i]] <- NA
      
    }
    
  }
  
  # collecting the distances per virus:
  mother_infant_distances_virus[[virusName]] <- as.numeric(na.omit(unname ( unlist(mother_infant_distances) )))
  unrelated_distances_virus[[virusName]] <- as.numeric(na.omit( unname( unlist(mother_unrelated_distances) ) ) )
  
  list_transmitted[list_transmitted$Virus==virusName,]$N_related_pairs_transmitted <- sum(mother_infant_distances_virus[[virusName]]==-1)
  list_transmitted[list_transmitted$Virus==virusName,]$N_unrelated_pairs_transmitted <- sum(unrelated_distances_virus[[virusName]]==-1)
  list_transmitted[list_transmitted$Virus==virusName,]$Perc_related_pairs_transmitted <- sum(mother_infant_distances_virus[[virusName]]==-1)/length(mother_infant_distances_virus[[virusName]])
  list_transmitted[list_transmitted$Virus==virusName,]$Perc_unrelated_pairs_transmitted <- sum(unrelated_distances_virus[[virusName]]==-1)/length(unrelated_distances_virus[[virusName]])
  list_transmitted[list_transmitted$Virus==virusName,]$N_related_pairs <- length(mother_infant_distances_virus[[virusName]])
  list_transmitted[list_transmitted$Virus==virusName,]$N_unrelated_pairs <- length(unrelated_distances_virus[[virusName]])
}

row.names(list_transmitted) <- list_transmitted$Virus
list_transmitted$p_value_kinship_transmission <- NA
for (i in list_transmitted$Virus) {
  
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

for_plot <- melt(list_transmitted[,c(52,56:57,62)])

pdf('./04.PLOTS/Perc_transmitted_in_pairs_maximized_Youden.pdf', width=15.7/2.54, height=24/2.54)
ggplot(for_plot, aes(value, ContigID_easy, fill=variable)) + 
  geom_col(position = 'dodge', alpha=0.8) +
  labs (y="Viruses", x="% transmitted in mother-infant samples pairs") + 
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
  geom_text(aes(label = significance_level_transmission, x = 1.05, y = ContigID_easy), size = 3, angle=270) 
dev.off()


pdf('./04.PLOTS/Individual_strain_variation.pdf', width=6/2.54, height=6/2.54)
ggplot(list_transmitted, aes('',y=cutpoint)) + 
  geom_violin()+
  geom_boxplot(width=0.1) +
  geom_sina() +
  labs (y="Kimura distance", x="Within-individuals") + 
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold"))
dev.off()


# Bacterial hosts of viruses that were predicted to have smaller within-family distances
host_assignment <- host_assignment[host_assignment$Virus %in% list_transmitted$Virus,]
host_assignment$species <- gsub('.*s__','',host_assignment$Host.taxonomy)
host_assignment$species <- sub(' ', '_', host_assignment$species)
host_assignment$present_metaphlan <- NA

for (i in unique(host_assignment$species)) {
  if ( length( grep(i, row.names(metaphlan_species)) )!=0 ) {
    host_assignment[host_assignment$species==i,8] <- T
  }
}
host_assignment <- host_assignment[host_assignment$species!='',]
unique(host_assignment[!is.na(host_assignment$present_metaphlan),]$species)

View(list_transmitted[list_transmitted$p_value_kinship_transmission_adj < 0.05,])

dist_example <- virus[["LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609"]]



#dist_example[dist_example<=0.818604e-05] <- 0

heatmap(as.matrix(dist_example))
#L85266_cov_2453.209609_tree <- nj(as.matrix(dist_example))
tree <- nj(as.matrix(dist_example))
plot(tree, cex=0.5)

tree <- full_join(tree, metadata, by='label')

ggtree(tree, options(ignore.negative.edge = T)) + 
  geom_tiplab(size=3, color="blue")

pdf('./04.PLOTS/PhTree_L85266_cov_2453.209609.pdf', width=29.7/2.54, height=11/2.54)
ggplot(tree, aes(x, y), options(ignore.negative.edge = TRUE) ) + 
  geom_tree() + 
  geom_tiplab(size=3, color="blue") +
  theme_tree() 
dev.off()

pdf('./04.PLOTS/PhTree_L87092_cov_9507.222710.pdf', width=29.7/2.54, height=11/2.54)
ggplot(L85266_cov_2453.209609_tree, aes(x, y), options(ignore.negative.edge=TRUE)) + 
  geom_tree() + 
  geom_tiplab(size=3, color="blue") +
  theme_tree() 
dev.off()

metaphlan_species <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header = T)
microbiome_species <- metaphlan_species[grep('s__', row.names(metaphlan_species) ) ,]
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_with_phenos.txt', sep='\t', header=T)
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels=c('P3','P7', 'B', 'M1', 'M2', 'M3', 'M6', 'M9','M12'), ordered=T)
MGS_metadata$Type <- factor(MGS_metadata$Type, levels=c('Mother', 'Infant'), ordered=T)
bacteroides_uniformis <- microbiome_species[grep('Bacteroides_uniformis', row.names(microbiome_species)),grep('0154B', colnames(microbiome_species)) ]
row.names(bacteroides_uniformis) <- 'Bacteroides uniformis'
library(reshape2)

plot_host <- melt(bacteroides_uniformis)
colnames(plot_host)[1] <- 'Short_sample_ID_bact'
plot_host <- merge(plot_host, MGS_metadata, by='Short_sample_ID_bact')

pdf('./04.PLOTS/Bacteroides_uniformis_FAM0154.pdf', width=14/2.54, height=8/2.54)
ggplot(plot_host, aes(Timepoint, value, color=Individual_ID, group=Individual_ID) ) + 
  geom_point() + 
  geom_line() +
  facet_wrap(~Type) +
  theme_bw()+
  scale_y_log10() + 
  labs (y="log-transformed abundance", x="Timepoint") +
  ggtitle('Bacteroides uniformis') + 
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(colour = "Type") 
dev.off()

##############################
# OUTPUT
##############################
write.table(host_assignment, '02.CLEAN_DATA/Table_predicted_hosts_shared_strains.txt', sep='\t', quote=F, row.names = F)
write.table(unique(host_assignment[!is.na(host_assignment$present_metaphlan),]$species), '02.CLEAN_DATA/List_predicted_hosts_shared_strains_260523.txt', sep='\t', quote=F, row.names = F, col.names = F)
write.table(metadata[(metadata$Alex_ID %in% colnames(virus[['LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609']])) &
                       metadata$source=='VLP', ]$NG_ID, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/L85266_LS0_VLP_positive_list.txt', sep='\t', quote=F, row.names = F, col.names = F)
write.table(metadata[(metadata$Alex_ID %in% colnames(virus[['LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609']])) &
                       metadata$source=='MGS', ]$NG_ID, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/L85266_LS0_MGS_positive_list.txt', sep='\t', quote=F, row.names = F, col.names = F)
write.table(metadata[(metadata$Alex_ID %in% colnames(virus[['LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438']])) &
                       metadata$source=='VLP', ]$NG_ID, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/L37775_LS1_VLP_positive_list.txt', sep='\t', quote=F, row.names = F, col.names = F)
write.table(metadata[(metadata$Alex_ID %in% colnames(virus[['LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438']])) &
                       metadata$source=='MGS', ]$NG_ID, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/L37775_LS1_MGS_positive_list.txt', sep='\t', quote=F, row.names = F, col.names = F)

