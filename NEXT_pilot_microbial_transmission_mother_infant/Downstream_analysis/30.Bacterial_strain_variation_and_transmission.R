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
library(cutpointr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(ggforestplot)
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

species_names <- read.table('02.CLEAN_DATA/List_bacterial_strains_reconstructed_and_tested.txt', sep='\t', header=T)

list_transmitted <- read.table('02.CLEAN_DATA/List_bacteria_results_checking_transmission.txt', sep='\t', header=T)
list_transmitted <- list_transmitted[list_transmitted$FDR<= 0.05,]
list_transmitted <- list_transmitted[order(list_transmitted$Host_species),]

bacterium_to_test <- bacterium[names(bacterium) %in% list_transmitted$Host_SGB]

selected_viruses <- read.table('02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt', sep='\t', header=T)

##############################
# ANALYSIS
##############################

# defining the inter-personal bacterium genetic distance variation

# creating lists to keep the respective distances:
within_infant_distances_bacterium <- list()
between_infants_distances_bacterium <- list()

within_mother_distances_bacterium <- list()
between_mothers_distances_bacterium <- list()

# loop to go over infant and maternal samples separately:
for (p in c("Infant", "Mother")) {
  
  # loop to go over all bacteria:
  for (n in 1:NROW(bacterium_to_test)) {
    
    # progress notification
    print(paste0(' > prepping within-infant distances for strain ', names(bacterium_to_test[n])))  
    
    # get the bacterium name and distance matrix
    bacteriumN <- bacterium_to_test[[n]]
    bacteriumName <- names(bacterium_to_test[n])
    
    # keep only infants or mothers in the matrix
    bacteriumN <- bacteriumN[grep(p, row.names(bacteriumN)), grep(p, colnames(bacteriumN))]
    
    # create empty lists to keep the within- and between-individual distances for the current bacterium selection:
    within_individual_distances <- list()
    between_individual_distances <- list()
    
    # collect pairwise distances per infant or mother:
    for (i in unique( gsub( paste0("(FAM0[0-9]{3}_", p, "_[A-Za-z]).*"), "\\1", colnames(bacteriumN) ) ) ) {
      
      # collect distances only for infants or mothers that have at least 2 longitudinal samples:
      if (length(grep(i, colnames(bacteriumN)))>1) {
        
        # get all pairwise distances to longitudinal samples of the infant or mother:
        A  <- bacteriumN[grep(i, colnames(bacteriumN)),grep(i, colnames(bacteriumN))]
        within_individual_distances[[i]] <- A[upper.tri(A)]
        
        ## if treating twins as unrelated individuals:
        #between_individual_distances[[i]] <- unlist(bacteriumN[grep(i, colnames(bacteriumN)),grep(i, colnames(bacteriumN), invert=T)])
        
        ## if treating twins as related individuals:
        # get all pairwise distances to samples of unrelated individuals by excluding columns with the samples of longitudinal samples/samples of the twin
        between_individual_distances[[i]] <- unlist(bacteriumN[grep(i, colnames(bacteriumN)), grep(substr(i, 1, 14), colnames(bacteriumN), invert=T)])
      } 
      
    }
    
    
    if (p=='Infant') {
      
      # keep the collected infant distances for the current bacterium selection:
      within_infant_distances_bacterium[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(within_individual_distances) )))
      between_infants_distances_bacterium[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(between_individual_distances) )))
      
    } else {
      
      # keep the collected infant distances for the current bacterium selection:
      within_mother_distances_bacterium[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(within_individual_distances) )))
      between_mothers_distances_bacterium[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(between_individual_distances) )))
      
    }
    
  }
  
}

#### Threshold for strain-sharing based on individual strain variation in both mothers and infants combined: 

cutpoints_all <- data.frame( matrix(NA, nrow = length(bacterium_to_test), ncol=3 )  )
row.names(cutpoints_all) <- names(bacterium_to_test)
colnames(cutpoints_all) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

bacterium_hists_data <- list()
for (bacteriumName in names(bacterium_to_test)) {
  
  bacteriumN <- data.frame( c(within_infant_distances_bacterium[[bacteriumName]],
                              between_infants_distances_bacterium[[bacteriumName]],
                              within_mother_distances_bacterium[[bacteriumName]],
                              between_mothers_distances_bacterium[[bacteriumName]]),
                            c(rep("Within", length(within_infant_distances_bacterium[[bacteriumName]])),
                              rep("Between", length(between_infants_distances_bacterium[[bacteriumName]])),
                              rep("Within", length(within_mother_distances_bacterium[[bacteriumName]])),
                              rep("Between", length(between_mothers_distances_bacterium[[bacteriumName]])) ) ) 
  
  colnames(bacteriumN) <- c('Distance', 'Variable')
  
  bacteriumN <- bacteriumN[bacteriumN$Distance!=Inf,]
  
  A <- bacterium_to_test[[bacteriumName]]
  A <- A[upper.tri(A)]
  
  bacteriumN$Distance <- bacteriumN$Distance/median( A, na.rm = T )
  
  Youden <- cutpointr(bacteriumN, Distance, Variable, 
                      pos_class = "Within",
                      neg_class = "Between", 
                      method = oc_youden_kernel, 
                      metric = youden, 
                      na.rm=T)
  cutpoints_all[bacteriumName,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_all[bacteriumName,"FDR_value"] <- quantile(bacteriumN[bacteriumN$Variable=="Between",]$Distance, probs = 0.05)
  
  cutpoints_all[bacteriumName,"N_within_comparisons"] <- length(c(within_infant_distances_bacterium[[bacteriumName]], within_mother_distances_bacterium[[bacteriumName]]))
  
  # lines that are needed for the correct threshold depiction at the plots later
  bacteriumN$Youden <- NA
  bacteriumN[1,"Youden"] <- Youden$optimal_cutpoint
  bacteriumN$FDR_value <- NA
  bacteriumN[1,"FDR_value"] <- quantile(bacteriumN[bacteriumN$Variable=="Between",]$Distance, probs = 0.05)
  bacteriumN$N_within_comparisons <- NA
  bacteriumN[1,"N_within_comparisons"] <- length(c(within_infant_distances_bacterium[[bacteriumName]], within_mother_distances_bacterium[[bacteriumName]]))
  
  bacterium_hists_data[[bacteriumName]] <- bacteriumN
}

# for alphabetic order:
bacterium_hists_data <- bacterium_hists_data[list_transmitted$Host_SGB]

distance_histograms_all <- list()

for (bacteriumName in names(bacterium_hists_data)) {
  
  distance_histograms_all[[bacteriumName]] <- ggplot(bacterium_hists_data[[bacteriumName]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(bacterium_hists_data[[bacteriumName]]$Distance) ), alpha=0.7) + 
    geom_density(aes(x=Distance, fill=Variable), alpha=0.2) +
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=Youden[1], color="Youden_index"), linetype="dashed") + 
    geom_vline(aes(xintercept=FDR_value[1], color="FDR_value"), linetype="dashed") +
    annotate(geom = "text", label=paste0("N=",bacterium_hists_data[[bacteriumName]]$N_within_comparisons[1]), x=Inf, y=Inf, hjust=+1.1,vjust=+2, size=3) +
    ggtitle( gsub('_', ' ', species_names$species[match(bacteriumName, species_names$Host_SGB)] ) ) +
    theme_bw() + 
    theme(title = element_text(size=7), 
          axis.title = element_text(size=9),
          axis.text = element_text(size=9),
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size=10)) +
    scale_color_manual(name="Threshold", 
                       labels=c(Youden_index="Youden index", FDR_value="5% FDR"),
                       values=c(Youden_index="black", FDR_value="red"), 
                       guide = guide_legend(order = 2)) + 
    scale_fill_manual(name="Distribution", 
                      labels=c(Between="Unrelated individual comparison", Within="Same individual comparison"),
                      values=c(Between="#F8766D", Within="#00BFC4"), 
                      guide = guide_legend(order = 1))
}


combined_plot <- distance_histograms_all[[1]]

for (h in 2:NROW(bacterium_hists_data)) {
  combined_plot <- combined_plot + distance_histograms_all[[h]]
}

combined_plot <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "") & theme(legend.position = "bottom")

pdf('./04.PLOTS/Cutpoints_within_mothers_and_infants_combined_bacterium_twins_as_related_oc_youden_kernel.pdf', width=30/2.54, height=24/2.54)
combined_plot
dev.off()

### Threshold calculation in infants 

# removing those that do not have between- and/or within- individual distances
within_infant_distances_bacterium <- within_infant_distances_bacterium[ lapply(within_infant_distances_bacterium, length) > 0  ]
between_infants_distances_bacterium <- between_infants_distances_bacterium[ lapply(between_infants_distances_bacterium, length) > 0 ]

identical( names(within_infant_distances_bacterium), names(between_infants_distances_bacterium))

cutpoints_infants <- data.frame( matrix(NA, nrow = length(within_infant_distances_bacterium), ncol=3 )  )
row.names(cutpoints_infants) <- names(within_infant_distances_bacterium)
colnames(cutpoints_infants) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

infants_bacterium_hists_data <- list()

for (bacteriumName in names(within_infant_distances_bacterium)) {
  
  bacteriumN <- data.frame( c(within_infant_distances_bacterium[[bacteriumName]],
                                 between_infants_distances_bacterium[[bacteriumName]]),
                               c(rep("Within", length(within_infant_distances_bacterium[[bacteriumName]])),
                                 rep("Between", length(between_infants_distances_bacterium[[bacteriumName]])) ) ) 
  
  colnames(bacteriumN) <- c('Distance', 'Variable')
  bacteriumN <- bacteriumN[bacteriumN$Distance!=Inf,]
  
  A <- bacterium_to_test[[bacteriumName]]
  # if normalazing on infant-only pairwise distances
  # A <- A[grep('Infant', row.names(A)), grep('Infant', row.names(A))]
  A <- A[upper.tri(A)]
  
  bacteriumN$Distance <- bacteriumN$Distance/median( A, na.rm=T )
  
  if ( sum(bacteriumN[bacteriumN$Variable=='Within',]$Distance, na.rm = T)!=0 & length(bacteriumN[bacteriumN$Variable=="Within",]$Distance) > 1 ) {
    Youden <- cutpointr(bacteriumN, Distance, Variable, 
                        pos_class = "Within",
                        neg_class = "Between", 
                        method = oc_youden_kernel, 
                        metric = youden, 
                        na.rm=T)
  } else {
    Youden <- cutpointr(bacteriumN, Distance, Variable, 
                        pos_class = "Within",
                        neg_class = "Between", 
                        method = maximize_metric, 
                        metric = youden, 
                        na.rm=T)
  }
  
  cutpoints_infants[bacteriumName,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_infants[bacteriumName,"FDR_value"] <- quantile(bacteriumN[bacteriumN$Variable=="Between",]$Distance, probs = 0.05, na.rm=T)
  
  cutpoints_infants[bacteriumName,"N_within_comparisons"] <- length(within_infant_distances_bacterium[[bacteriumName]])
  
  # lines that are needed for the correct threshold depiction at the plots later
  bacteriumN$Youden <- NA
  bacteriumN[1,"Youden"] <- Youden$optimal_cutpoint
  bacteriumN$FDR_value <- NA
  bacteriumN[1,"FDR_value"] <- quantile(bacteriumN[bacteriumN$Variable=="Between",]$Distance, probs = 0.05)
  bacteriumN$N_within_comparisons <- NA
  bacteriumN[1,"N_within_comparisons"] <- length(within_infant_distances_bacterium[[bacteriumName]])
  
  infants_bacterium_hists_data[[bacteriumName]] <- bacteriumN
  
}

# for alphabetic order:
infants_bacterium_hists_data <- infants_bacterium_hists_data[list_transmitted$Host_SGB]
infants_bacterium_hists_data <- infants_bacterium_hists_data[ lapply(infants_bacterium_hists_data, length) > 0]

distance_histograms <- list()

for (bacteriumName in names(infants_bacterium_hists_data)) {
  
  distance_histograms[[bacteriumName]] <- ggplot(infants_bacterium_hists_data[[bacteriumName]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(infants_bacterium_hists_data[[bacteriumName]]$Distance) ), alpha=0.7) + 
    geom_density(aes(x=Distance, fill=Variable), alpha=0.2) + 
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=Youden[1], color="Youden_index"), linetype="dashed") + 
    geom_vline(aes(xintercept=FDR_value[1], color="FDR_value"), linetype="dashed") +
    annotate(geom = "text", label=paste0("N=",infants_bacterium_hists_data[[bacteriumName]]$N_within_comparisons[1]), x=Inf, y=Inf, hjust=+1.1,vjust=+2, size=3) +
    ggtitle( gsub('_', ' ', species_names$species[match(bacteriumName, species_names$Host_SGB)] ) ) +
    theme_bw() + 
    theme(title = element_text(size=7), 
          axis.title = element_text(size=9),
          axis.text = element_text(size=9),
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size=10)) +
    scale_color_manual(name="Threshold", 
                       labels=c(Youden_index="Youden index", FDR_value="5% FDR"),
                       values=c(Youden_index="black", FDR_value="red"), 
                       guide = guide_legend(order = 2)) + 
    scale_fill_manual(name="Distribution", 
                      labels=c(Between="Unrelated individual comparison", Within="Same individual comparison"),
                      values=c(Between="#F8766D", Within="#00BFC4"), 
                      guide = guide_legend(order = 1))
  
}


combined_plot <- distance_histograms[[1]]

for (h in 2:NROW(infants_bacterium_hists_data)) {
  combined_plot <- combined_plot + distance_histograms[[h]]
}


combined_plot <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "") & theme(legend.position = "bottom")

pdf('./04.PLOTS/Cutpoints_within_infants_bacterium_twins_as_related_oc_youden_kernel.pdf', width=29/2.54, height=24/2.54)
combined_plot
dev.off()

#### Threshold in mothers
within_mother_distances_bacterium <- within_mother_distances_bacterium[ lapply(within_mother_distances_bacterium, length) > 0  ]
between_mothers_distances_bacterium <- between_mothers_distances_bacterium[ lapply(between_mothers_distances_bacterium, length) > 0 ]

identical( names(within_mother_distances_bacterium), names(between_mothers_distances_bacterium))

cutpoints_mothers <- data.frame( matrix(NA, nrow = length(within_mother_distances_bacterium), ncol=3 )  )
row.names(cutpoints_mothers) <- names(within_mother_distances_bacterium)
colnames(cutpoints_mothers) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

mothers_bacterium_hists_data <- list()
for (bacteriumName in names(within_mother_distances_bacterium)) {
  
  bacteriumN <- data.frame( c(within_mother_distances_bacterium[[bacteriumName]],
                                 between_mothers_distances_bacterium[[bacteriumName]]),
                               c(rep("Within", length(within_mother_distances_bacterium[[bacteriumName]])),
                                 rep("Between", length(between_mothers_distances_bacterium[[bacteriumName]])) ) ) 
  
  colnames(bacteriumN) <- c('Distance', 'Variable')
  bacteriumN <- bacteriumN[bacteriumN$Distance!=Inf,]
  
  A <- bacterium_to_test[[bacteriumName]]
  # if normalazing on maternal-only pairwise distances
  # A <- A[grep('Mother', row.names(A)), grep('Mother', row.names(A))]
  A <- A[upper.tri(A)]
  
  bacteriumN$Distance <- bacteriumN$Distance/median( A, na.rm=T)
  
  Youden <- cutpointr(bacteriumN, Distance, Variable, 
                      pos_class = "Within",
                      neg_class = "Between", 
                      method = oc_youden_kernel, 
                      metric = youden, 
                      na.rm=T)
  cutpoints_mothers[bacteriumName,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_mothers[bacteriumName,"FDR_value"] <- quantile(bacteriumN[bacteriumN$Variable=="Between",]$Distance, probs = 0.05, na.rm=T)
  
  
  cutpoints_mothers[bacteriumName,"N_within_comparisons"] <- length(within_mother_distances_bacterium[[bacteriumName]])
  
  # lines that are needed for the correct threshold depiction at the plots later
  bacteriumN$Youden <- NA
  bacteriumN[1,"Youden"] <- Youden$optimal_cutpoint
  bacteriumN$FDR_value <- NA
  bacteriumN[1,"FDR_value"] <- quantile(bacteriumN[bacteriumN$Variable=="Between",]$Distance, probs = 0.05)
  bacteriumN$N_within_comparisons <- NA
  bacteriumN[1,"N_within_comparisons"] <- length(within_mother_distances_bacterium[[bacteriumName]])
  
  mothers_bacterium_hists_data[[bacteriumName]] <- bacteriumN
}


# for alphabetic order:
mothers_bacterium_hists_data <- mothers_bacterium_hists_data[list_transmitted$Host_SGB]
mothers_bacterium_hists_data <- mothers_bacterium_hists_data[ lapply(mothers_bacterium_hists_data, length) > 0]

distance_histograms_mothers <- list()

for (bacteriumName in names(mothers_bacterium_hists_data)) {
  
  distance_histograms_mothers[[bacteriumName]] <- ggplot(mothers_bacterium_hists_data[[bacteriumName]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(mothers_bacterium_hists_data[[bacteriumName]]$Distance) ), alpha=0.7) + 
    geom_density(aes(x=Distance, fill=Variable), alpha=0.2) +  
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=Youden[1], color="Youden_index"), linetype="dashed") + 
    geom_vline(aes(xintercept=FDR_value[1], color="FDR_value"), linetype="dashed") +
    annotate(geom = "text", label=paste0("N=",mothers_bacterium_hists_data[[bacteriumName]]$N_within_comparisons[1]), x=Inf, y=Inf, hjust=+1.1,vjust=+2, size=3) +
    ggtitle( gsub('_', ' ', species_names$species[match(bacteriumName, species_names$Host_SGB)] ) ) +
    theme_bw() + 
    theme(title = element_text(size=7), 
          axis.title = element_text(size=9),
          axis.text = element_text(size=9),
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size=10)) +
    scale_color_manual(name="Threshold", 
                       labels=c(Youden_index="Youden index", FDR_value="5% FDR"),
                       values=c(Youden_index="black", FDR_value="red"), 
                       guide = guide_legend(order = 2)) + 
    scale_fill_manual(name="Distribution", 
                      labels=c(Between="Unrelated individual comparison", Within="Same individual comparison"),
                      values=c(Between="#F8766D", Within="#00BFC4"), 
                      guide = guide_legend(order = 1))
  
}


combined_plot_mothers <- distance_histograms_mothers[[1]]

for (h in 2:NROW(mothers_bacterium_hists_data)) {
  combined_plot_mothers <- combined_plot_mothers + distance_histograms_mothers[[h]]
}

combined_plot_mothers <- combined_plot_mothers +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "") & theme(legend.position = "bottom")

pdf('./04.PLOTS/Cutpoints_within_mothers_bacterium_oc_youden_kernel.pdf', width=29/2.54, height=24/2.54)
combined_plot_mothers
dev.off()

##### comparison of cutpoints: cutpoints calculated for maternal samples only (max 6 months apart) and for mothers and babies together (max 12 months apart)
cutpoints_compare <- rbind(cutpoints_all, cutpoints_mothers)
cutpoints_compare$source <- c(rep("Combined", length(cutpoints_all$Youden_index)), rep("Mothers", length(cutpoints_mothers$Youden_index)))
cutpoints_compare$Host_SGB <- row.names(cutpoints_compare)
# is Youden index different? (No, p-value: 0.2238)
wilcox.test(cutpoints_compare[cutpoints_compare$Host_SGB!="SGB17247",]$Youden_index ~ cutpoints_compare[cutpoints_compare$Host_SGB!="SGB17247",]$source, paired=T)
# is 5% FDR value different? (No, p-value: 0.2873), but the N of within comparisons is different as well
wilcox.test(cutpoints_compare[cutpoints_compare$Host_SGB!="SGB17247",]$FDR_value ~ cutpoints_compare[cutpoints_compare$Host_SGB!="SGB17247",]$source, paired=T)

#### Chose to use the combined mother-infant cutpoints
list_transmitted$cutpoint <- NA

for (i in list_transmitted$Host_SGB) {
  list_transmitted[list_transmitted$Host_SGB==i,]$cutpoint <- ifelse( cutpoints_all[i,"Youden_index"] >= cutpoints_all[i,"FDR_value"],  cutpoints_all[i,"FDR_value"],  cutpoints_all[i,"Youden_index"])
}

selected_viruses$Host_Youden_index <- cutpoints_all$Youden_index[match(selected_viruses$Host_SGB, row.names(cutpoints_all))]
selected_viruses$Host_FDR_ipv_Youden <- cutpoints_all$FDR_value[match(selected_viruses$Host_SGB, row.names(cutpoints_all))]
selected_viruses$N_within_comparisons_host <- cutpoints_all$N_within_comparisons[match(selected_viruses$Host_SGB, row.names(cutpoints_all))]

# based on the inter-individual strain variation, how often is the bacterium transmitted between related and unrelated individuals?

transmission_bacterium <- bacterium_to_test[names(bacterium_to_test) %in% list_transmitted$Host_SGB]
mother_infant_distances_bacterium <- list()
unrelated_distances_bacterium <- list()

row.names(list_transmitted) <- list_transmitted$Host_SGB
list_transmitted$N_related_pairs_transmitted <- NA
list_transmitted$N_unrelated_pairs_transmitted <- NA
list_transmitted$Perc_related_pairs_transmitted <- NA
list_transmitted$Perc_unrelated_pairs_transmitted <- NA
list_transmitted$N_related_pairs <- NA
list_transmitted$N_unrelated_pairs <- NA

for (bacteriumName in names(transmission_bacterium)) {
  # get the virus name and distance matrix

  bacteriumN <- transmission_bacterium[[bacteriumName]]
  
  A <- bacteriumN[upper.tri(bacteriumN)]
  
  bacteriumN <- bacteriumN/median( unname(unlist(A)), na.rm = T )

  bacteriumN[bacteriumN <= list_transmitted[list_transmitted$Host_SGB==bacteriumName,]$cutpoint ] <- -1
  bacteriumN <- bacteriumN[grep('Mother', row.names(bacteriumN)), grep('Infant', colnames(bacteriumN))]
  
  # two lists to store pair-wise distances for the current bacterium per family
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
  
  list_transmitted[bacteriumName,]$N_related_pairs_transmitted <- sum(mother_infant_distances_bacterium[[bacteriumName]]==-1)
  list_transmitted[bacteriumName,]$N_unrelated_pairs_transmitted <- sum(unrelated_distances_bacterium[[bacteriumName]]==-1)
  list_transmitted[bacteriumName,]$Perc_related_pairs_transmitted <- sum(mother_infant_distances_bacterium[[bacteriumName]]==-1)/length(mother_infant_distances_bacterium[[bacteriumName]])
  list_transmitted[bacteriumName,]$Perc_unrelated_pairs_transmitted <- sum(unrelated_distances_bacterium[[bacteriumName]]==-1)/length(unrelated_distances_bacterium[[bacteriumName]])
  list_transmitted[bacteriumName,]$N_related_pairs <- length(mother_infant_distances_bacterium[[bacteriumName]])
  list_transmitted[bacteriumName,]$N_unrelated_pairs <- length(unrelated_distances_bacterium[[bacteriumName]])
}


selected_viruses$Host_Transmitted_in_N_related <- list_transmitted$N_related_pairs_transmitted[match(selected_viruses$Host_SGB, list_transmitted$Host_SGB)]
selected_viruses$Host_Transmitted_in_N_unrelated <- list_transmitted$N_unrelated_pairs_transmitted[match(selected_viruses$Host_SGB, list_transmitted$Host_SGB)]

list_transmitted$p_value_kinship_transmission <- NA

testing_enrichment_transmission <- list()
for (i in list_transmitted$Host_SGB) {
  
  testing_enrichment_transmission[[i]] <- matrix( c(list_transmitted[i,"N_related_pairs_transmitted"], 
                                                    (list_transmitted[i, "N_related_pairs"] - list_transmitted[i,"N_related_pairs_transmitted"]),
                                                    list_transmitted[i,"N_unrelated_pairs_transmitted"],
                                                    (list_transmitted[i, "N_unrelated_pairs"] - list_transmitted[i,"N_unrelated_pairs_transmitted"])),
                                                  nrow=2,
                                                  dimnames = list(Transmitted=c('Related', 'Unrelated'),
                                                                  Not_transmitted=c('Related', 'Unrelated'))) 
  
  check_association_with_kinship <- fisher.test(testing_enrichment_transmission[[i]], alternative="greater")
  
  list_transmitted[i,]$p_value_kinship_transmission <- check_association_with_kinship$p.value
  
}

list_transmitted$p_value_kinship_transmission_adj <- p.adjust(list_transmitted$p_value_kinship_transmission, method = "BH")

selected_viruses$Host_FDR_enirched_related_transmission <- list_transmitted$p_value_kinship_transmission_adj[match(selected_viruses$Host_SGB, list_transmitted$Host_SGB)]
selected_viruses$Host_Transmission_enriched_in_related <- ifelse( (selected_viruses$Host_FDR_enirched_related_transmission <= 0.05), "YES", "NO" )

list_transmitted$significance_level_transmission <- NA
list_transmitted[list_transmitted$p_value_kinship_transmission_adj > 0.05,]$significance_level_transmission <- 'ns'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.05,]$significance_level_transmission <- '*'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.01,]$significance_level_transmission <- '**'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.001,]$significance_level_transmission <- '***'
list_transmitted$species_names <- gsub('_', ' ',list_transmitted$Host_species)

# calculate how many unique mother-infant pairs are positive for the bacterium:
# create empty lists to keep the N unique mother-infant pairs sharing the bacterium:
related_positive_bacterium <- list()
unrelated_positive_bacterium <- list()

for (bacteriumName in names(bacterium_to_test) ) {
  
  bacteriumN <- bacterium_to_test[[bacteriumName]]
  # reformat symmetric distance matrix to only contain mother-infant distances (mothers in rows and infants in columns)
  bacteriumN <- bacteriumN[grep('Mother', row.names(bacteriumN)), grep('Infant', colnames(bacteriumN))]
  
  # two lists to store N positive pairs
  related_positive_families <- list()
  unrelated_positive_families <- list()
  
  # loop over all families
  for (i in unique( substr(row.names( bacteriumN ), 1, 7)) ) { # getting the FAM_IDs of all bacterium-positive families from the distance matrix 
    
    # Check if within the selected family both mother and infant have virus-positive samples
    if (any(grep(paste0(i, '_Mother'),  row.names(bacteriumN) ))==T & any(grep(paste0(i, '_Infant'),  colnames(bacteriumN) ))==T) {
      
      bacteriumNM <- bacteriumN[grep(paste0(i, '_Mother'),  row.names(bacteriumN) ),, drop=F]
      # Parsing pair-wise distances between maternal and infant samples from the selected family
      related_positive_families[[i]] <- length(unique( substr(grep(paste0(i, '_Infant'),  colnames(bacteriumNM), value=T ), 1, 16) ))
      
      # Parsing pair_wise distances between maternal samples and samples of unrelated individuals:
      # Assigning all pair-wise distances between samples of related individuals as NA for the selected family
      
      unrelated_positive_families[[i]] <- length(unique( substr(grep(paste0(i, '_Infant'),  colnames(bacteriumNM), value=T, invert=T ), 1, 16) ))
      
    } else {
      
      related_positive_families[[i]] <- NA
      unrelated_positive_families[[i]] <- NA
      
    }
    
  }
  
  # collecting the distances per virus:
  related_positive_bacterium[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(related_positive_families) )))
  unrelated_positive_bacterium[[bacteriumName]] <- as.numeric(na.omit( unname( unlist(unrelated_positive_families) ) ) )
}

list_transmitted$N_unique_positive_related_pairs <- NA
for (i in list_transmitted$Host_SGB) {
  
  list_transmitted[list_transmitted$Host_SGB==i,]$N_unique_positive_related_pairs <- sum(related_positive_bacterium[[i]])
}

list_transmitted$N_unique_positive_related_pairs_from_32 <- paste0(list_transmitted$N_unique_positive_related_pairs, '/32')

list_transmitted$ord <- sprintf("%02i", 26:1)

for_plot <- melt(list_transmitted[,c("Host_SGB", 
                                     "species_names", 
                                     "Perc_related_pairs_transmitted", 
                                     "Perc_unrelated_pairs_transmitted", 
                                     "significance_level_transmission", 
                                     "N_unique_positive_related_pairs_from_32",
                                     "ord")])


p1 <- ggplot(for_plot, aes(value, ord, fill=variable)) + 
  geom_col(position = 'dodge', alpha=0.8) +
  labs (y="Bacteria", x="% mother-infant sample\npairs with shared bacterium") + 
  scale_y_discrete(labels = setNames(as.character(for_plot$species_name), for_plot$ord)) +
  theme_bw()+
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=6,face="bold"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=7, face="bold"),
        legend.position = "bottom") +
  labs(fill="Kinship") + 
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  geom_text(aes(label = significance_level_transmission, x = 1.05, y = ord), size = 1.6, angle=270) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE, title.position = 'left'))

 
p2 <- ggplot(list_transmitted, aes(N_unique_positive_related_pairs_from_32, ord)) + 
  geom_stripes(odd = "#33333333", even = "#00000000") + 
  labs (y="Viruses", x="N positive\nmother-infant\ndyads") +
  scale_y_discrete(labels = setNames(as.character(for_plot$species_name), for_plot$ord)) +
  geom_text(aes(label = N_unique_positive_related_pairs_from_32, x=1.5, y = ord, fontface = "bold"), color="#17B971",size=2    ) + 
  theme_bw() + 
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=6, face="bold"),
        legend.position = "none") +
  scale_x_discrete(limits=c("0/32", "0/32")) 

transmission_enrichment <- p1+p2 + 
  plot_layout(ncol = 2, nrow = 1, guides="collect", widths = c(5, 1.2)) + 
  plot_annotation(title = "") & theme(legend.position = 'bottom', 
                                      legend.title = element_text(size=7, face="bold"),
                                      legend.text = element_text(size=6)) 

pdf('./04.PLOTS/Perc_transmitted_in_pairs_maximized_Youden_bacteria_with_N.pdf', width=9/2.54, height=12/2.54)
transmission_enrichment
dev.off()

pdf('./04.PLOTS/Perc_transmittedin_pairs_maximized_Youden_bacteria.pdf', width=7/2.54, height=12/2.54)
p1
dev.off()

##############################
# OUTPUT
##############################
write.table(selected_viruses, "02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt", sep='\t', row.names = F, quote=F)
write.table(list_transmitted[,-c(1:6)], '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Perc_transmitted_bacteria_in_pairs_maximized_Youden_with_N.txt', sep='\t', quote=F)
