setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore inter-individual variation of virus 
# 'strains'/consensuses reconstructed from both VLP and MGS 
# samples in order to estimate the number of strain sharing
# events.
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


selected_viruses <- read.table('02.CLEAN_DATA/List_viruses_selected_transmission_metadata.txt', sep='\t', header=T)

list_transmitted <- read.table('02.CLEAN_DATA/List_viruses_results_checking_transmission.txt', sep='\t', header=T)
list_transmitted <- list_transmitted[list_transmitted$FDR<= 0.05,]

folder = "02.CLEAN_DATA/dist.matrices/"  
file_list = list.files(path=folder, pattern="*.txt")  
virus=lapply(paste0(folder, file_list), function(x) read.table(x, sep='\t', header=T))
names(virus) <- gsub(".dist.txt", '', file_list)

virus <- virus[names(virus) %in% list_transmitted$Virus]

host_assignment <- read.csv('01.RAW_DATA/Host_assignment/MERGED_Host_prediction_to_genome_m90.csv')

metaphlan <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)

##############################
# ANALYSIS
##############################

# collecting the within- and between-individual pairwise distances for infants and mothers separately

# creating lists to keep the respective distances:
within_infant_distances_virus <- list()
between_infants_distances_virus <- list()

within_mother_distances_virus <- list()
between_mothers_distances_virus <- list()


# loop to go over infant and maternal samples separately:
for (p in c("Infant", "Mother")) {
  
  # loop to go over all viruses:
  for (n in 1:NROW(virus)) {
    
    # progress notification
    print(paste0(' > prepping within-infant distances for strain ', names(virus[n])))  
    
    # get the virus name and distance matrix
    virusName <- names(virus[n])
    virusN <- virus[[n]]
    
    # keep only infants or mothers in the matrix
    virusN <- virusN[grep(p, row.names(virusN)), grep(p, colnames(virusN))]
    
    # create empty lists to keep the within- and between-individual distances for the current virus selection:
    within_individual_distances <- list()
    between_individual_distances <- list()
    
    # collect pairwise distances per infant or mother:
    for (i in unique( gsub( paste0("(FAM0[0-9]{3}_", p, "_[A-Za-z]).*"), "\\1", colnames(virusN) ) ) ) {
      
      # collect distances only for infants or mothers that have at least 2 longitudinal samples:
      if (length(grep(i, colnames(virusN)))>1) {
        
        # get all pairwise distances to longitudinal samples of the infant or mother:
        A  <- virusN[grep(i, colnames(virusN)),grep(i, colnames(virusN))]
        within_individual_distances[[i]] <- A[upper.tri(A)]
        
        ## if treating twins as unrelated individuals:
        #between_individual_distances[[i]] <- unlist(virusN_infants[grep(i, colnames(virusN_infants)),grep(i, colnames(virusN_infants), invert=T)])
        
        ## if treating twins as related individuals:
        # get all pairwise distances to samples of unrelated individuals by excluding columns with the samples of longitudinal samples/samples of the twin
        between_individual_distances[[i]] <- unlist(virusN[grep(i, colnames(virusN)), grep(substr(i, 1, 14), colnames(virusN), invert=T)])
        
      } 
      
    }
    
    # keep the collected infant distances for the current virus selection:
    
    if (p=='Infant') {
      
      within_infant_distances_virus[[virusName]] <- as.numeric(na.omit(unname ( unlist(within_individual_distances) )))
      between_infants_distances_virus[[virusName]] <- as.numeric(na.omit(unname ( unlist(between_individual_distances) )))
      
    } else {
      
      # keep the collected infant distances for the current virus selection:
      within_mother_distances_virus[[virusName]] <- as.numeric(na.omit(unname ( unlist(within_individual_distances) )))
      between_mothers_distances_virus[[virusName]] <- as.numeric(na.omit(unname ( unlist(between_individual_distances) )))
      
    }
    
  }
  
}

#### Threshold merged individuals (infants and mothers together): 

cutpoints_all <- data.frame( matrix(NA, nrow = length(virus), ncol=3 )  )
row.names(cutpoints_all) <- names(virus)
colnames(cutpoints_all) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

virus_hists_data <- list()
for (n in 1:NROW(virus)) {
  
  virusN <- data.frame( c(within_infant_distances_virus[[n]],
                          between_infants_distances_virus[[n]],
                          within_mother_distances_virus[[n]],
                          between_mothers_distances_virus[[n]]),
                        c(rep("Within", length(within_infant_distances_virus[[n]])),
                          rep("Between", length(between_infants_distances_virus[[n]])),
                          rep("Within", length(within_mother_distances_virus[[n]])),
                          rep("Between", length(between_mothers_distances_virus[[n]])) ) ) 
  
  colnames(virusN) <- c('Distance', 'Variable')
  
  virusN$Distance <- virusN$Distance/median( unname(unlist(virus[[n]])) )
  
  Youden <- cutpointr(virusN, Distance, Variable, 
                      pos_class = "Within",
                      neg_class = "Between", 
                      method = oc_youden_kernel, 
                      metric = youden)
  cutpoints_all[n,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_all[n,"FDR_value"] <- quantile(virusN[virusN$Variable=="Between",]$Distance, probs = 0.05)
  
  cutpoints_all[n,"N_within_comparisons"] <- length(c(within_infant_distances_virus[[n]], within_mother_distances_virus[[n]]))
  
  virus_hists_data[[names(virus[n])]] <- virusN
}

distance_histograms_all <- list()

for (n in 1:NROW(virus_hists_data)) {
  
  distance_histograms_all[[names(virus_hists_data[n])]] <- ggplot(virus_hists_data[[n]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(virus_hists_data[[n]]$Distance) ), alpha=0.7) + 
    geom_density(aes(y = (after_stat(count)/sum(after_stat(count)))*2000), alpha=0.2   ) + 
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=cutpoints_all[ names(virus_hists_data[n])  ,"Youden_index"]), linetype="dashed") + 
    geom_vline(aes(xintercept=cutpoints_all[ names(virus_hists_data[n])  ,"FDR_value"]), linetype="dashed", color="red") +
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

pdf('./04.PLOTS/Cutpoints_within_mothers_and_infants_combined_viruses_twins_as_related_oc_youden_kernel.pdf', width=27/2.54, height=24/2.54)
combined_plot
dev.off()


### Threshold calculation in infants:

cutpoints_infants <- data.frame( matrix(NA, nrow = length(virus), ncol=3 )  )
row.names(cutpoints_infants) <- names(virus)
colnames(cutpoints_infants) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')

infants_virus_hists_data <- list()
for (n in 1:NROW(virus)) {
  
  infants_virus <- data.frame( c(within_infant_distances_virus[[n]],
                                 between_infants_distances_virus[[n]]),
                               c(rep("Within", length(within_infant_distances_virus[[n]])),
                                 rep("Between", length(between_infants_distances_virus[[n]])) ) ) 
  
  colnames(infants_virus) <- c('Distance', 'Variable')
  
  infants_virus$Distance <- infants_virus$Distance/median( unname(unlist(virus[[n]])) )
  
  # cannot use oc_youden_kernel if the positive class has only zeros
  # cannot use oc_youden_kernel if there is only 1 observation in the class
  if ( sum(infants_virus[infants_virus$Variable=='Within',]$Distance, na.rm = T)!=0 & length(infants_virus[infants_virus$Variable=="Within",]$Distance) > 1) {
    Youden <- cutpointr(infants_virus, Distance, Variable, 
                        pos_class = "Within",
                        neg_class = "Between", 
                        method = oc_youden_kernel, 
                        metric = youden)
  } else {
    Youden <- cutpointr(infants_virus, Distance, Variable, 
                        pos_class = "Within",
                        neg_class = "Between", 
                        method = maximize_metric, 
                        metric = youden)
  }
  
  
  cutpoints_infants[n,"Youden_index"] <- Youden$optimal_cutpoint
  
  cutpoints_infants[n,"FDR_value"] <- quantile(infants_virus[infants_virus$Variable=="Between",]$Distance, probs = 0.05)
  
  cutpoints_infants[n,"N_within_comparisons"] <- length(within_infant_distances_virus[[n]])
  
  infants_virus_hists_data[[names(virus[n])]] <- infants_virus
}

distance_histograms_infants <- list()

for (n in 1:NROW(infants_virus_hists_data)) {
  
  distance_histograms_infants[[names(infants_virus_hists_data[n])]] <- ggplot(infants_virus_hists_data[[n]], aes(x=Distance, fill=Variable)) + 
                                  geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(infants_virus_hists_data[[n]]$Distance) ), alpha=0.7) + 
                                  geom_density(aes(y = (after_stat(count)/sum(after_stat(count)))*2000), alpha=0.2   ) + 
                                  labs(x="Normalized Distance", y="proportion (%)") +
                                  geom_vline(aes(xintercept=cutpoints_infants[ names(infants_virus_hists_data[n])  ,"Youden_index"]), linetype="dashed") + 
                                  geom_vline(aes(xintercept=cutpoints_infants[ names(infants_virus_hists_data[n])  ,"FDR_value"]), linetype="dashed", color="red") +
                                  ggtitle( gsub("\\...\\K\\d+", "", gsub(".*_length", "L"  , names(infants_virus_hists_data[n]) ), perl=T) ) +
                                  theme_bw() + 
                                  theme(title = element_text(size=8))
  
}


combined_plot_infants <- distance_histograms_infants[[1]]

for (h in 2:NROW(infants_virus_hists_data)) {
  combined_plot_infants <- combined_plot_infants + distance_histograms_infants[[h]]
}


combined_plot_infants <- combined_plot_infants +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "")
pdf('./04.PLOTS/Cutpoints_within_infants_viruses_twins_as_related_oc_youden_kernel.pdf', width=27/2.54, height=24/2.54)
combined_plot_infants
dev.off()


#### Threshold in mothers

# removing those that do not have between-individual distances
between_mothers_distances_virus <- between_mothers_distances_virus[ lapply(between_mothers_distances_virus, length) > 0 ]
within_mother_distances_virus <- within_mother_distances_virus[ lapply(within_mother_distances_virus, length) > 0 ]

if( identical(names(between_mothers_distances_virus), names(between_mothers_distances_virus)) ) {
  
  cutpoints_mothers <- data.frame( matrix(NA, nrow = length(within_mother_distances_virus), ncol=3 )  )
  row.names(cutpoints_mothers) <- names(within_mother_distances_virus)
  colnames(cutpoints_mothers) <- c('Youden_index', 'FDR_value', 'N_within_comparisons')
  
  mothers_virus_hists_data <- list()
  for (n in 1:NROW(within_mother_distances_virus)) {
    
    mothers_virus <- data.frame( c(within_mother_distances_virus[[n]],
                                   between_mothers_distances_virus[[n]]),
                                 c(rep("Within", length(within_mother_distances_virus[[n]])),
                                   rep("Between", length(between_mothers_distances_virus[[n]])) ) ) 
    
    colnames(mothers_virus) <- c('Distance', 'Variable')
    
    
    if (median( unname(unlist(within_mother_distances_virus[[n]])) ) == 0) {
      
      mothers_virus$Distance <- mothers_virus$Distance/mean( unname(unlist(within_mother_distances_virus[[n]])) )
    
      } else {
        
      mothers_virus$Distance <- mothers_virus$Distance/median( unname(unlist(within_mother_distances_virus[[n]])) )
    
      }
    
    
    if ( sum(mothers_virus[mothers_virus$Variable=='Within',]$Distance, na.rm = T)!=0 & length(mothers_virus[mothers_virus$Variable=="Within",]$Distance) > 1  ) {
      
      Youden <- cutpointr(mothers_virus, Distance, Variable, 
                          pos_class = "Within",
                          neg_class = "Between", 
                          method = oc_youden_kernel, 
                          metric = youden)
      
    } else {
      
      Youden <- cutpointr(mothers_virus, Distance, Variable, 
                          pos_class = "Within",
                          neg_class = "Between", 
                          method = maximize_metric, 
                          metric = youden)
    
      }

    cutpoints_mothers[n,"Youden_index"] <- Youden$optimal_cutpoint
    
    
    cutpoints_mothers[n,"FDR_value"] <- quantile(mothers_virus[mothers_virus$Variable=="Between",]$Distance, probs = 0.05)
    
    cutpoints_mothers[n,"N_within_comparisons"] <- length(within_mother_distances_virus[[n]])
    
    mothers_virus_hists_data[[names(within_mother_distances_virus[n])]] <- mothers_virus
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
  
  
} else {
  
  print("Check the between and within distances")
  
}

pdf('./04.PLOTS/Cutpoints_within_mothers_viruses_oc_youden_kernel.pdf', width=27/2.54, height=24/2.54)
combined_plot_mothers
dev.off()

#### Chose to use the combined mother-infant cutpoints
list_transmitted$cutpoint <- NA

for (i in list_transmitted$Virus) {
  list_transmitted[list_transmitted$Virus==i,]$cutpoint <- ifelse( cutpoints_all[i,"Youden_index"] >= cutpoints_all[i,"FDR_value"],  cutpoints_all[i,"FDR_value"],  cutpoints_all[i,"Youden_index"])
}

selected_viruses$Youden_index <- cutpoints_all$Youden_index[match(selected_viruses$Virus, row.names(cutpoints_all))]
selected_viruses$FDR_ipv_Youden <- cutpoints_all$FDR_value[match(selected_viruses$Virus, row.names(cutpoints_all))]
selected_viruses$N_within_comparisons <- cutpoints_all$N_within_comparisons[match(selected_viruses$Virus, row.names(cutpoints_all))]

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

selected_viruses$Transmitted_in_N_related <- list_transmitted$N_related_pairs_transmitted[match(selected_viruses$Virus, list_transmitted$Virus)]
selected_viruses$Transmitted_in_N_unrelated <- list_transmitted$N_unrelated_pairs_transmitted[match(selected_viruses$Virus, list_transmitted$Virus)]

row.names(list_transmitted) <- list_transmitted$Virus
list_transmitted$p_value_kinship_transmission <- NA

testing_enrichment_transmission <- list()
for (i in list_transmitted$Virus) {
  
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

selected_viruses$FDR_enirched_related_transmission <- list_transmitted$p_value_kinship_transmission_adj[match(selected_viruses$Virus, list_transmitted$Virus)]
selected_viruses$Transmission_enriched_in_related <- ifelse( (selected_viruses$FDR_enirched_related_transmission <= 0.05), "YES", "NO" )


list_transmitted$significance_level_transmission <- NA
list_transmitted[list_transmitted$p_value_kinship_transmission_adj > 0.05,]$significance_level_transmission <- 'ns'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.05,]$significance_level_transmission <- '*'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.01,]$significance_level_transmission <- '**'
list_transmitted[list_transmitted$p_value_kinship_transmission_adj <= 0.001,]$significance_level_transmission <- '***'

list_transmitted$ContigID_easy <- selected_viruses$ContigID_easy[match(list_transmitted$Virus, selected_viruses$Virus)]


# calculate how many unique mother-infant pairs are positive for the virus:
for (n in 1:NROW(virus) ) {
  
  # get the virus name and distance matrix
  virusN <- virus[[n]]
  virusName <- names(virus[n])
  # reformat symmetric distance matrix to only contain mother-infant distances (mothers in rows and infants in columns)
  virusN <- virusN[grep('Mother', row.names(virusN)), grep('Infant', colnames(virusN))]
  
  # two lists to store N positive pairs
  related_positive_families <- list()
  unrelated_positive_families <- list()
  
  # loop over all families
  for (i in unique( substr(row.names( virusN ), 1, 7)) ) { # getting the FAM_IDs of all virus-positive families from the distance matrix 
    
    # Check if within the selected family both mother and infant have virus-positive samples
    if (any(grep(paste0(i, '_Mother'),  row.names(virusN) ))==T & any(grep(paste0(i, '_Infant'),  colnames(virusN) ))==T) {
      
      virusNM <- virusN[grep(paste0(i, '_Mother'),  row.names(virusN) ),, drop=F]
      # Parsing pair-wise distances between maternal and infant samples from the selected family
      related_positive_families[[i]] <- length(unique( substr(grep(paste0(i, '_Infant'),  colnames(virusNM), value=T ), 1, 16) ))
      
      # Parsing pair_wise distances between maternal samples and samples of unrelated individuals:
      # Assigning all pair-wise distances between samples of related individuals as NA for the selected family
      
      unrelated_positive_families[[i]] <- length(unique( substr(grep(paste0(i, '_Infant'),  colnames(virusNM), value=T, invert=T ), 1, 16) ))
      
    } else {
      
      related_positive_families[[i]] <- NA
      unrelated_positive_families[[i]] <- NA
      
    }
    
  }
  
  # collecting the distances per virus:
  related_positive_virus[[virusName]] <- as.numeric(na.omit(unname ( unlist(related_positive_families) )))
  unrelated_positive_virus[[virusName]] <- as.numeric(na.omit( unname( unlist(unrelated_positive_families) ) ) )
}

list_transmitted$N_unique_positive_related_pairs <- NA
for (i in list_transmitted$Virus) {
  
  list_transmitted[list_transmitted$Virus==i,]$N_unique_positive_related_pairs <- sum(related_positive_virus[[i]])
}

list_transmitted$N_unique_positive_related_pairs_from_32 <- paste0(list_transmitted$N_unique_positive_related_pairs, '/32')

for_plot <- melt(list_transmitted[,c("ContigID_easy", "Perc_related_pairs_transmitted", "Perc_unrelated_pairs_transmitted","significance_level_transmission", "N_unique_positive_related_pairs_from_32")])

#### Plot for proportion of related-unrelated pairs sharing the virus
p1 <- ggplot(for_plot, aes(value, ContigID_easy, fill=variable)) + 
  geom_col(position = 'dodge', alpha=1) +
  labs (y="Viruses", x="% mother-infant sample\npairs with shared virus") + 
  theme_bw()+
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=6,face="bold"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=7, face="bold"), 
        legend.position = 'bottom') +
  labs(fill="Kinship") + 
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  geom_text(data=for_plot[for_plot$variable=='Perc_related_pairs_transmitted',], aes(label = significance_level_transmission, x = 1.05, y = ContigID_easy), size = 2, angle=270)


related_positive_virus <- list()
unrelated_positive_virus <- list()

p2 <- ggplot(list_transmitted, aes(N_unique_positive_related_pairs_from_32, ContigID_easy)) + 
  geom_stripes(odd = "#33333333", even = "#00000000") + 
  labs (y="Viruses", x="N positive\nmother-infant\ndyads") +
  geom_text(aes(label = N_unique_positive_related_pairs_from_32, x=1.5, y = ContigID_easy, fontface = "bold"), color="#17B971",size=3    ) + 
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
  plot_annotation(title = "") & theme(legend.position = 'bottom') 

pdf('./04.PLOTS/Perc_transmitted_in_pairs_maximized_Youden_with_N.pdf', width=8/2.54, height=14/2.54)
transmission_enrichment
dev.off()

pdf('./04.PLOTS/Perc_transmitted_viruses_in_pairs_maximized_Youden.pdf', width=7/2.54, height=14/2.54)
p1
dev.off()

# Bacterial hosts of viruses that were predicted to have smaller within-family distances
host_assignment <- host_assignment[host_assignment$Virus %in% list_transmitted$Virus,]
host_assignment$genus <- gsub('.*g__','', sapply(strsplit(host_assignment$Host.taxonomy, '\\;'), "[", 6))
host_assignment$species <- gsub('.*s__','',host_assignment$Host.taxonomy)
host_assignment$species <- sub(' ', '_', host_assignment$species)
host_assignment$species <- gsub('_([[:upper:]])', '',host_assignment$species, perl=T)

host_assignment$present_metaphlan <- NA

for (i in unique(host_assignment$species)) {
  if ( length( grep(i, row.names(metaphlan)) )!=0 ) {
    host_assignment[host_assignment$species==i,"present_metaphlan"] <- T
  }
}
host_assignment <- host_assignment[host_assignment$species!='',]

unique(host_assignment[!is.na(host_assignment$present_metaphlan),]$species)

selected_viruses <- merge(selected_viruses, host_assignment, by='Virus', all.x=T)
colnames(selected_viruses)[length(selected_viruses)] <- "Host_in_metaphlan"

##############################
# OUTPUT
##############################
write.table(list_transmitted[,-c(1:6)], '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Perc_transmitted_viruses_in_pairs_maximized_Youden_with_N.txt', sep='\t', quote=F)

write.table(host_assignment, '02.CLEAN_DATA/Table_predicted_hosts_shared_strains_maximized_Youden_combined_wilcox_less_050623.txt', sep='\t', quote=F, row.names = F)
write.table(unique(host_assignment[!is.na(host_assignment$present_metaphlan),]$species), '02.CLEAN_DATA/List_predicted_hosts_shared_strains_050623.txt', sep='\t', quote=F, row.names = F, col.names = F)

write.table(selected_viruses, '02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt', sep='\t', quote = F, row.names=F)
#write.table(metadata[(metadata$Alex_ID %in% colnames(virus[['LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609']])) &
#                       metadata$source=='VLP', ]$NG_ID, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/L85266_LS0_VLP_positive_list.txt', sep='\t', quote=F, row.names = F, col.names = F)
#write.table(metadata[(metadata$Alex_ID %in% colnames(virus[['LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609']])) &
#                       metadata$source=='MGS', ]$NG_ID, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/L85266_LS0_MGS_positive_list.txt', sep='\t', quote=F, row.names = F, col.names = F)
#write.table(metadata[(metadata$Alex_ID %in% colnames(virus[['LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438']])) &
#                       metadata$source=='VLP', ]$NG_ID, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/L37775_LS1_VLP_positive_list.txt', sep='\t', quote=F, row.names = F, col.names = F)
#write.table(metadata[(metadata$Alex_ID %in% colnames(virus[['LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438']])) &
#                       metadata$source=='MGS', ]$NG_ID, '03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/L37775_LS1_MGS_positive_list.txt', sep='\t', quote=F, row.names = F, col.names = F)

