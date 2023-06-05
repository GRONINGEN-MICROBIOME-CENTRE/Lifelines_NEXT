setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we compare genetic pairwise distances of bacterial 
# strains from MGS samples between related and unrelated 
# mother-infant pairs
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(ggplot2)
library(ggforce)
library(ggforestplot)
library(reshape2)
library(patchwork)
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

## since my there are strains that are not present in Infants 
bacterium_to_test <- list()
for (n in 1:NROW(bacterium)) {
  
  bacteriumN <- bacterium[[n]]
  bacteriumName <- names(bacterium[n])
  
  if ( length(grep('Infant', colnames(bacteriumN)))!=0 & length(grep('Mother', colnames(bacteriumN)))!=0) {
    bacterium_to_test[[bacteriumName]] <- bacterium[[bacteriumName]]
  }
  
}

metaphlan <- read.table('01.RAW_DATA/Metaphlan4_all_samples/LLNEXT_metaphlan_4_complete_10_02_2023.txt', sep = '\t', header = T)

species_names <- data.frame(matrix(NA, nrow=length(bacterium), ncol=2))
colnames(species_names) <- c('Host_SGB', 'Host_Taxonomy')
species_names$Host_SGB <- names(bacterium)

for (i in species_names$Host_SGB) {
  species_names[species_names$Host_SGB==i,"Host_Taxonomy"] <- metaphlan$clade_name[grep(i, metaphlan$clade_name)]
}

species_names$species <- gsub('.*s__','', sapply(strsplit(species_names$Host_Taxonomy, '\\|'), "[", 7))

species_names$Present_in_mothers_and_infants <- ifelse(species_names$Host_SGB %in% names(bacterium_to_test), "YES", "NO")

selected_viruses <- read.table('02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt', sep='\t', header = T)
selected_viruses$Host_in_strainphlan <- ifelse( selected_viruses$species %in% species_names$species, 'YES', 'NO'  )
selected_viruses[is.na(selected_viruses$Host_in_metaphlan),]$Host_in_strainphlan <- NA

selected_viruses <- merge(selected_viruses, species_names, by='species', all.x=T)
colnames(selected_viruses)[c(1, 65)] <- c('Host_species', 'Host_strain_present_in_mothers_and_infants')
##############################
# ANALYSIS
##############################

# two lists to store pair-wise distances per bacterium:
mother_infant_distances_bacterium <- list()
unrelated_distances_bacterium <- list()

# loop through each bacterium from the selection
for (n in 1:NROW(bacterium_to_test) ) {
  
  print(paste0(' > prepping within-mother-infant distances for strain ', names(bacterium[n])))  
  
  # get the bacterium name and distance matrix
  bacteriumN <- bacterium_to_test[[n]]
  bacteriumName <- names(bacterium_to_test[n])
  # reformat symmetric distance matrix to only contain mother-infant distances (mothers in rows and infants in columns)
  bacteriumN <- bacteriumN[grep('Mother', row.names(bacteriumN)), grep('Infant', colnames(bacteriumN))]
  
  # two lists to store pair-wise distances for the current bacterium[[n]] per family
  mother_infant_distances <- list()
  mother_unrelated_distances <- list()
  
  # loop over all families
  for (i in unique( substr(row.names( bacteriumN ), 1, 7)) ) { # getting the FAM_IDs of all bacterium-positive families from the distance matrix 
    
    # Check if within the selected family both mother and infant have bacterium-positive samples
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
  
  # collecting the distances per bacterium:
  mother_infant_distances_bacterium[[bacteriumName]] <- as.numeric(na.omit(unname ( unlist(mother_infant_distances) )))
  unrelated_distances_bacterium[[bacteriumName]] <- as.numeric(na.omit( unname( unlist(mother_unrelated_distances) ) ) )
}


# testing if distances in mother-infant pairs are smaller than between unrelated individuals

p_value_real <- as.data.frame(matrix(NA, nrow=length(bacterium_to_test), ncol=3))
colnames(p_value_real)[c(1:3)] <- c("Host_SGB", "N_related_distances_host", "p_value")

plot_distances <- data.frame()

for (n in 1:NROW(bacterium_to_test)) {
  p_value_real[n,"Host_SGB"] <- names(bacterium_to_test[n])
  p_value_real[n,"N_related_distances_host"] <- length(mother_infant_distances_bacterium[[n]])
  
  if (length(unrelated_distances_bacterium[[n]])!=0) {
    
    vector4analysis = c(mother_infant_distances_bacterium[[n]], unrelated_distances_bacterium[[n]])
    factor4analysis = c(rep("Related",length(mother_infant_distances_bacterium[[n]])),
                        rep("Unrelated",length(unrelated_distances_bacterium[[n]])))
    
    wilcoxon_real = wilcox.test(vector4analysis ~ factor4analysis, alternative='less', paired=F)
    p_value_real[n,"p_value"] <- wilcoxon_real$p.value
    
    #table for plot
    bacterium_name <- c( rep( names(bacterium[n]), length(factor4analysis) ) )
    plot_distances <- rbind(plot_distances, data.frame(bacterium_name, factor4analysis, vector4analysis))
    
  } else {
    
    p_value_real[n,"p_value"] <- "Comparison not possible"
    
  }
  
}


selected_viruses <- merge(selected_viruses, p_value_real[,c("Host_SGB", "N_related_distances_host","p_value")], by='Host_SGB', all.x = T)
colnames(selected_viruses)[length(selected_viruses)] <- 'Host_Related_unrelated_comparison_possible'
selected_viruses$Host_Related_unrelated_comparison_possible <- ifelse(selected_viruses$Host_Related_unrelated_comparison_possible=='Comparison not possible',
                                                                 "NO",
                                                                 "YES")

p_value_real <- p_value_real[p_value_real$p_value!='Comparison not possible',]
bacterium_to_test <- bacterium_to_test[names(bacterium_to_test) %in% p_value_real$Host_SGB]

species_names$Related_unrelated_compared <- ifelse(species_names$Host_SGB %in% p_value_real$Host_SGB, "YES", "NO")
species_names$Comments_comparison <- NA
species_names[species_names$Host_SGB=="SGB15317",]$Comments_comparison <- "No related pairs"
species_names[species_names$Host_SGB=="SGB1846",]$Comments_comparison <- "No related pairs"
species_names[species_names$Host_SGB=="SGB1853",]$Comments_comparison <- "No unrelated pairs"

# storing F-statistics for permuted tables
p_value_perm <- list()

# loop over all bacteria
for (n in 1:NROW(bacterium_to_test) ) {
  
  # creating a vector for storing F-statistics for the bacterium[[n]]
  p_value <- as.numeric(c())
  
  # loop over 1000 permutations
  for (i in 1:1000) {
    
    bacteriumN <- bacterium_to_test[[n]]
    bacteriumName <- names(bacterium_to_test[n])
    bacteriumN <- bacteriumN[grep('Mother', row.names(bacteriumN)), grep('Infant', colnames(bacteriumN)), drop=F]
    
    # randomly permuting the samples:
    # first get the random permutation of columns
    sample_permut = sample( 1:ncol(bacteriumN) )
    # second get the random permutation of columns
    sample_permut_row = sample( 1:nrow(bacteriumN) )
    
    # reorder the table of distances for the chosen bacterium according to the permutation
    #perm_table = bacteriumN[sample_permut,sample_permut]
    perm_table = as.data.frame(bacteriumN[sample_permut_row,sample_permut])
    # rename columns
    colnames(perm_table) = colnames(bacteriumN)
    row.names(perm_table) = row.names(bacteriumN)
    
    mother_infant_distances <- list()
    mother_unrelated_distances <- list()
    
    # creating a vector for storage of pair-wise distances between related samples for this iteration of permutation
    mother_infant_distances_bacterium_perm <- c()
    # creating a vector for storage of pair-wise distances between unrelated samples for this iteration of permutation
    mother_unrelated_distances_perm <- c()
    
    for (j in unique( substr(row.names( perm_table ), 1, 7)) ) {
      
      if (any(grep(paste0(j, '_Mother'),  row.names(perm_table) ))==T & any(grep(paste0(j, '_Infant'),  colnames(perm_table) ))==T) {
        
        mother_infant_distances[[j]] <- unlist( perm_table[grep( paste0(j, '_Mother') , row.names(perm_table) ), grep( paste0(j, '_Infant'), colnames( perm_table ) ) ] )
        
        mother_unrelated_distances[[j]] <- unlist( perm_table[grep( paste0(j, '_Mother') , row.names(perm_table) ), 
                                                              grep( paste0(j, '_Infant'), colnames( perm_table ), invert = T ) ] )
        
      } else {
        mother_infant_distances[[j]] <- NA
        mother_unrelated_distances[[j]] <- NA
      }
      
    }
    
    mother_infant_distances_bacterium_perm <- as.numeric(na.omit(unname ( unlist(mother_infant_distances) )))
    mother_unrelated_distances_perm  <- as.numeric(na.omit( unname( unlist(mother_unrelated_distances) ) ) )
    
    #calculating p-value for permuted table
    
    if (length(mother_unrelated_distances_perm)!=0) {
      
      vector4analysis <- c(mother_infant_distances_bacterium_perm, mother_unrelated_distances_perm )
      factor4analysis <- c( rep('Related', length(mother_infant_distances_bacterium_perm) ),
                            rep('Unrelated', length(mother_unrelated_distances_perm ) ) )
      wilcox_perm = wilcox.test(vector4analysis ~ factor4analysis, alternative='less', paired=F)
      
      # storing p-value for this permutation
      p_value[i] <- wilcox_perm$p.value
      
    } else {
      p_value[i] <- "Comparison not possible"
    }
    
  }
  
  # storing p-value of permutations for every bacterium
  p_value_perm[names(bacterium_to_test[n])] <- as.data.frame(p_value)
}


# calculation of p-value based on permutations:

p_value_real$p_value_adj <- NA

# loop goes over all viruses
for (i in p_value_real$Host_SGB) {
  
  # calculating the probability of event if the effect is random
  p_value_real[p_value_real$Host_SGB==i,"p_value_adj"] <- sum(p_value_perm[[i]] <= as.numeric(p_value_real[p_value_real$Host_SGB==i,'p_value']))/1000
}

# calculating FDR of pairwise distances comparison between samples of related and unrelated individuals
p_value_real$FDR <- p.adjust(p_value_real$p_value_adj, method = "BH")

family_bacteria <- p_value_real

selected_viruses <- merge(selected_viruses, p_value_real[,c("Host_SGB","FDR")], by='Host_SGB', all.x = T)
colnames(selected_viruses)[length(colnames(selected_viruses))] <- "Host_FDR_dist_comparison"
selected_viruses$Host_Distances_Related_lower <- ifelse(selected_viruses$Host_FDR_dist_comparison<=0.05, "YES", "NO")

selected_viruses$N_unrelated_distances_host <- NA

for (h in unique(selected_viruses[!is.na(selected_viruses$Host_SGB),]$Host_SGB) ) {
  
  selected_viruses[selected_viruses$Host_SGB==h & !is.na(selected_viruses$Host_SGB),]$N_unrelated_distances_host <- length(unrelated_distances_bacterium[[h]])
  
}

selected_viruses$Host_easy_name <- paste0(selected_viruses$Host_SGB, '_', selected_viruses$Host_species)
selected_viruses[is.na(selected_viruses$Host_SGB),]$Host_easy_name <- NA
##### PLOTS ####
# adding the smallest non-zero Kimura distance to all distances (to use logarithmic scale in the plot)
plot_distances$vector4analysis <- plot_distances$vector4analysis + min(plot_distances[plot_distances$vector4analysis!=0,]$vector4analysis)
# showing only those viruses that have more than 5 pair-wise distances for related samples
plot_distances_select <- plot_distances[ plot_distances$bacterium_name %in% family_bacteria[family_bacteria$N_related_distances>5,]$Host_SGB,]


# color-coding the y-axis titles depending on statistical significance of the differnece:
myPalette <- family_bacteria
myPalette$color <- NA
myPalette[myPalette$FDR>0.05,]$color <- 'grey'
myPalette[myPalette$FDR<=0.05,]$color <- 'black'
myPalette <- myPalette[myPalette$Host_SGB %in% unique(plot_distances_select$bacterium_name),]

# renaming contigs for easier perception:
plot_distances_select$easy_name <- species_names$species[match(plot_distances_select$bacterium_name, species_names$Host_SGB)]
plot_distances_select$easy_name <- paste0(plot_distances_select$bacterium_name, '_', plot_distances_select$easy_name)


# all virus strains
p1 <- ggplot(plot_distances_select, aes(vector4analysis, easy_name, fill=factor4analysis)) + 
  labs (y="Bacterial strains", x="Log-scaled Kimura distance") + 
  geom_sina(aes(fill=factor4analysis), size=1.3, alpha=0.8, shape=21, stroke=0) +
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  theme_bw()+
  scale_x_log10() +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=8,face="bold"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(fill="Kinship", color='') + 
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) + 
  theme(axis.text.y = element_text(colour=myPalette$color), 
        legend.position = 'none')

N_pairwise_distance <- melt(selected_viruses[,c("Host_SGB", "Host_easy_name", "N_related_distances_host", "N_unrelated_distances_host")])
N_pairwise_distance <- N_pairwise_distance[!is.na(N_pairwise_distance$Host_SGB),]
N_pairwise_distance <- N_pairwise_distance[N_pairwise_distance$Host_SGB %in% plot_distances_select$bacterium_name,]

p2 <- ggplot(N_pairwise_distance, aes(value, Host_easy_name, fill=variable) ) + 
  geom_bar(stat='identity', position='dodge', color="black", alpha=0.5) +
  scale_x_log10(breaks=c(10^0, 10^1, 10^3)) + 
  labs (y="", x="Log-scaled\nN distances") +
  theme_bw()+
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=8,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  labs(fill="Kinship")+
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  guides(fill = "none")

combined_plot <- p1 + p2 +
  plot_layout(ncol = 2, nrow = 1, guides="collect", widths = c(4, 1.5)) + 
  plot_annotation(title = "") & theme(legend.position = 'bottom') 


pdf('./04.PLOTS/Infant_bacterium_strain_all_wilcox_less_final.pdf', width=15/2.54, height=17/2.54)
combined_plot
dev.off()

# those where distances are significantly different:



myPalette <- myPalette[myPalette$FDR<0.05 & myPalette$N_related_distances_host >=7,]

plot_distances_select <- plot_distances_select[plot_distances_select$bacterium_name %in% myPalette$Host_SGB,]
N_pairwise_distance <- N_pairwise_distance[N_pairwise_distance$Host_SGB %in% myPalette$Host_SGB,]
N_pairwise_distance$species_names <- gsub('_',' ',species_names$species[match(N_pairwise_distance$Host_SGB, species_names$Host_SGB)])

species_order <- data.frame(sort(unique((N_pairwise_distance$species_name))))
colnames(species_order) <- "Host_species"
species_order$ord <- sprintf("%02i", 19:1)

N_pairwise_distance$ord <- species_order$ord[match(N_pairwise_distance$species_names, species_order$Host_species)]

plot_distances_select$species_names <- gsub('_', ' ', species_names$species[match(plot_distances_select$bacterium_name, species_names$Host_SGB)])
plot_distances_select$ord <- species_order$ord[match(plot_distances_select$species_names, species_order$Host_species)]


p3 <- ggplot(plot_distances_select, aes(vector4analysis, ord, fill=factor4analysis)) + 
  labs (y="Bacterial strains", x="Log-scaled Kimura distance") + 
  geom_sina(aes(fill=factor4analysis), size=0.8, alpha=0.8, shape=21, stroke=0) +
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  theme_bw()+
  scale_y_discrete(labels = setNames(as.character(plot_distances_select$species_name), plot_distances_select$ord)) +
  scale_x_log10(breaks=c(1e-01, 1e+00, 1e+01, 1e+02), labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=6,face="bold"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8, face="bold"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  labs(fill="Kinship", color="") + 
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) + 
  theme(legend.position = "none")

p4 <- ggplot(N_pairwise_distance, aes(value, ord, fill=variable) ) + 
  geom_bar(stat='identity', position='dodge', color="black", alpha=0.5) +
  scale_y_discrete(labels = setNames(as.character(N_pairwise_distance$species_name), N_pairwise_distance$ord)) +
  scale_x_log10(breaks=c(10^0, 10^1, 10^3)) + 
  labs (y="", x="Log-scaled\nN distances") +
  theme_bw()+
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=6,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8, face="bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  labs(fill="Kinship")+
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  guides(fill = "none")

combined_plot_significant <- p3 + p4 +
  plot_layout(ncol = 2, nrow = 1, guides="collect", widths = c(5, 1.5)) + 
  plot_annotation(title = "") & theme(legend.position = 'bottom') 

pdf('./04.PLOTS/Infant_bacterium_strains_significant_wilcox_less_final.pdf', width=9/2.54, height=11/2.54)
combined_plot_significant
dev.off()

family_bacteria$Host_species <- species_names$species[match(family_bacteria$Host_SGB, species_names$Host_SGB)]
##############################
# OUTPUT
##############################
write.table(plot_distances, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Bacterium_distances_related_vs_unrelated.txt', sep='\t', quote=F, row.names=F)
write.table(family_bacteria, '02.CLEAN_DATA/List_bacteria_results_checking_transmission.txt', sep='\t', quote=F, row.names = F)
write.table(species_names, '02.CLEAN_DATA/List_bacterial_strains_reconstructed_and_tested.txt', sep='\t', quote=F, row.names=F)
write.table(selected_viruses, '02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt', sep='\t', quote=F, row.names = F)
