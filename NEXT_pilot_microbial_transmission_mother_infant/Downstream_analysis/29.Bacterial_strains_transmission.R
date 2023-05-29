setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore inter-individual variation of bacterial 
# strains that are predicted to be infected by the viruses of
# interest
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

library(scales)

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
colnames(species_names) <- c('SGB', 'Taxonomy')
species_names$SGB <- names(bacterium)

for (i in species_names$SGB) {
  species_names[species_names$SGB==i,"Taxonomy"] <- metaphlan$clade_name[grep(i, metaphlan$clade_name)]
}

species_names$Species <- gsub('.*s__','', sapply(strsplit(species_names$Taxonomy, '\\|'), "[", 7))

species_names$Present_in_mothers_and_infants <- ifelse(species_names$SGB %in% names(bacterium_to_test), "YES", "NO")

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
colnames(p_value_real)[c(1:3)] <- c("Bacterium", "N_related_distances", "p_value")

plot_distances <- data.frame()

for (n in 1:NROW(bacterium_to_test)) {
  p_value_real[n,1] <- names(bacterium_to_test[n])
  p_value_real[n,2] <- length(mother_infant_distances_bacterium[[n]])
  if (length(unrelated_distances_bacterium[[n]])!=0) {
    vector4analysis = c(mother_infant_distances_bacterium[[n]], unrelated_distances_bacterium[[n]])
    factor4analysis = c(rep("Related",length(mother_infant_distances_bacterium[[n]])),
                        rep("Unrelated",length(unrelated_distances_bacterium[[n]])))
    wilcoxon_real = wilcox.test(vector4analysis ~ factor4analysis, alternative='less', paired=F)
    p_value_real[n,3] <- wilcoxon_real$p.value
    
    #table for plot
    bacterium_name <- c( rep( names(bacterium[n]), length(factor4analysis) ) )
    plot_distances <- rbind(plot_distances, data.frame(bacterium_name, factor4analysis, vector4analysis))
  } else {
    p_value_real[n,3] <- "Comparison not possible"
  }
  
}

p_value_real <- p_value_real[p_value_real$p_value!='Comparison not possible',]
bacterium_to_test <- bacterium_to_test[names(bacterium_to_test) %in% p_value_real$Bacterium]

species_names$Related_unrelated_compared <- ifelse(species_names$SGB %in% p_value_real$Bacterium, "YES", "NO")
species_names$Comments_comparison <- NA
species_names[species_names$SGB=="SGB15317",]$Comments_comparison <- "No related pairs"
species_names[species_names$SGB=="SGB1846",]$Comments_comparison <- "No related pairs"
species_names[species_names$SGB=="SGB1853",]$Comments_comparison <- "No unrelated pairs"

# storing F-statistics for permuted tables
p_value_perm <- list()

# loop over all bacteriumes
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
      
      # storing p-vaule for this permutation
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
for (i in p_value_real$Bacterium) {
  
  # calculating the probability of event if the effect is random
  p_value_real[p_value_real$Bacterium==i,"p_value_adj"] <- sum(p_value_perm[[i]] <= as.numeric(p_value_real[p_value_real$Bacterium==i,3]))/1000
}

# calculating FDR
family_bacteria <- p_value_real[!is.na(p_value_real$p_value_adj),]
family_bacteria$FDR <- p.adjust(family_bacteria$p_value_adj, method = "BH")


##### PLOTS ####
# adding the smallest non-zero Kimura distance to all distances (to use logarithmic scale in the plot)
plot_distances$vector4analysis <- plot_distances$vector4analysis + min(plot_distances[plot_distances$vector4analysis!=0,]$vector4analysis)
# showing only those viruses that have more than 5 pair-wise distances for related samples
plot_distances_select <- plot_distances[ plot_distances$bacterium_name %in% family_bacteria[family_bacteria$N_related_distances>5,]$Bacterium,]


# color-coding the y-axis titles depending on statistical significance of the differnece:
myPalette <- family_bacteria
myPalette$color <- NA
myPalette[myPalette$FDR>0.05,]$color <- 'grey'
myPalette[myPalette$FDR<=0.05,]$color <- 'black'
myPalette <- myPalette[myPalette$Bacterium %in% unique(plot_distances_select$bacterium_name),]

# renaming contigs for easier perception:
plot_distances_select$easy_name <- species_names$Species[match(plot_distances_select$bacterium_name, species_names$SGB)]
plot_distances_select$easy_name <- paste0(plot_distances_select$bacterium_name, '_', plot_distances_select$easy_name)


# all virus strains
pdf('./04.PLOTS/Infant_bacteria_strain_all_wilcoxon_less.pdf', width=29.7/2.54, height=21/2.54)
ggplot(plot_distances_select, aes(vector4analysis,easy_name, fill=factor4analysis)) + 
  geom_boxplot(outlier.shape = NA,alpha=0.5) +
  labs (y="Bacterial strains", x="Log-scaled Kimura distance") + 
  geom_sina(aes(fill=factor4analysis), size=0.6,alpha=0.5) +
  #geom_jitter(,size=0.6,alpha=0.5) +
  theme_bw()+
  scale_x_log10() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(fill="Kinship") + geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) + 
  theme(axis.text.y = element_text(colour=myPalette$color) )
dev.off()


# those where distances are significantly different:

myPalette <- myPalette[myPalette$FDR<0.05,]
plot_distances_select <- plot_distances_select[plot_distances_select$bacterium_name %in% myPalette$Bacterium,]


pdf('./04.PLOTS/Infant_bacteria_strains_significant_wilcoxon_less.pdf', width=29.7/2.54, height=21/2.54)
ggplot(plot_distances_select, aes(vector4analysis, easy_name, fill=factor4analysis)) + 
  geom_boxplot(outlier.shape = NA,alpha=0.5) +
  labs (y="Bacterial strains", x="Log-scaled Kimura distance") + 
  geom_sina(aes(fill=factor4analysis), size=0.6,alpha=0.5) +
  #geom_jitter(aes(),size=0.6,alpha=0.5) +
  theme_bw()+
  scale_x_log10() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  labs(fill="Kinship") + geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) + 
  theme(axis.text.y = element_text(colour=myPalette$color) )
dev.off()

family_bacteria$Species <- species_names$Species[match(family_bacteria$Bacterium, species_names$SGB)]
##############################
# OUTPUT
##############################

write.table(family_bacteria, '02.CLEAN_DATA/List_bacteria_results_checking_transmission.txt', sep='\t', quote=F, row.names = F)
write.table(species_names, '02.CLEAN_DATA/List_bacterial_strains_reconstructed_and_tested.txt', sep='\t', quote=F, row.names=F)

