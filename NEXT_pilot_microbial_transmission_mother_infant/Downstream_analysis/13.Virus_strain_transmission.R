setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we compare genetic pairwise distances of virus
# consensuses reconstructed from VLP and MGS samples between
# related and unrelated mother-infant pairs
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

#library(scales)
##############################
# Input data
##############################

selected_viruses <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/Viruses_shared_min_5_families_UPD_final.txt', sep='\t', header=T)
# remove columns that are not informative for the further analysis:
exlude_columns <- c("excluding", "final_exclusion", "Order$", "Percent_of_votes", "Family$", "Percent_of_votes.1", "geNomad_taxonomy$")
selected_viruses <- selected_viruses[,grep(paste(exlude_columns,collapse="|"), 
                                          colnames(selected_viruses), invert = T)]
colnames(selected_viruses)[1] <- 'Virus'

folder = "02.CLEAN_DATA/dist.matrices/"  
file_list = list.files(path=folder, pattern="*.txt")  
virus=lapply(paste0(folder, file_list), function(x) read.table(x, sep='\t', header=T))
names(virus) <- gsub(".dist.txt", '', file_list)

##############################
# ANALYSIS
##############################

# two lists to store pair-wise distances per virus:
mother_infant_distances_virus <- list()
unrelated_distances_virus <- list()

# loop through each virus from the selection
for (n in 1:NROW(virus) ) {

  print(paste0(' > prepping within-mother-infant distances for strain ', names(virus[n])))  
  
  # get the virus name and distance matrix
  virusN <- virus[[n]]
  virusName <- names(virus[n])
  # reformat symmetric distance matrix to only contain mother-infant distances (mothers in rows and infants in columns)
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
                                                        grep( paste0(i, '_Infant'), colnames(virusN), invert = T ) ] )
      
    } else {
      
      mother_infant_distances[[i]] <- NA
      mother_unrelated_distances[[i]] <- NA
      
    }

  }
  
  # collecting the distances per virus:
  mother_infant_distances_virus[[virusName]] <- as.numeric(na.omit(unname ( unlist(mother_infant_distances) )))
  unrelated_distances_virus[[virusName]] <- as.numeric(na.omit( unname( unlist(mother_unrelated_distances) ) ) )
}



# testing if distances in mother-infant pairs are smaller than between unrelated individuals

p_value_real <- as.data.frame(matrix(NA, nrow=length(virus), ncol=3))
colnames(p_value_real)[c(1:3)] <- c("Virus", "N_related_distances", "p_value")

plot_distances <- data.frame()

for (n in 1:NROW(virus)) {
  p_value_real[n,"Virus"] <- names(virus[n])
  p_value_real[n,"N_related_distances"] <- length(mother_infant_distances_virus[[n]])
  
  if (length(unrelated_distances_virus[[n]])!=0) {
    
    vector4analysis = c(mother_infant_distances_virus[[n]], unrelated_distances_virus[[n]])
    factor4analysis = c(rep("Related",length(mother_infant_distances_virus[[n]])),
                        rep("Unrelated",length(unrelated_distances_virus[[n]])))
    
    wilcoxon_real = wilcox.test(vector4analysis ~ factor4analysis, alternative='less', paired=F)
    p_value_real[n,"p_value"] <- wilcoxon_real$p.value
    
    #table for plot
    virus_name <- c( rep( names(virus[n]), length(factor4analysis) ) )
    plot_distances <- rbind(plot_distances, data.frame(virus_name, factor4analysis, vector4analysis))
    
  } else {
    # if there were no distances to unrelated infants for the mother:
    p_value_real[n,3] <- "Comparison not possible"
  
    }
  
}

selected_viruses <- merge(selected_viruses, p_value_real[,c("Virus", "N_related_distances","p_value")], by='Virus', all.x = T)
colnames(selected_viruses)[length(selected_viruses)] <- 'Related_unrelated_comparison_possible'
selected_viruses$Related_unrelated_comparison_possible <- ifelse(selected_viruses$Related_unrelated_comparison_possible=='Comparison not possible',
                                                                 "NO",
                                                                 "YES")

#permutations test

# storing F-statistics for permuted tables
p_value_perm <- list()

# loop over all viruses
for (n in 1:NROW(virus) ) {
  
  # creating a vector for storing F-statistics for the virus[[n]]
  p_value <- as.numeric(c())
  
  # loop over 1000 permutations
  for (i in 1:1000) {
    
    virusN <- virus[[n]]
    virusName <- names(virus[n])
    virusN <- virusN[grep('Mother', row.names(virusN)), grep('Infant', colnames(virusN))]
    
    # randomly permuting the samples:
    # first get the random permutation of columns
    sample_permut = sample( 1:ncol(virusN) )
    # second get the random permutation of columns
    sample_permut_row = sample( 1:nrow(virusN) )
    
    # reorder the table of distances for the chosen virus according to the permutation
    #perm_table = virusN[sample_permut,sample_permut]
    perm_table = virusN[sample_permut_row,sample_permut]
    # rename columns
    colnames(perm_table) = colnames(virusN)
    row.names(perm_table) = row.names(virusN)
    
    mother_infant_distances <- list()
    mother_unrelated_distances <- list()
    
    # creating a vector for storage of pair-wise distances between related samples for this iteration of permutation
    mother_infant_distances_virus_perm <- c()
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
    
    mother_infant_distances_virus_perm <- as.numeric(na.omit(unname ( unlist(mother_infant_distances) )))
    mother_unrelated_distances_perm  <- as.numeric(na.omit( unname( unlist(mother_unrelated_distances) ) ) )

    #calculating p-value for permuted table
    
    if (length(mother_unrelated_distances_perm)!=0) {
      
      vector4analysis <- c(mother_infant_distances_virus_perm, mother_unrelated_distances_perm )
      factor4analysis <- c( rep('Related', length(mother_infant_distances_virus_perm) ),
                            rep('Unrelated', length(mother_unrelated_distances_perm ) ) )
      wilcox_perm = wilcox.test(vector4analysis ~ factor4analysis, alternative='less', paired=F)
      
      # storing p-vaule for this permutation
      p_value[i] <- wilcox_perm$p.value
      
    } else {
      p_value[i] <- "Comparison not possible"
    }
    
  }
  
  # storing p-value of permutations for every virus
  p_value_perm[names(virus[n])] <- as.data.frame(p_value)
}

# calculation of p-value based on permutations:
p_value_real$p_value_adj <- NA

# loop goes over all viruses
for (i in p_value_real$Virus) {
  
  # calculating the probability of event if the effect is random
  p_value_real[p_value_real$Virus==i,"p_value_adj"] <- sum(p_value_perm[[i]] <= as.numeric(p_value_real[p_value_real$Virus==i,3]))/1000
}

# calculating FDR
family_viruses <- p_value_real[!is.na(p_value_real$p_value_adj),]
family_viruses$FDR <- p.adjust(family_viruses$p_value_adj, method = "BH")

selected_viruses <- merge(selected_viruses, family_viruses[,c("Virus","FDR")], by='Virus', all.x = T)
colnames(selected_viruses)[length(colnames(selected_viruses))] <- "FDR_dist_comparison"
selected_viruses$Distances_Related_lower <- ifelse(selected_viruses$FDR_dist_comparison<=0.05, "YES", "NO")

if (identical(selected_viruses$Virus, names(unrelated_distances_virus))) {
  
  selected_viruses$N_unrelated_distances <- as.numeric(lapply(unrelated_distances_virus, function(x) {
    length(x)
  }))
  
}


##### PLOTS ####
# adding the smallest non-zero Kimura distance to all distances (to use logarithmic scale in the plot)
plot_distances$vector4analysis <- plot_distances$vector4analysis + min(plot_distances[plot_distances$vector4analysis!=0,]$vector4analysis)
# showing only those viruses that have more than 5 pair-wise distances for related samples
plot_distances_select <- plot_distances[ plot_distances$virus_name %in% names(virus)[lengths(mother_infant_distances_virus)>5] ,]

# renaming contigs for easier perception:
plot_distances_select$easy_name <- selected_viruses$ContigID_easy[match(plot_distances_select$virus_name, selected_viruses$Virus)]

# color-coding the y-axis titles depending on statistical significance of the differnece:
myPalette <- family_viruses
myPalette$color <- NA
myPalette[myPalette$FDR>0.05,]$color <- 'grey'
myPalette[myPalette$FDR<=0.05,]$color <- 'black'
myPalette <- myPalette[myPalette$Virus %in% unique(plot_distances_select$virus_name),]
myPalette$easy_name <-  selected_viruses$ContigID_easy[match(myPalette$Virus, selected_viruses$Virus)]
myPalette <- myPalette[order(myPalette$easy_name),]

# all virus strains
p1 <- ggplot(plot_distances_select, aes(vector4analysis,easy_name, fill=factor4analysis)) + 
  labs (y="Viruses", x="Log-scaled\nKimura distance") + 
  geom_sina(aes(fill=factor4analysis), size=1.3, alpha=0.8, shape=21, stroke=0) +
  geom_boxplot(outlier.shape = NA, alpha=0.3) +
  theme_bw()+
  scale_x_log10() +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=8,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  labs(fill="Kinship", color="") + geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) + 
  theme(axis.text.y = element_text(colour=myPalette$color),
        legend.position = "none")


N_pairwise_distance <- melt(selected_viruses[,c("Virus", "ContigID_easy", "N_related_distances", "N_unrelated_distances")])
N_pairwise_distance$value <- N_pairwise_distance$value+1


p2 <- ggplot(N_pairwise_distance, aes(value, ContigID_easy, fill=variable) ) + 
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

pdf('./04.PLOTS/Infant_virus_strain_all_more_than_7_fams_wilcox_less_final.pdf', width=9/2.54, height=17/2.54)
combined_plot
dev.off()

# those where distances are significantly different:

myPalette$N_positive_families <- selected_viruses$N_positive_families[match(myPalette$Virus, selected_viruses$Virus)]
myPalette <- myPalette[myPalette$FDR<0.05 & myPalette$N_positive_families >=7,]
plot_distances_select <- plot_distances_select[plot_distances_select$virus_name %in% myPalette$Virus,]

N_pairwise_distance <- N_pairwise_distance[N_pairwise_distance$Virus %in% myPalette$Virus,]

p3 <- ggplot(plot_distances_select, aes(vector4analysis, easy_name, fill=factor4analysis)) + 
  labs (y="Viruses", x="Log-scaled\nKimura distance") + 
  geom_sina(aes(fill=factor4analysis), size=0.9, alpha=0.8, shape=21, stroke=0) +
  geom_boxplot(outlier.shape = NA, alpha=0.3) +
  theme_bw()+
  scale_x_log10() +
  theme(axis.text=element_text(size=6), 
        axis.title=element_text(size=6,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8, face="bold"),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  labs(fill="Kinship", color="") + geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) + 
  theme(axis.text.y = element_text(colour=myPalette$color),
        legend.position = "none")

p4 <- ggplot(N_pairwise_distance, aes(value, ContigID_easy, fill=variable) ) + 
  geom_bar(stat='identity', position='dodge', color="black", alpha=0.5) +
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

pdf('./04.PLOTS/Infant_virus_strains_significant_more_than_7_fams_wilcox_less_final.pdf', width=7.7/2.54, height=14/2.54)
combined_plot_significant
dev.off()
##############################
# OUTPUT
##############################
write.table(plot_distances, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Virus_distances_related_vs_unrelated.txt', sep='\t', quote=F, row.names=F)
write.table(selected_viruses, '02.CLEAN_DATA/List_viruses_selected_transmission_metadata.txt', sep='\t', quote=F, row.names=F)
write.table(family_viruses, '02.CLEAN_DATA/List_viruses_results_checking_transmission.txt', sep='\t', quote=F, row.names = F)

