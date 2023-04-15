setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the bacterial genera prevalence in 
# infant and adult gut
#############################################################

##############################
# Functions
##############################
taxonomy_abundance <- function(taxonomy_table) {
  ##Function to calculate mean excluding 0 values
  nzmean <- function(a){
    mean(a[a!=0])
  }
  ##Function to calculate nº of 0
  zsum <- function(a){
    sum (a==0)
  }
  ##Function to calculate nº of non-0
  nsum <- function(a){
    sum (a!=0)
  }
  ## Loop for each column (taxonomy) in the taxonomy table
  for (i in 1:ncol(taxonomy_table)) {
    #Calculate mean for each column
    aa = mean(taxonomy_table[,i])
    #Calculate number of non-zeros (individuals)
    bb = nsum(taxonomy_table[,i])
    #Calculate mean without taking into account the 0
    cc = nzmean(taxonomy_table[,i])
    #Calculate number of zeros 
    dd = zsum(taxonomy_table[,i])
    ee= (dd/(dd+bb))*100
    my_results[i,1] = aa
    my_results[i,2] = bb
    my_results[i,3] = cc
    my_results[i,4] = dd
    my_results[i,5] = ee
  }
  return(my_results)
}
##############################
# Loading libraries
##############################
library(ggplot2)
library(stringr)
library(data.table)
##############################
# Input data
##############################

metaphlan_genera <- read.table('02.CLEAN_DATA/Microbiome_genera_unfiltred.txt', sep='\t', header=T)
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_with_phenos.txt', sep='\t', header=T)

##############################
# ANALYSIS
##############################
my_results=matrix(ncol = 5, nrow=ncol( as.data.frame( t(metaphlan_genera[,colnames(metaphlan_genera) %in% MGS_metadata[MGS_metadata$Type=='Infant',]$Short_sample_ID_bact,]) ) ) )
my_results=as.data.frame( taxonomy_abundance(as.data.frame( t(metaphlan_genera[,colnames(metaphlan_genera) %in% MGS_metadata[MGS_metadata$Type=='Infant',]$Short_sample_ID_bact,]) ) ) )
rownames(my_results) = row.names(metaphlan_genera)
colnames(my_results) = c("Mean","Prevalence", "nz_mean", "N_of_0", "perc_missing") 

my_results_infants <- my_results

my_results=matrix(ncol = 5, nrow=ncol( as.data.frame( t(metaphlan_genera[,colnames(metaphlan_genera) %in% MGS_metadata[MGS_metadata$Type=='Mother',]$Short_sample_ID_bact,]) ) ) )
my_results=as.data.frame( taxonomy_abundance(as.data.frame( t(metaphlan_genera[,colnames(metaphlan_genera) %in% MGS_metadata[MGS_metadata$Type=='Mother',]$Short_sample_ID_bact,]) ) ) )
rownames(my_results) = row.names(metaphlan_genera)
colnames(my_results) = c("Mean","Prevalence", "nz_mean", "N_of_0", "perc_missing") 

my_results_mothers <- my_results

my_results_mothers_top20 <- my_results_mothers[order(my_results_mothers$Prevalence, decreasing = F),]
my_results_mothers_top20 <- my_results_mothers_top20[c(584:603),]
my_results_mothers_top20$Type <- 'Mother'
my_results_mothers_top20$taxa <- row.names(my_results_mothers_top20)

my_results_infants_top20 <- my_results_infants[order(my_results_infants$Prevalence, decreasing = F),]
my_results_infants_top20 <- my_results_infants_top20[c(584:603),]
my_results_infants_top20$Type <- 'Infant'
my_results_infants_top20$taxa <- row.names(my_results_infants_top20)


my_results <- rbind(my_results_mothers_top20,my_results_infants_top20)
my_results$Type <- factor(my_results$Type, levels=c('Mother','Infant'), ordered=T)
my_results$genus <- gsub('.*g__', '', my_results$taxa)
my_results$phylum <- gsub('.*p__','', sapply(strsplit(my_results$taxa, '\\|'), "[", 2))

my_results <- as.data.table(my_results)
my_results[, ord := sprintf("%02i", frank(my_results, Type, Prevalence, ties.method = "first"))]


pdf('./04.PLOTS/Bacterial_genera_prevalence.pdf', width=15.7/2.54, height=21/2.54)
ggplot(my_results, aes(x=Prevalence, y=ord, fill=phylum)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Type, ncol=1, scales = 'free_y', drop=T) +
  labs (y="Bacterial genus", x="Detected in\nN samples") +
  scale_y_discrete(labels = setNames(as.character(my_results$genus), my_results$ord)) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_blank()) +
  scale_fill_manual(values=c("#D69A44", "#44729D","#9E3729","#9B58A2"))
dev.off()


##############################
# OUTPUT
##############################
