setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we merge virus data to taxonomic levels of their 
# hosts
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(vegan)

library(dplyr)
library(lme4)
library(lmerTest)
library(tibble)

library(reshape2)
library(ggplot2)
library(ggforce)
library(ggforestplot)
##############################
# Input data
##############################
host_assignment <- read.csv('01.RAW_DATA/Host_assignment/MERGED_Host_prediction_to_genus_m90.csv', header=T)
host_assignment$Kingdom <- gsub('.*d__','', sapply(strsplit(host_assignment$Host.genus, '\\;'), "[", 1))
host_assignment$Phylum <- gsub('.*p__','', sapply(strsplit(host_assignment$Host.genus, '\\;'), "[", 2))

microbiome_phyla <- read.table('02.CLEAN_DATA/Microbiome_phyla_unfiltred.txt', sep='\t', header=T)
microbiome_phyla$Kingdom <- gsub('.*k__','', sapply(strsplit(row.names(microbiome_phyla), '\\|'), "[", 1))
microbiome_phyla$Phylum <- gsub('.*p__','', sapply(strsplit(row.names(microbiome_phyla), '\\|'), "[", 2))

microbiome_phyla_filt_CLR <- read.table('02.CLEAN_DATA/Microbiome_phyla_filtred_CLR_trans.txt', sep='\t', header=T)
colnames(microbiome_phyla_filt_CLR) <- gsub('.*p__', '',colnames(microbiome_phyla_filt_CLR))
microbiome_phyla_filt_CLR <- microbiome_phyla_filt_CLR[row.names(microbiome_phyla_filt_CLR) %in% paste0(VLP_metadata$Universal_fecal_ID, 'B'),]

## to resolve bacterial taxonomy reassignment:
microbiome <- read.table('01.RAW_DATA/Metaphlan4_all_samples/LLNEXT_metaphlan_4_complete_10_02_2023.txt', sep='\t', header=T)
tmp_taxa <- microbiome$clade_name
microbiome_phyla$Class <- gsub('.*p__','', sapply(strsplit(row.names(microbiome_phyla), '\\|'), "[", 2))
tmp_taxa <- as.data.frame(tmp_taxa)
tmp_taxa$Phylum <- gsub('.*p__','', sapply(strsplit(tmp_taxa$tmp_taxa, '\\|'), "[", 2))

unique(tmp_taxa[!tmp_taxa$Phylum %in% microbiome_phyla$Phylum,]$Phylum )

RPKM_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_with_phenos.txt', sep='\t', header=T)
row.names(VLP_metadata) <- VLP_metadata$Short_sample_ID
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_with_phenos.txt', sep='\t', header=T)
##############################
# ANALYSIS
##############################
# modifying the Phylum names in host_assignment:
host_assignment$Phylum_new <- host_assignment$Phylum

#### some manual work:
unique(host_assignment[ ! host_assignment$Phylum_new %in% microbiome_phyla$Phylum,]$Phylum)

host_assignment[grep('Firmicutes', host_assignment$Phylum),]$Phylum_new <- 'Firmicutes'
host_assignment[host_assignment$Phylum=='Bacteroidota',]$Phylum_new <- 'Bacteroidetes'
host_assignment[host_assignment$Phylum=='Verrucomicrobiota',]$Phylum_new <- 'Verrucomicrobia'
host_assignment[host_assignment$Phylum=='Actinobacteriota',]$Phylum_new <- 'Actinobacteria'
host_assignment[host_assignment$Phylum=='Synergistota',]$Phylum_new <- 'Synergistetes'
host_assignment[host_assignment$Phylum=='Fusobacteriota',]$Phylum_new <- 'Fusobacteria'
host_assignment[host_assignment$Phylum=='Spirochaetota',]$Phylum_new <- 'Spirochaetes'
host_assignment[host_assignment$Phylum=='Planctomycetota',]$Phylum_new <- 'Planctomycetes'
host_assignment[host_assignment$Phylum=='Elusimicrobiota',]$Phylum_new <- 'Elusimicrobia'
host_assignment[host_assignment$Phylum=='Thermoplasmatota',]$Phylum_new <- 'Candidatus_Thermoplasmatota'
# in metaphlan, Euryarchaeota is a phylum:
host_assignment[host_assignment$Phylum=='Methanobacteriota',]$Phylum_new <- 'Euryarchaeota'
# in metaphlan, Campylobacterota is Proteobacteria:
host_assignment[host_assignment$Phylum=='Campylobacterota',]$Phylum_new <- 'Proteobacteria'

# since Cyanobacteria and Deinococcus_Thermus are present in the metaphlan database, but not in our samples:
host_assignment <- host_assignment[-grep('Cyanobacteria', host_assignment$Phylum),]
host_assignment <- host_assignment[-grep('Deinococcota', host_assignment$Phylum),]
# single predictions that are not in metaphlan:
host_assignment <- host_assignment[-grep('Thermotogota', host_assignment$Phylum),]
host_assignment <- host_assignment[-grep('Dependentiae', host_assignment$Phylum),]
host_assignment <- host_assignment[-grep('Aquificota', host_assignment$Phylum),]
host_assignment <- host_assignment[-grep('Myxococcota', host_assignment$Phylum),]
# grouping the rest:
#host_assignment[host_assignment$Phylum %in% unique(host_assignment[ ! host_assignment$Phylum_new %in% microbiome_phyla$Phylum,]$Phylum),]$Phylum_new <- 'Other'

host_assignment_phylum <- host_assignment
host_assignment_phylum <- host_assignment [!duplicated(host_assignment[c('Virus','Phylum_new')]),]
host_assignment_phylum <- host_assignment_phylum[order(host_assignment_phylum$Virus, host_assignment_phylum$Confidence.score ), ]
host_assignment_phylum_unique <- host_assignment_phylum[ !duplicated(host_assignment_phylum$Virus), ]


## Creating Phylum guild:
RPKM_VLP_by_BacPhylum <- data.frame(matrix(ncol = ncol(RPKM_VLP), nrow=length(unique(host_assignment_phylum_unique$Phylum_new)) ) )
colnames(RPKM_VLP_by_BacPhylum) <- colnames(RPKM_VLP)
row.names(RPKM_VLP_by_BacPhylum) <- unique(host_assignment_phylum_unique$Phylum_new)

for (i in colnames(RPKM_VLP_by_BacPhylum) ) {
  
  for (j in row.names(RPKM_VLP_by_BacPhylum) ) {
    
    RPKM_VLP_by_BacPhylum[row.names(RPKM_VLP_by_BacPhylum)==j,i] <- sum( RPKM_VLP[ c( host_assignment_phylum_unique[ host_assignment_phylum_unique$Phylum_new==j, ]$Virus ) , i] )
    
  }
  
}

## testing differences between mothers and babies:
# filtering for presence in more than 5% of microbiome samples:
RPKM_VLP_by_BacPhylum_filt <- RPKM_VLP_by_BacPhylum[(rowSums(RPKM_VLP_by_BacPhylum!=0) > 0.05*ncol(RPKM_VLP_by_BacPhylum)),  ]
RPKM_VLP_by_BacPhylum_filt <- as.data.frame(t(RPKM_VLP_by_BacPhylum_filt))
# CLR-transformation
my_pseudocount_normal=min(RPKM_VLP_by_BacPhylum_filt[RPKM_VLP_by_BacPhylum_filt!=0])/2
RPKM_VLP_by_BacPhylum_filt_CLR<-decostand(RPKM_VLP_by_BacPhylum_filt, "clr", pseudocount=my_pseudocount_normal)

df<-merge(VLP_metadata, RPKM_VLP_by_BacPhylum_filt_CLR, by="row.names")
row.names(df)<-df$Row.names
df$Row.names=NULL

Prevalent= names (RPKM_VLP_by_BacPhylum_filt_CLR)
names (VLP_metadata)

pheno_list= "Type"

Overall_result =tibble() 

for (Bug in Prevalent){
  if (! Bug %in% colnames(df)){ next }
  #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
  # print (c(Bug, Prevalence))
  Bug2 = paste(c("`",Bug, "`"), collapse="")
  for ( pheno in pheno_list){
    pheno2 = paste(c("`",pheno, "`"), collapse="")
    df$NG_ID[is.na(df[colnames(df) == pheno]) == F] -> To_keep
    df_pheno = filter(df, NG_ID %in% To_keep )
    Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|NEXT_ID)"), collapse="" )) 
    lmer(Model0, df_pheno) -> resultmodel0
    base_model=resultmodel0
    Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",pheno2, "+ (1|NEXT_ID)"), collapse="" ))
    lmer(Model2, df_pheno, REML = F) -> resultmodel2
    M = "Mixed"
    as.data.frame(anova(resultmodel2, base_model))[2,8]->p_simp
    as.data.frame(summary(resultmodel2)$coefficients)[3,1:5] -> Summ_simple
    Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
    rbind(Overall_result, temp_output) -> Overall_result
  }
}

p =as.data.frame(Overall_result)
p$FDR<-p.adjust(p$P, method = "BH")
p$significance_level <- NA
p[p$FDR > 0.05,]$significance_level <- 'ns'
p[p$FDR <= 0.05,]$significance_level <- '*'
p[p$FDR <= 0.01,]$significance_level <- '**'
p[p$FDR <= 0.001,]$significance_level <- '***'

#### ALL DETECTED
RPKM_VLP_by_BacPhylum_plot <- as.data.frame(t(RPKM_VLP_by_BacPhylum_filt_CLR))
RPKM_VLP_by_BacPhylum_plot$Phylum <- row.names(RPKM_VLP_by_BacPhylum_plot)
RPKM_VLP_by_BacPhylum_plot_melt <- melt(RPKM_VLP_by_BacPhylum_plot)
RPKM_VLP_by_BacPhylum_plot_melt$Type <- VLP_metadata$Type[match(RPKM_VLP_by_BacPhylum_plot_melt$variable, VLP_metadata$Short_sample_ID)]
RPKM_VLP_by_BacPhylum_plot_melt$Kingdom <- host_assignment_phylum_unique$Kingdom[match(RPKM_VLP_by_BacPhylum_plot_melt$Phylum, host_assignment_phylum_unique$Phylum_new)]
RPKM_VLP_by_BacPhylum_plot_melt<- RPKM_VLP_by_BacPhylum_plot_melt[order(RPKM_VLP_by_BacPhylum_plot_melt$Phylum),]
RPKM_VLP_by_BacPhylum_plot_melt$significance_level <- p$significance_level[match(RPKM_VLP_by_BacPhylum_plot_melt$Phylum, p$Bug )]

pdf('./04.PLOTS/Phyla_virus_guild_difference.pdf', width=15/2.54, height=20/2.54)
ggplot(RPKM_VLP_by_BacPhylum_plot_melt, aes(value, reorder(Phylum, desc(Phylum)), fill=Type)) +
  geom_boxplot(outlier.shape = NA,alpha=0.5) +
  geom_sina(aes(color=Type), size=0.8,alpha=0.5) +
  theme_bw()+
  theme(axis.text=element_text(size=16), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size=13),
        legend.title = element_text(size=16, face="bold")) +
  geom_stripes(odd = "#33333333", even = "#00000000") + 
  facet_grid(rows = vars(Kingdom), scales = "free_y", space='free_y') + 
  geom_text(aes(label = significance_level, x = 11, y = Phylum), size = 5, angle=270) + 
  labs (y="Predicted host phylum", x="CLR abundance of phages")
dev.off()


##### IS THERE CORRELATION BETWEEN PREVALENCE AND ABUNDANCE OF BACTERIA AND BACTERIOPHAGES?
##### Does it depend on Type?
RPKM_VLP_by_BacPhylum_filt_CLR <- RPKM_VLP_by_BacPhylum_filt_CLR[row.names(RPKM_VLP_by_BacPhylum_filt_CLR) %in% paste0(substr(row.names(microbiome_phyla_filt_CLR), 1,7), 'V'),]
colnames(RPKM_VLP_by_BacPhylum_filt_CLR) <- paste0(colnames(RPKM_VLP_by_BacPhylum_filt_CLR), '_phages')
# list of bacterial phyla that have phages predicted to infect them
Common_phyla <- colnames(microbiome_phyla_filt_CLR)[colnames(microbiome_phyla_filt_CLR) %in% colnames(RPKM_VLP_by_BacPhylum_filt)]
# necessary to merge:
row.names(microbiome_phyla_filt_CLR) <- paste0(substr(row.names(microbiome_phyla_filt_CLR), 1,7), 'V')

bacteria_and_phages <- merge(microbiome_phyla_filt_CLR[,Common_phyla], RPKM_VLP_by_BacPhylum_filt_CLR[,paste0(Common_phyla, '_phages')], by='row.names')
colnames(bacteria_and_phages)[1] <- 'Short_sample_ID'
bacteria_and_phages <- merge(bacteria_and_phages, VLP_metadata, by='Short_sample_ID')

df <- bacteria_and_phages
Prevalent= Common_phyla
pheno_list= "Type"

Overall_result =tibble() 

for (Bug in Prevalent){
  if (! Bug %in% colnames(df)){ next }
  #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
  # print (c(Bug, Prevalence))
  Bug2 = paste(c("`",Bug, "`"), collapse="")
  phylum_phage <- paste0(Bug, '_phages')
  for ( pheno in pheno_list){
    pheno2 = paste(c("`",pheno, "`"), collapse="")
    df$NG_ID[is.na(df[colnames(df) == pheno]) == F] -> To_keep
    df_pheno = filter(df, NG_ID %in% To_keep )
    Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",phylum_phage, "+ (1|NEXT_ID)"), collapse="" )) 
    lmer(Model0, df_pheno) -> resultmodel0
    base_model=resultmodel0
    
    linearModel0=as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",phylum_phage), collapse="" ))
    summary(lm(linearModel0, df_pheno))$adj.r.squared -> R_squared_M0
    
    Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",phylum_phage,"+",pheno2, "+ (1|NEXT_ID)"), collapse="" ))
    lmer(Model2, df_pheno, REML = F) -> resultmodel2
    linearModel2=as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads +",phylum_phage, "+",pheno2), collapse="" ))
    summary(lm(linearModel2, df_pheno))$adj.r.squared -> R_squared_M2
    M = "Mixed"
    as.data.frame(anova(resultmodel2, base_model))[2,8]->p_simp
    as.data.frame(summary(resultmodel2)$coefficients)[3,1:5] -> Summ_simple
    Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple", R_squared_M0=R_squared_M0, R_squared_M2=R_squared_M2) -> temp_output
    rbind(Overall_result, temp_output) -> Overall_result
  }
}

p_bacteria_and_phages =as.data.frame(Overall_result)
p_bacteria_and_phages$FDR<-p.adjust(p_bacteria_and_phages$P, method = "BH")
p_bacteria_and_phages$significance_level <- NA
p_bacteria_and_phages[p_bacteria_and_phages$FDR > 0.05,]$significance_level <- 'ns'
p_bacteria_and_phages[p_bacteria_and_phages$FDR <= 0.05,]$significance_level <- '*'
p_bacteria_and_phages[p_bacteria_and_phages$FDR <= 0.01,]$significance_level <- '**'
p_bacteria_and_phages[p_bacteria_and_phages$FDR <= 0.001,]$significance_level <- '***'

#### as a result, adding Type to the correlation of phylum abundance vs phylum phage abundance increased the explained variance

##### plot of correlations:

correlation_plot <- rbind(data.frame(unname(bacteria_and_phages[,c("Actinobacteria","Actinobacteria_phages","Type")])), 
                          data.frame(unname(bacteria_and_phages[,c("Proteobacteria", "Proteobacteria_phages", "Type")])),
                          data.frame(unname(bacteria_and_phages[,c("Euryarchaeota","Euryarchaeota_phages","Type")])), 
                          data.frame(unname(bacteria_and_phages[,c("Verrucomicrobia", "Verrucomicrobia_phages", "Type")])))
colnames(correlation_plot)[c(1:3)] <- c('Phylum', 'Phages', 'Type')
correlation_plot$Source <- c(rep("Actinobacteria", nrow(bacteria_and_phages)),
                      rep("Proteobacteria", nrow(bacteria_and_phages)),
                      rep("Euryarchaeota", nrow(bacteria_and_phages)),
                      rep("Verrucomicrobia", nrow(bacteria_and_phages)))

pdf('./04.PLOTS/Phyla_virus_guild_correlation.pdf', width=15/2.54, height=12/2.54)
ggplot(correlation_plot, aes(Phylum, Phages, color=Type)) +
  geom_point(alpha=0.6, size=0.5) +
  geom_smooth(method = 'lm') + 
  facet_wrap(~Source, ncol=2, nrow=2, scales = 'free') +
  theme_bw()+
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14,face="bold"),
        strip.text = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14, face="bold")) +
  labs (y="CLR abundance of phages", x="CLR abundance of phyla")
dev.off()


### not used graph
RPKM_VLP_by_BacPhylum['Unknown',] <- colSums(RPKM_VLP) - colSums(RPKM_VLP_by_BacPhylum)
RPKM_VLP_by_BacPhylum_RA <- as.data.frame( t(t(RPKM_VLP_by_BacPhylum)/colSums(RPKM_VLP_by_BacPhylum) )  )
RPKM_VLP_by_BacPhylum_RA$BacPhylum <- row.names(RPKM_VLP_by_BacPhylum_RA)
RPKM_VLP_by_BacPhylum_RA_melt <- melt(RPKM_VLP_by_BacPhylum_RA)

ggplot(RPKM_VLP_by_BacPhylum_RA_melt, aes(BacPhylum, Variable))

ggplot(RPKM_VLP_by_BacPhylum_RA_melt, aes(variable, value, fill=BacPhylum)) + 
  geom_bar(stat = 'identity') + 
  #facet_grid(~Type, scales = "free", space = "free_x") + 
  xlab(label = "Bacteria phyla") + 
  ylab(label = "Number of predicted bacteriophages") + 
  labs(fill='Phyla') +
  theme_bw() + 
  #scale_fill_manual(values=colorBlindGrey8) +
  theme(axis.text.y=element_text(size=16), 
        axis.text.x=element_blank(),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size = 12))

###### OUTPUT #####
write.table(p, '03a.RESULTS/Phylum_phage_Mixed_models.txt', sep='\t', quote=F, row.names=F)
