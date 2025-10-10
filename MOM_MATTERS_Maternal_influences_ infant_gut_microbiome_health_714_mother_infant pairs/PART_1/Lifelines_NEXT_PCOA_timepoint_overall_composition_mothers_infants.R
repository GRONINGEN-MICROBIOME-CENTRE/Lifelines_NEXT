#################### Associating PC's with Timepoint ####################
# Author: Trishla Sinha 
# Last update: 31st of July, 2024

mixed_models_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      Model0 = as.formula(paste( c(Bug2,  " ~ dna_conc + clean_reads_FQ_1 + BATCH_NUMBER+(1|NEXT_ID)"), collapse="" )) 
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      Model2 = as.formula(paste( c(Bug2,  " ~ dna_conc  + clean_reads_FQ_1 + BATCH_NUMBER+",pheno2, "+ (1|NEXT_ID)"), collapse="" ))
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}

# Load metadata and phenotypes
metadata<-read.delim("~/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")

# Merging files and selecting infant relevant data
metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)],  #convert character columns to factors
                                                   as.factor)
row.names(metadata)<-metadata$NG_ID
metadata$BATCH_NUMBER<-as.factor(metadata$BATCH_NUMBER)
metadata_infants<-metadata[metadata$Type=="infant", ] # n=2939
metadata_mothers<-metadata[metadata$Type=="mother", ] # n=1638
metadata_mothers <- subset(metadata_mothers, Timepoint_categorical != "M1" & Timepoint_categorical != "M2")
metadata_mothers$Timepoint_categorical <- factor(metadata_mothers$Timepoint_categorical, 
                                                 levels = c("P12", "P28", "B", "M3"))

taxa<-read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
mother_taxa<-taxa[row.names(taxa)%in% rownames(metadata_mothers),] 
mother_taxa<-mother_taxa[match(row.names(metadata_mothers),row.names(mother_taxa)),]

# Filtering data 
mother_NEXT_ID<-metadata_mothers %>%
  select(NEXT_ID)
mother_taxa_all<-merge(mother_NEXT_ID,mother_taxa, by="row.names" )
row.names(mother_taxa_all)<-mother_taxa_all$Row.names
mother_taxa_all$Row.names=NULL
unique_counts <- sapply(mother_taxa_all, function(x) length(unique(mother_taxa_all$NEXT_ID[x >0.001]))) # Here I am counting the unique elements in the NEXT_ID column where the corresponding value in each column (i.e., x) is greater than the given cut-off. 
mother_taxa_all_filt <- mother_taxa_all[, unique_counts >= 0.3*length(unique(mother_taxa_all$NEXT_ID)) ] # Setting a 20% cut-off on prevalence (3rd august version) 
mother_taxa_all_filt$NEXT_ID=NULL

mother_taxa_SGB<-mother_taxa_all_filt[,grep("t__",colnames(mother_taxa_all_filt))]
my_pseudocount_normal=min(mother_taxa_SGB[mother_taxa_SGB!=0])/2# 

distance_mother=vegdist(mother_taxa_SGB, method = "aitchison", pseudocount=my_pseudocount_normal) 
mypcoa_CLR=cmdscale(distance_mother, k = 20, eig = T)
my_var_CLR=round(mypcoa_CLR$eig*100/sum(mypcoa_CLR$eig),2)[1:20]
barplot (my_var_CLR)
mypcoa_df_CLR=as.data.frame(mypcoa_CLR$points)
names(mypcoa_df_CLR) <- c('PC1','PC2','PC3','PC4','PC5', 'PC6','PC7','PC8','PC9','PC10',
                          'PC11','PC12','PC13','PC14','PC15', 'PC16','PC17','PC18','PC19','PC20')

time_on_pc <- mixed_models_taxa(metadata_mothers, 
                                        "NG_ID", 
                                        mypcoa_df_CLR, 
                                        c("Timepoint_categorical"))
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/results/beta_diversity/")
write.table (time_on_pc, "timepoint_beta_diversity.txt", sep = "\t", row.names = F)


