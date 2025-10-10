########################## NEXT CAZYMES DESCRIPTIVES TIMEPOINT MOTHERS AND INFANTS ###########################################
### AUTHOR:  TRISHLA SINHA, ALEXANDER KURULSHIKOV
### ORIGINAL SCRIPT: 9th June, 2025
### LAST UPDATE: 13th Septmber , 2025 


# Load libraries 
library(tidyverse)
library(data.table)
library(ggplot2)
library(dplyr)
library(lmerTest)
library(pheatmap)
library(vegan)
library(wesanderson)

# Load functions 
Prepare_cazymes = function(File = "/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes/NEXT_cayman_all_17_03_2025.txt", Prevalence_min = 0, Transform = "log" ){
  caz = read_tsv(File)
  caz %>% filter(! feature %in% c("total_reads", "filtered_reads")) -> caz
  #Get unique counts + ambigous counts normalized
  caz %>% select(ID, feature, combined_rpkm ) %>% spread(key = feature, value = combined_rpkm) -> caz
  #Change ID
  caz$ID = caz$ID 
  #Get summary
  summary_as_df <- function(x) {
    result <- summary(x)
    data.frame(stat = names(result), value = as.vector(result))
  }
  summary_df <- caz %>%
    select(-ID) %>%
    purrr::map(summary_as_df) %>%
    bind_rows(.id = "Cazyme") %>% spread(stat, value) %>% as_tibble()
  summary_df = summary_df %>% mutate( `NA's` = ifelse( is.na(`NA's`), 0,  `NA's` ), N = dim(caz)[1], Present= N - `NA's`)
  summary_df = summary_df %>% mutate(Present_perc = Present/N ) %>% mutate(Keep = ifelse(Present_perc>=Prevalence_min, T, F ) )
  #Make NA into 0
  caz %>% mutate_all(~ replace(., is.na(.), 0)) -> caz
  #Remove Cazymes with high NA prop
  caz %>% select(c("ID", summary_df$Cazyme[summary_df$Keep==T] )) -> caz_anal
  if (Transform == "log"){
    Log_transform(caz_anal)  ->caz_anal_tranf
  } else if (Transform == "clr"){
    CLR_transformation(caz_anal)  -> caz_anal_tranf
  } else { caz_anal_tranf = NULL}
  
  return(list("caz" = caz, "caz_filtered" = caz_anal, "caz_tranf"=caz_anal_tranf  ,"Summary_caz" = summary_df ) )
}


Log_transform = function(df){
  colnames(df)[1] = "ID"
  PS = Pseudocount(df)
  df %>% select(-ID) -> df2
  log10(df2 + PS) %>% as_tibble() %>% mutate(ID = df$ID, .before=1) %>% return()
}
Pseudocount = function(df_tara){
  df_tara %>% select(-ID) %>% as.matrix() -> PS
  min(PS[PS!=0])/2 %>% return()
}
CLR_transformation = function( df_tara ){
  df_2 = select(df_tara, -ID)
  df_2 = df_2 + Pseudocount(df_tara)
  matrix_data = as.matrix(df_2)
  
  clr_matrix <- apply(matrix_data, 1, function(row){
    exp(mean(log(row)) ) -> geometric_mean
    log(row / geometric_mean)
  } ) %>% t()
  clr_matrix %>% as_tibble() %>% mutate(ID =df_tara$ID  , .before=1) %>% return()
  
}
Filter_prevalence = function(df_tara, Prevalence = 0.1 ){
  df_tara %>% select(-ID) %>% apply(2, function(x){ mean(x!=0) }  ) -> Prevalences
  Prevalences[Prevalences>=Prevalence] %>% names() -> KEEP
  df_tara %>% select(c("ID",KEEP) ) %>% return()
}


# Load data 
setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
metadata<-read.delim("/Users/trishlasinha/Desktop/LLNEXT/Analysis/metadata/LLNEXT_metadata_15_04_2024.txt")
metadata<- metadata %>%
  filter(!(Type == "mother" & Timepoint_categorical %in% c("M1", "M2")))
metadata$Timepoint_categorical=factor(metadata$Timepoint_categorical, levels = c("P12","P28","B", "W2", "M1", "M2", "M3", "M6", 'M9', "M12"))
metadata$Type=factor(metadata$Type, levels = c("mother", "infant"))
metadata_infants<-metadata[metadata$Type=="infant", ]
#names (metadata_infants)[2]<-"next_id_infant"
metadata_infants$Timepoint_categorical=factor(metadata_infants$Timepoint_categorical, levels = c("W2", "M1", "M2", "M3", "M6", 'M9', "M12"))
metadata_infants$Timepoint_early_late <- dplyr::case_when(
  metadata_infants$Timepoint_categorical %in% c("W2", "M1", "M2", "M3") ~ "early",
  metadata_infants$Timepoint_categorical %in% c("M6", "M9", "M12") ~ "late",
  TRUE ~ NA_character_
)

# PCOA plot CAZYMES MOTHER AND INFANT TIMEPOINTS 
caz_raw<-cazymes_clean[[1]]
caz_raw<-as.data.frame(caz_raw)
row.names(caz_raw)<-caz_raw$ID
caz_raw$ID=NULL
my_pseudocount_normal=min(caz_raw[caz_raw!=0])/2#
distance=vegdist(caz_raw, method = "aitchison", pseudocount=my_pseudocount_normal) 
mypcoa_CLR=cmdscale(distance, k = 20, eig = T)
my_var_CLR=round(mypcoa_CLR$eig*100/sum(mypcoa_CLR$eig),2)[1:20]
barplot (my_var_CLR)
mypcoa_df_CLR=as.data.frame(mypcoa_CLR$points)
names(mypcoa_df_CLR) <- c('PC1','PC2','PC3','PC4','PC5', 'PC6','PC7','PC8','PC9','PC10',
                          'PC11','PC12','PC13','PC14','PC15', 'PC16','PC17','PC18','PC19','PC20')

timepoint_colors <- c("#f90404", "#f78310", "#fbd123","#b5dd88","#41c0b4", "#4397bb", "#eca4c9", "#cb4563","#a42097", "#390962" )

phenos_pc_CLR=merge(metadata,mypcoa_df_CLR, by="row.names")
row.names (phenos_pc_CLR)<-phenos_pc_CLR$Row.names


overall <- ggplot(phenos_pc_CLR , aes(PC1, PC2, color = Type)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(group = Type, fill = Type, color = Type),
               type = "norm", linetype = 2, geom = "polygon", alpha = 0.05, show.legend = F) +
  xlab(paste("PCo1=", round(my_var_CLR[1], digits = 2), "%", sep = "")) +
  ylab(paste("PCo2=", round(my_var_CLR[2], digits = 2), "%", sep = "")) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5)) +
  scale_color_manual(name=NULL,
                     breaks = c("mother", "infant"),
                     labels = c("mother", "infant"),
                     values = wes_palette("Moonrise3", n = 2)) +
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

overall <- ggExtra::ggMarginal(overall, type = "histogram", groupColour = F, groupFill = TRUE,
                               xparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'),
                               yparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'))

overall

overall_timepoint <- ggplot(phenos_pc_CLR, aes(PC1, PC2, color = Timepoint_categorical)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(group = Timepoint_categorical, fill = Timepoint_categorical, color = Timepoint_categorical),
               type = "norm", linetype = 2, geom = "polygon", alpha = 0.05, show.legend = F) +
  xlab(paste("PCo1=", round(my_var_CLR[1], digits = 2), "%", sep = "")) +
  ylab(paste("PCo2=", round(my_var_CLR[2], digits = 2), "%", sep = "")) +
  scale_fill_manual(values = timepoint_colors) +  
  scale_color_manual(values = timepoint_colors) +  
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

overall_timepoint <- ggExtra::ggMarginal(overall_timepoint, type = "histogram", groupColour = F, groupFill = TRUE,
                               xparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'),
                               yparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'))

overall_timepoint

phenos_pc_CLR_infants <-phenos_pc_CLR[phenos_pc_CLR$Type=="infant",]
infant_colors<-c("#b5dd88","#41c0b4", "#4397bb", "#eca4c9", "#cb4563","#a42097", "#390962")
overall_timepoint_infants <- ggplot(phenos_pc_CLR_infants, aes(PC1, PC2, color = Timepoint_categorical)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(group = Timepoint_categorical, fill = Timepoint_categorical, color = Timepoint_categorical),
               type = "norm", linetype = 2, geom = "polygon", alpha = 0.05, show.legend = F) +
  xlab(paste("PCo1=", round(my_var_CLR[1], digits = 2), "%", sep = "")) +
  ylab(paste("PCo2=", round(my_var_CLR[2], digits = 2), "%", sep = "")) +
  scale_fill_manual(values = infant_colors) +  
  scale_color_manual(values = infant_colors) +  
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key = element_rect(fill = NA),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

overall_timepoint_infants <- ggExtra::ggMarginal(overall_timepoint_infants, type = "histogram", groupColour = F, groupFill = TRUE,
                                         xparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'),
                                         yparams = list(bins = 60, alpha = 0.8, position = 'identity', color = 'white'))

overall_timepoint_infants


# Loading & processing cazyme data 

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
caz_tables = Prepare_cazymes(File = "NEXT_cayman_infant.txt", Prevalence_min = 0.3, Transform = "NULL" )
caz_filtered = as.data.frame (caz_tables[["caz_filtered"]])
row.names(caz_filtered)<-caz_filtered$ID
caz_filtered$ID<-NULL

# Loading microbiome data 
metaphlan = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLEAN_10_07_2023.txt")
metaphlan_clr = read.delim("~/Desktop/LLNEXT/Analysis/taxa/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt")
metaphlan_infants = metaphlan[match(as.character(metadata_infants$NG_ID),rownames(metaphlan)),grep("t__",colnames(metaphlan))]
colnames(metaphlan_infants) = sub(".*s__","",colnames(metaphlan_infants))
metaphlan_infants = metaphlan_infants[,colnames(metaphlan_clr)] # Filtering to only the columns in the filtered file 
metaphlan_infants <- metaphlan_infants[rownames(metaphlan_infants) %in% rownames(caz_filtered), ]




# Data transformations 
cazyme_alr = do_clr_externalWeighting(caz_filtered,metaphlan_infants)
cazyme_alr<-as.data.frame(cazyme_alr)
cazyme_alr_null <-nullify_zeros(cazyme_alr,caz_filtered)


# Association of timepoint with cazyme profile in infants 
Timepoint_Cazyme_infants <- mixed_models_without_time_correction(metadata_infants, 
                                              "NG_ID", 
                                              cazyme_alr_null, 
                                              c("Timepoint_categorical", "Timepoint_early_late"))

setwd("/Users/trishlasinha/Desktop/LLNEXT/Analysis/cazymes")
write.table(Timepoint_Cazyme_infants, "association_early_late_timepoint_infants.txt", sep = "\t")


wide_data <- reshape2::dcast(Timepoint_Cazyme_infants, 
                             Bug ~ Feature, 
                             value.var = "t value")


colnames(wide_data) <- gsub("Timepoint_categorical", "", colnames(wide_data))


filtered_all_wide <- wide_data %>%
  select(Bug, M1, M2, M3, M6, M9, M12)

wide_data_matrix <- as.matrix(filtered_all_wide[, -1])  
row.names(wide_data_matrix) <- filtered_all_wide$Bug

custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)


pdf("cazyme_infant_timepoints.pdf", width = 5, height = 12)
pheatmap(
  wide_data_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  color = custom_colors,
  main = "",
  fontsize_col = 12,
  fontsize_row = 3,
  fontsize_number = 12,
  angle_col = 0  
)
dev.off()


early_late_colors <- c(
  "mother"   = "white", # soft green
  "early"  = "#4397bb", # light-medium green
  "late" = "#390962"  # very deep green
)

