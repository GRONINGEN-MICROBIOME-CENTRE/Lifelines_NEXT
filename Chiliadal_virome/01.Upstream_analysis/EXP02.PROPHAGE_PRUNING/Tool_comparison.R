setwd("~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/04.RAW_DATA/")

##########################################
# Benchmarking prophage pruning tools
##########################################

##############################
# Loading libraries
##############################
library(stringr)
library(reshape2)
library(ggplot2)
library(ggforce)
library(ggtext)
##############################
# Functions
##############################

##############################
# Input data
##############################

real_coordinates <- read.table('04.Virome_discovery/prophage_benchmark/Prophage_benchmark_coordinates.txt', sep='\t', header=T)

# Cenote-Taker3
CT3_summary <- read.table('04.Virome_discovery/prophage_benchmark/CenoteTaker3_virus_summary.tsv', sep='\t', header = T)
CT3_summary$contig_id <- CT3_summary$contig

CT3_prune <- read.table('04.Virome_discovery/prophage_benchmark/CenoteTaker3_prune_summary.tsv', sep='\t', header=T)
CT3_prune$contig_id <- gsub('@$','',gsub('NaN', '', paste0(CT3_prune$contig, "@", CT3_prune$chunk_name)))

CT3 <- merge(CT3_summary[,c("contig_id", "input_name")], CT3_prune, by='contig_id')
CT3$contig_id <- NULL
colnames(CT3)[c(1,6,7)] <- c("contig_id", "Begin", "End")

# CheckV
CHV_prune <- read.table('04.Virome_discovery/prophage_benchmark/CHV_prophage_coordinates.txt', sep='\ ', header=F)
CHV_prune$contig_id <- gsub('_1$', '', CHV_prune$V1)
CHV_prune$Begin <- gsub('-.*','',CHV_prune$V2)
CHV_prune$End <- gsub('.*-|\\/.*','',CHV_prune$V2)

CHV_quality <- read.table('04.Virome_discovery/prophage_benchmark/quality_summary.tsv', sep='\t', header=T)

CHV <- merge(CHV_prune[,c("contig_id", "Begin", "End")], 
             CHV_quality[,c("contig_id", "contig_length", "provirus", "proviral_length", 
                            "gene_count", "viral_genes", "host_genes", "checkv_quality")], by='contig_id', all.y = T)
CHV <- CHV[CHV$provirus=='Yes',]
CHV$Begin <- as.numeric(CHV$Begin)
CHV$End <- as.numeric(CHV$End)
# geNomad
GND_prune <- read.table('04.Virome_discovery/prophage_benchmark/all_positive_check_provirus.tsv', sep='\t', header=T)

GND_summary <- read.table('04.Virome_discovery/prophage_benchmark/all_positive_check_virus_summary.tsv', sep='\t', header = T)

GND <- GND_prune[,c("source_seq", "start", "end", "length")]

colnames(GND)[c(1:3)] <- c("contig_id", "Begin", "End")

# PhageBoost
PHB_prune <- read.table('04.Virome_discovery/prophage_benchmark/PhageBoost_output.txt', sep='\t', header=F)
PHB_prune$contig <- str_extract(PHB_prune$V1, "phage\\d+")
PHB_prune$Begin <- sapply(str_split(gsub('.*phage[0-9]+_','',PHB_prune$V1), "_"), `[[`, 1)
PHB_prune$End <- sapply(str_split(gsub('.*phage[0-9]+_','',PHB_prune$V1), "_"), `[[`, 2)
PHB_prune$tmp <- gsub('_phage[0-9]+_.*','',PHB_prune$V1)

CHV_quality$tmp <- gsub('\\..*','',CHV_quality$contig_id)

PHB <- merge(PHB_prune, CHV_quality[,c("tmp", "contig_id")], by='tmp', all.x=T)
PHB$tmp <- NULL
PHB$V1 <- NULL
PHB <- PHB[,c("contig_id", "contig", "Begin", "End")]

PHB$Begin <- as.numeric(PHB$Begin)
PHB$End <- as.numeric(PHB$End)

# Phigaro
PHI <- read.table('04.Virome_discovery/prophage_benchmark/all_positive_check.phigaro.tsv', sep='\t', header=T)
colnames(PHI)[c(1,3,4)] <- c("contig_id", "Begin", "End")

# VIBRANT
VIB <- read.table('04.Virome_discovery/prophage_benchmark/VIBRANT_integrated_prophage_coordinates_all_positive_check.AA.tsv', sep = '\t', header=T)
VIB <- VIB[,c("scaffold", "fragment", "nucleotide.start", "nucleotide.stop", "nucleotide.length")]

colnames(VIB)[c(1, 3, 4)] <- c("contig_id", "Begin", "End")

# VirSorter2
VS2 <- read.table('04.Virome_discovery/prophage_benchmark/final-viral-boundary.tsv', sep='\t', header = T)

colnames(VS2)[c(1, 4, 5)] <- c("contig_id", "Begin", "End")

# pilot metadata, to check lytic viruses pruning by CheckV:
pilot_metadata <- read.table('~/Desktop/Projects_2022/NEXT_pilot_FUP/02.CLEAN_DATA/VLP_viral_contigs_metadata_noNA.txt', sep='\t', header = T)

##############################
# ANALYSIS
##############################
to_compare <- c('CHV', 'CT3', 'GND', 'PHB', 'PHI', 'VIB', 'VS2')
comparison <- real_coordinates
comparison$length <- comparison$End - comparison$Begin

# loop to parse the coordinates and length of detected prophages
for (tool in to_compare) {
  
  comparison[ , paste0( tool, '_Begin' ) ] <- NA
  comparison[ , paste0( tool, '_End' ) ] <- NA
  comparison[ , paste0( tool, '_length' ) ] <- NA
  
  for (i in real_coordinates$contig_id) {
    
    if ( length( get(tool)[get(tool)$contig_id==i,]$Begin ) != 1 ) {
      
      if ( length( get(tool)[get(tool)$contig_id==i,]$Begin ) == 0) {
        
        comparison[comparison$contig_id==i,paste0( tool, '_Begin' )] <- NA
        comparison[comparison$contig_id==i,paste0( tool, '_End' )] <- NA
        comparison[comparison$contig_id==i,paste0( tool, '_length' )] <- NA
        
      } else {
        
        chooser <- get(tool)[get(tool)$contig_id==i,]$Begin
        index <- which.min( abs(chooser - get(tool)[get(tool)$contig_id==i,]$Begin) )
        
        comparison[comparison$contig_id==i,paste0( tool, '_Begin' )] <- get(tool)[get(tool)$contig_id==i,]$Begin[index]
        comparison[comparison$contig_id==i,paste0( tool, '_End' )] <- get(tool)[get(tool)$contig_id==i,]$End[index]
        comparison[comparison$contig_id==i,paste0( tool, '_length' )] <- (get(tool)[get(tool)$contig_id==i,]$End[index] - get(tool)[get(tool)$contig_id==i,]$Begin[index])
      }
      
    } else {
      
      comparison[comparison$contig_id==i,paste0( tool, '_Begin' )] <- get(tool)[get(tool)$contig_id==i,]$Begin
      comparison[comparison$contig_id==i,paste0( tool, '_End' )] <- get(tool)[get(tool)$contig_id==i,]$End
      comparison[comparison$contig_id==i,paste0( tool, '_length' )] <- get(tool)[get(tool)$contig_id==i,]$End - get(tool)[get(tool)$contig_id==i,]$Begin
    }
    
  }
  
  
}

comparison_dif <- cbind(comparison$Candidate.prophage, 
                        comparison[,grep('_Begin', colnames(comparison))] - comparison$Begin,
                        comparison[,grep('_End', colnames(comparison))] - comparison$End)

apply(comparison_dif[ , grep('Begin', colnames(comparison_dif)) ],2,median, na.rm=T)
apply(comparison_dif[ , grep('End', colnames(comparison_dif)) ],2,median, na.rm=T)

dif_melt <- melt(comparison_dif)
dif_melt$Tail <- gsub('.*_', '', dif_melt$variable)
dif_melt$Tool <- gsub('_.*', '', dif_melt$variable)
dif_melt$value_normalized <- dif_melt$value/1000
dif_melt$value_log <- log10(abs(dif_melt$value)+1)

# What is the difference in prophage coordinates (start and end)?
pdf('../05.PLOTS/02.EXP_PROPHAGE_PRUNING/Dif_positive_control_coordiantes.pdf', width=12/2.54, height=8/2.54)
ggplot(dif_melt, aes(Tool, abs(value), color=Tool)) + 
  geom_boxplot(aes(fill=Tool), outlier.colour = NA, width=0.5, alpha=0.4) +
  geom_sina(aes(color=Tool), maxwidth=0.5, alpha=1, size=0.6) +
  facet_wrap(~Tail) +
  labs(fill='Tool', color='Tool', y="Difference in bp", x="", title = "Difference in confirmed* prophages coordinates") +
  theme_bw() +
  theme(axis.text =element_text(size=7),
        axis.title = element_text(size=9, face="bold"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.background = element_rect(color="black", fill=NA), 
        plot.title = element_text(size=9, face="bold", hjust = 0.5),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.position = 'none') + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the difference in prophage coordinates (start and end) expressed in log10?
pdf('../05.PLOTS/02.EXP_PROPHAGE_PRUNING/Dif_positive_control_coordiantes_log10.pdf', width=12/2.54, height=8/2.54)
ggplot(dif_melt, aes(Tool, abs(value_log), color=Tool)) + 
  geom_boxplot(aes(fill=Tool), outlier.colour = NA, width=0.5, alpha=0.4) +
  geom_sina(aes(color=Tool), maxwidth=0.5, alpha=1, size=0.6) +
  facet_wrap(~Tail) +
  labs(fill='Tool', color='Tool', y="Log<sub>10</sub> difference in bp", x="", 
       title = "Difference in confirmed* prophages coordinates") +
  theme_bw() +
  theme(axis.text =element_text(size=7),
        axis.title.y = element_markdown(),
        axis.title = element_text(size=9, face="bold"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.background = element_rect(color="black", fill=NA), 
        plot.title = element_text(size=9, face="bold", hjust = 0.5),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.position = 'none') + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()


comparison_length <- cbind(comparison$Candidate.prophage, comparison[,grep('_length', colnames(comparison))] - comparison$length)
apply(comparison_length [ , grep('length', colnames(comparison_length )) ],2,median, na.rm=T)

length_melt <- melt(comparison_length)
length_melt$Tool <- gsub('_.*', '', length_melt$variable)

# What is the difference in prophage length expressed in bp?
pdf('../05.PLOTS/02.EXP_PROPHAGE_PRUNING/Dif_positive_control_length.pdf', width=8/2.54, height=8/2.54)
ggplot(length_melt, aes(Tool, abs(value), color=Tool)) + 
  geom_boxplot(aes(fill=Tool), outlier.colour = NA, width=0.5, alpha=0.4) +
  geom_sina(aes(color=Tool), maxwidth=0.5, alpha=1, size=0.6) +
  labs(fill='Tool', color='Tool', y="Difference in bp", x="", title = "Difference in confirmed* prophages length") +
  theme_bw() +
  theme(axis.text =element_text(size=7),
        axis.title = element_text(size=9, face="bold"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.background = element_rect(color="black", fill=NA), 
        plot.title = element_text(size=9, face="bold", hjust = 0.5),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.position = 'none') + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the difference in prophage length expressed in log10?
pdf('../05.PLOTS/02.EXP_PROPHAGE_PRUNING/Dif_positive_control_length_log10.pdf', width=8/2.54, height=8/2.54)
ggplot(length_melt, aes(Tool, log(abs(value) + 1), color=Tool)) + 
  geom_boxplot(aes(fill=Tool), outlier.colour = NA, width=0.5, alpha=0.4) +
  geom_sina(aes(color=Tool), maxwidth=0.5, alpha=1, size=0.6) +
  labs(fill='Tool', color='Tool', y="Log<sub>10</sub> difference in bp", x="", 
       title = "Difference in confirmed* prophages length") +
  theme_bw() +
  theme(axis.text =element_text(size=7),
        axis.title.y = element_markdown(),
        axis.title = element_text(size=9, face="bold"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.background = element_rect(color="black", fill=NA), 
        plot.title = element_text(size=9, face="bold", hjust = 0.5),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.position = 'none') + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# Sensitivity vs specificity
detection <- data.frame(matrix(NA, nrow=length(to_compare), ncol=3))
colnames(detection) <- c('Tool', 'N_detected_CP', 'N_detected_FP')
detection$Tool <- to_compare

for (tool in to_compare) {
  
  detection[detection$Tool==tool,"N_detected_CP"] <- length(na.omit(comparison[,paste0(tool, '_length')]))
  detection[detection$Tool==tool,"N_detected_FP"] <- -dim(get(tool))[1]
  
}


detection_long <- tidyr::pivot_longer(detection, cols = starts_with("N_detected_"))

pdf('../05.PLOTS/02.EXP_PROPHAGE_PRUNING/Compare_N_detected_CP_vs_FP.pdf', width=8/2.54, height=8/2.54)
ggplot(detection_long, aes(x = Tool, y = value, fill = name, color=name)) +
  geom_bar(stat = "identity") +
  labs(fill='Detected prophages', y="N detected prophages", x="", 
       title = "Detection of confirmed* vs novel** prophages") +
  scale_fill_manual(values = c("N_detected_CP" = "#A3C9AA", "N_detected_FP" = "#C68484"), labels=c('Confirmed', 'All')) +
  scale_color_manual(values = c("N_detected_CP" = "#4a5c4e", "N_detected_FP" = "#754f4f")) +
  scale_y_continuous(breaks = c(20, 0, -20, -40),
                     labels = c(20, 0, 20, 40)) +
  geom_hline(aes(yintercept=17), linewidth=1, linetype="dashed") +
  geom_hline(aes(yintercept=-17), linewidth=1, linetype="dashed") +
  theme_bw() +
  theme(axis.text =element_text(size=7),
        axis.title.y = element_markdown(),
        axis.title = element_text(size=9, face="bold"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.background = element_rect(color="black", fill=NA), 
        plot.title = element_text(size=9, face="bold", hjust = 0.5),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.position = 'bottom') + 
  guides(color = "none") + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What about predicted complete phage genomes that are not predicted to be temperate?
nrow(pilot_metadata[!is.na(pilot_metadata$Pred_phatyp) & 
                     pilot_metadata$checkv_quality=='Complete' & 
                     pilot_metadata$Pred_phatyp=='virulent' & 
                     pilot_metadata$provirus_checkV=='Yes',])
# 12 lytic (according to PhaTYPE) phages that were contaminated by bacterial sequences according to CheckV
nrow(pilot_metadata[!is.na(pilot_metadata$Pred_phatyp) & 
                      pilot_metadata$checkv_quality=='Complete' & 
                      pilot_metadata$Pred_phatyp=='virulent',])
# in total, there are 314 complete lytic viruses => 4% might be overtrimmed

nrow(pilot_metadata[pilot_metadata$checkv_quality=='Complete' & 
                      pilot_metadata$temperate==0 & 
                      pilot_metadata$provirus_checkV=='Yes',])
# 3 lytic (according to temperate pVOGs) phages that were contaminated by bacterial sequences according to CheckV
nrow(pilot_metadata[pilot_metadata$checkv_quality=='Complete' & 
                      pilot_metadata$temperate==0,])
# in total, there are 247 complete lytic viruses => 1% might be overtrimmed
nrow(pilot_metadata[pilot_metadata$checkv_quality=='Complete' & 
                      pilot_metadata$CrAss==1 & 
                      pilot_metadata$provirus_checkV=='Yes',])

nrow(pilot_metadata[pilot_metadata$checkv_quality=='Complete' & 
                      pilot_metadata$CrAss==1,])

Diff_predictions <- data.frame(
  ADD = c('PhaTYPE', 'PhaTYPE', 'pVOG', 'pVOG'),
  Category = c("Provirus", "Pure", "Provirus", "Pure"),
  Count = c(12, 314 - 12, 3, 247) # 12 red, rest blue
)

pdf('../05.PLOTS/02.EXP_PROPHAGE_PRUNING/Compare_N_trimmed_complete_lytic_CheckV.pdf', width=6/2.54, height=10/2.54)
ggplot(Diff_predictions, aes(x=ADD, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#7E1717", "#068DA9")) +
  labs(x = "Tool predicitng lifestyle", 
       y = "N complete lytic viruses", 
       title = "Putative host comtamination\nin lytic sequences",
       fill = "Host contamination") +
  theme_bw() +
  theme(axis.text =element_text(size=7),
        axis.title.y = element_markdown(),
        axis.title = element_text(size=9, face="bold"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.background = element_rect(color="black", fill=NA), 
        plot.title = element_text(size=9, face="bold", hjust = 0.5),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.position = 'bottom') + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))

dev.off()

#### ! Sometimes, CheckV overtrims seemingly complete genomes (those with high confidence DTRs). 
#### Current solution: since we keep all extended contigs & the overall completeness drop is not substantial,
#### keep overtrimmed sequences (even though they are often VC representatives) & later, in case it is needed
#### to follow-up something with them, try to curate these genomes, there are not that many of them (10-20 of those)
prepruned_quality <- read.table('04.RAW_DATA/04.Virome_discovery/prophage_benchmark/prepruned_contigs_quality_summary.tsv', sep='\t', header=T)
extended_tof <- extended_tof <- read.table("04.RAW_DATA/04.Virome_discovery/Extended_table_of_origin", sep='\t', header=T)
prepruned_quality$new_quality <- extended_tof$checkv_quality[match(prepruned_quality$contig_id, extended_tof$POST_CBR_CID)]
prepruned_quality$new_completeness <- extended_tof$completeness[match(prepruned_quality$contig_id, extended_tof$POST_CBR_CID)]
