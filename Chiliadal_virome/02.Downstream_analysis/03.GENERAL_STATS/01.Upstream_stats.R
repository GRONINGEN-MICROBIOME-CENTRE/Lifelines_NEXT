setwd('~/Desktop/Projects_2021/NEXT_virome/09.DATA_ANALYSIS/')

##########################################
# Raw and quality-trimmed reads stats
##########################################

##############################
# Loading libraries
##############################
library(dplyr)
library(purrr)
library(reshape2)
library(ggplot2)
library(ggforce)
library(MetBrewer)
library(tidyr)
library(tibble)
library(UpSetR)
library(scales)
##############################
# Functions
##############################
my_theme <-   theme_bw() + 
        theme(axis.title = element_text(size=12,face="bold"),
        axis.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=12,face="bold"),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-3,-10))


my_theme2 <- theme_bw() +
        theme(axis.text =element_text(size=10),
        axis.title = element_text(size=10, face="bold"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.background = element_rect(color="black", fill=NA),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.position = 'bottom')

my_theme3 <- theme_bw() +
        theme(axis.text.x =element_text(size=10), 
        axis.text.y =element_text(size=10), 
        axis.title=element_blank(),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.background = element_rect(color="black", fill=NA),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8, face="bold"),
        legend.position = 'bottom') 
##############################
# Input data
##############################
Demuth <- met.brewer("Demuth")
Monet <- met.brewer('Monet')

metadata <- read.delim('01.METADATA/Chiliadal_metadata_ver_01_06042023.txt', sep='\t')

N_raw_reads <- read.table('04.RAW_DATA/01.QC_tuning/N_raw_reads', sep='\ ', header=F)
colnames(N_raw_reads) <- c('Sequencing_ID_VLP', 'N_raw_reads')
N_raw_reads$N_raw_reads <- N_raw_reads$N_raw_reads*2
N_raw_reads$Type <- metadata$Type_VLP[match(N_raw_reads$Sequencing_ID_VLP, metadata$Sequencing_ID_VLP)]

N_postkn_reads <- read.table('04.RAW_DATA/02.Read_stats/N_post_kneaddata_reads', sep='\t', header=F)
colnames(N_postkn_reads) <- c('Sequencing_ID_VLP', 'N_postkn_reads')
N_postkn_reads$Type <- metadata$Type_VLP[match(N_postkn_reads$Sequencing_ID_VLP, metadata$Sequencing_ID_VLP)]

ViromeQC <- read.table('04.RAW_DATA/02.Read_stats/viromeqc_all_stat.txt', sep='\t', header=T)
ViromeQC$Type <- metadata$Type_VLP[match(ViromeQC$Sample, metadata$Sequencing_ID_VLP)]
colnames(ViromeQC)[1] <- 'Sequencing_ID_VLP'

N_duplicated <- read.table('04.RAW_DATA/02.Read_stats/N_duplicated_reads', sep='\t', header=F)
colnames(N_duplicated) <- c('Sequencing_ID_VLP', 'N_duplicates')
N_duplicated$Type <- metadata$Type_VLP[match(N_duplicated$Sequencing_ID_VLP, metadata$Sequencing_ID_VLP)]

N_human <- read.table('04.RAW_DATA/02.Read_stats/N_humane_reads', sep='\t', header=F)
colnames(N_human) <- c('Sequencing_ID_VLP', 'N_human_reads')
N_human$Type <- metadata$Type_VLP[match(N_human$Sequencing_ID_VLP, metadata$Sequencing_ID_VLP)]

Bases <- read.table('04.RAW_DATA/02.Read_stats/N_bases', sep = ':', header=F)
sum(Bases$V3) #6.435068e+12

contigs_stat <- read.table('04.RAW_DATA/02.Read_stats/Contigs_stat', sep='\t', header=T)
contigs_stat$Type <- metadata$Type_VLP[match(contigs_stat$Assembly, metadata$Sequencing_ID_VLP)]
colnames(contigs_stat)[1] <- "Sequencing_ID_VLP"

table_of_origin <- read.table("04.RAW_DATA/04.Virome_discovery/table_of_origin", sep='\t', header=T)
table_of_origin$length <- as.numeric(gsub("^.*_NODE_\\d+_length_|_cov_\\d+\\.\\d+$", "", table_of_origin$V1))
table_of_origin$SAMPLE_ID <- gsub("_NODE.*", "", table_of_origin$V1)

N_disc_viruses <- data.frame(table(table_of_origin$SAMPLE_ID))
colnames(N_disc_viruses) <- c("Sequencing_ID_VLP", "N_disc_viruses")

Length_disc_viruses <- aggregate(table_of_origin[,"length"], list(table_of_origin[,"SAMPLE_ID"]), FUN = sum)
colnames(Length_disc_viruses) <- c("Sequencing_ID_VLP", "Length_disc_viruses")

alignment_to_all <- read.table('04.RAW_DATA/02.Read_stats/bowtie2.perc.alignment.all.txt', sep='\t', header=F)
colnames(alignment_to_all) <- c("Sequencing_ID_VLP", "all_perc_align")

alignment_to_1kbp <- read.table('04.RAW_DATA/02.Read_stats/bowtie2.perc.alignment.1kbp.txt', sep='\t', header=F)
colnames(alignment_to_1kbp) <- c("Sequencing_ID_VLP", "1kbp_perc_align")

alignment_to_own_vir <- read.table('04.RAW_DATA/02.Read_stats/bowtie2.perc.alignment.own.vir.txt', sep='\t', header=F)
colnames(alignment_to_own_vir) <- c("Sequencing_ID_VLP", "own_vir_perc_align")

pilot_vlp_vqc <- read.table('04.RAW_DATA/02.Read_stats/viromeqc_vlp.txt', sep='\t', header=T)
pilot_fsk_vqc <- read.table('04.RAW_DATA/02.Read_stats/viromeqc_total_microbiome.txt', sep='\t', header=T)
pilot_exp_vqc <- rbind(pilot_vlp_vqc, pilot_fsk_vqc)
ViromeQC_pilot <- melt(pilot_exp_vqc[,c(1,4:10)])

extended_tof <- read.table("04.RAW_DATA/04.Virome_discovery/Extended_table_of_origin", sep='\t', header=T)

alignment_to_ext_vir <- read.table('04.RAW_DATA/02.Read_stats/bowtie2.perc.alignment.ext.vir.txt', sep='\t', header=F)
colnames(alignment_to_ext_vir) <- c("Sequencing_ID_VLP", "ext_vir_perc_align")

alignment_to_ext_prun_vir <- read.table('04.RAW_DATA/02.Read_stats/bowtie2.perc.alignment.ext.prun.vir.txt', sep='\t', header=F)
colnames(alignment_to_ext_prun_vir) <- c("Sequencing_ID_VLP", "ext_vir_prun_perc_align")

alignment_to_w_neg_der95 <- read.table('04.RAW_DATA/02.Read_stats/bowtie2.perc.alignment.w_neg_der95.txt', sep='\t', header = F)
colnames(alignment_to_w_neg_der95) <- c("Sequencing_ID_VLP", "w_neg_der95")

Kandinsky <- met.brewer("Kandinsky")
##############################
# ANALYSIS
##############################

# upstream_metadata:
df_list <- list(metadata, N_raw_reads[,c("Sequencing_ID_VLP", "N_raw_reads")], 
                N_postkn_reads[,c("Sequencing_ID_VLP", "N_postkn_reads")],
                N_human[,c("Sequencing_ID_VLP", "N_human_reads")],
                N_duplicated[,c("Sequencing_ID_VLP", "N_duplicates")],
                ViromeQC[,c("Sequencing_ID_VLP","Reads_HQ", "SSU.rRNA.alignment.rate", "LSU.rRNA.alignment.rate", "Bacterial_Markers.alignment.rate", "total.enrichmnet.score")],
                contigs_stat[,c("Sequencing_ID_VLP", "contigs.....0.bp.", "contigs.....1000.bp.",
                                "Total.length.....1000.bp.", "N50")],
                alignment_to_all,
                N_disc_viruses, Length_disc_viruses,
                alignment_to_1kbp, alignment_to_own_vir, alignment_to_ext_vir, alignment_to_ext_prun_vir,
                alignment_to_w_neg_der95) 

VLP_metadata <- df_list %>% reduce(full_join, by='Sequencing_ID_VLP')
VLP_metadata[VLP_metadata$Sequencing_ID_VLP %in% c("CHV199906F12","CHV200013F12"),"Type_VLP"] <- 'NC'

# adding alternative enrichment calculation (given the assumption that all extra 
# alignment to SSU and LSU come from multiple copy ribosomal genes and their transcripts,
# and using only single-copy bacterial markers alignment rates):
VLP_metadata$sc.enrichmnet.score <- 0.700892560801689/VLP_metadata$Bacterial_Markers.alignment.rate

VLP_metadata$perc_viral_contigs <- VLP_metadata$N_disc_viruses/VLP_metadata$contigs.....1000.bp.*100
VLP_metadata$perc_viral_length <- VLP_metadata$Length_disc_viruses/VLP_metadata$Total.length.....1000.bp.*100

# N extended contigs ps

extended_tof$Sequencing_ID_VLP <- gsub("_.*", "", extended_tof$Original_CID)
N_extended_ps <- as.data.frame(table(extended_tof[extended_tof$COB_status=="E1",c("Sequencing_ID_VLP", "COB_status")]))
colnames(N_extended_ps)[grep('Freq', colnames(N_extended_ps))] <- "N_extended_contigs"

VLP_metadata <- reduce(list(VLP_metadata, N_extended_ps[,c("Sequencing_ID_VLP", "N_extended_contigs")]), full_join, by="Sequencing_ID_VLP")
VLP_metadata[is.na(VLP_metadata$N_extended_contigs),]$N_extended_contigs <- 0

Length_ext_viruses <- aggregate(extended_tof[,"POST_CBR_length"], list(extended_tof[,"Sequencing_ID_VLP"]), FUN = sum)
colnames(Length_ext_viruses) <- c("Sequencing_ID_VLP", "Length_ext_viruses")

VLP_metadata <- reduce(list(VLP_metadata, Length_ext_viruses), full_join, by="Sequencing_ID_VLP")

N_post_COB_viruses <- data.frame(table(extended_tof$Sequencing_ID_VLP))
colnames(N_post_COB_viruses) <- c("Sequencing_ID_VLP", "N_post_COB_viruses")

VLP_metadata <- reduce(list(VLP_metadata, N_post_COB_viruses), full_join, by="Sequencing_ID_VLP")

# N pruned contigs ps

N_pruned_ps <- as.data.frame(table(extended_tof[extended_tof$CHV_pruned=="Yes",c("Sequencing_ID_VLP", "CHV_pruned")]))
colnames(N_pruned_ps)[grep('Freq', colnames(N_pruned_ps))] <- "N_pruned_contigs"

VLP_metadata <- reduce(list(VLP_metadata, N_pruned_ps[,c("Sequencing_ID_VLP", "N_pruned_contigs")]), full_join, by="Sequencing_ID_VLP")
VLP_metadata[is.na(VLP_metadata$N_pruned_contigs),]$N_pruned_contigs <- 0

Length_pru_viruses <- aggregate(extended_tof[,"POST_CHV_length"], list(extended_tof[,"Sequencing_ID_VLP"]), FUN = sum)
colnames(Length_pru_viruses) <- c("Sequencing_ID_VLP", "Length_pru_viruses")

VLP_metadata <- reduce(list(VLP_metadata, Length_pru_viruses), full_join, by="Sequencing_ID_VLP")

# Is there enrichment of extended contigs among pruned contigs?
VLP_metadata$delta_ext <- VLP_metadata$ext_vir_perc_align - VLP_metadata$own_vir_perc_align
VLP_metadata$delta_pru <- VLP_metadata$ext_vir_prun_perc_align - VLP_metadata$ext_vir_perc_align

map_delta_corr <- cor.test(VLP_metadata$delta_pru, VLP_metadata$delta_ext, method = "spearman")
map_delta_corr$p.value

N_E1_P1 <- length(extended_tof[grep('E1_P1', extended_tof$New_CID),]$New_CID) # 441
N_E1_P0 <- length(extended_tof[grep('E1', extended_tof$New_CID),]$New_CID) - N_E1_P1 # 10,637
N_E0_P1 <- length(extended_tof[grep('E0_P1', extended_tof$New_CID),]$New_CID)
N_E0_P0 <- length(extended_tof[grep('E0_P0', extended_tof$New_CID),]$New_CID)

length_delta_corr <- cor.test(extended_tof[grep('E1_P1', extended_tof$New_CID), "POST_CBR_length"]  - extended_tof[grep('E1_P1', extended_tof$New_CID), "Original_length"],
                              extended_tof[grep('E1_P1', extended_tof$New_CID), "POST_CHV_length"] -  extended_tof[grep('E1_P1', extended_tof$New_CID), "POST_CBR_length"],
                              method = "spearman")
length_delta_corr$p.value

ext_enrich <- matrix(c(N_E1_P1, N_E0_P1, N_E1_P0, N_E0_P0), nrow=2, dimnames = list(c('Extended', 'Not_extended'),
                                                                                    c('Pruned', 'Not_pruned')))

chisq.test(ext_enrich)$p.value

# Is there a correlation between read alignment rate to VIR_DB and the ViromeQC enrichment score?

# There are a few very high (artificially) total enrichment scores, so we exclude them by introducing a cut-off
cor_maping_vir_sc <- cor.test(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,]$w_neg_der95, 
         log10(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,]$sc.enrichmnet.score), 
         method = "spearman")

cor_maping_vir_tot <- cor.test(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,]$w_neg_der95, 
                              log10(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,]$total.enrichmnet.score), 
                              method = "spearman")

cor_maping_vir_ssu <- cor.test(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,]$w_neg_der95, 
                               log10(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,]$SSU.rRNA.alignment.rate), 
                               method = "spearman")


##### plots:

# What is the distribution of postkneaddata reads in mothers and infants?
round(median(na.rm = T, N_postkn_reads$N_postkn_reads/1000000), 1)
round(sd(na.rm = T, N_postkn_reads$N_postkn_reads/1000000), 1)
pdf('05.PLOTS/03.GENERAL_STATS/N_postkn_reads.pdf', width=16/2.54, height=10/2.54)
ggplot() + 
  geom_histogram(data = N_postkn_reads[N_postkn_reads$Type=="M",], aes(N_postkn_reads/1000000, fill = "M"),  alpha = 0.2) + 
  geom_histogram(data = N_postkn_reads[N_postkn_reads$Type=="K",], aes(N_postkn_reads/1000000, fill = "K"), alpha = 0.2) +
  geom_density(data=N_postkn_reads, aes(x=N_postkn_reads/1000000, y = (after_stat(count)/sum(after_stat(count)))*19000,color=Type), alpha=0.2) + 
  labs(x="Clean reads, M", y="N samples", fill="Type") +
  geom_vline(aes(xintercept=round(median(N_postkn_reads[N_postkn_reads$Type=="M",]$N_postkn_reads/1000000, na.rm = T),1), color="M"),linewidth=1, linetype="dashed") + 
  geom_vline(aes(xintercept=round(median(N_postkn_reads[N_postkn_reads$Type=="K",]$N_postkn_reads/1000000, na.rm = T),1), color="K"),linewidth=1, linetype="dashed") +
  geom_vline(aes(xintercept=round(35489695/1000000,1), color="NC"), linewidth=1,linetype="dotted") +
  geom_vline(aes(xintercept=round(22196016/1000000,1), color="NC"), linewidth=1,linetype="dotted") +
  scale_fill_manual(breaks=c('M', 'K'), values=c("#00235B", "#E21818")) +
  scale_color_manual(breaks=c('M', 'K', 'NC'), values=c("#00235B", "#E21818", "black")) +
  my_theme +
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5))
dev.off()


# What is the distribution of clean reads in mothers and infants?
round(median(na.rm = T, ViromeQC$Reads_HQ/1000000), 1)
round(sd(na.rm = T, ViromeQC$Reads_HQ/1000000), 1)
pdf('05.PLOTS/03.GENERAL_STATS/N_clean_reads.pdf', width=16/2.54, height=10/2.54)
ggplot() + 
  geom_histogram(data = ViromeQC[ViromeQC$Type=="M",], aes(ViromeQC[ViromeQC$Type=="M",]$Reads_HQ/1000000, fill = "M"),  alpha = 0.2) + 
  geom_histogram(data = ViromeQC[ViromeQC$Type=="K",], aes(ViromeQC[ViromeQC$Type=="K",]$Reads_HQ/1000000, fill = "K"), alpha = 0.2) +
  geom_density(data=ViromeQC, aes(x=ViromeQC$Reads_HQ/1000000, y = (after_stat(count)/sum(after_stat(count)))*19000,color=Type), alpha=0.2) + 
  labs(x="Clean reads, M", y="N samples", fill="Type") +
  geom_vline(aes(xintercept=round(median(ViromeQC[ViromeQC$Type=="M",]$Reads_HQ/1000000, na.rm = T),1), color="M"),linewidth=1, linetype="dashed") + 
  geom_vline(aes(xintercept=round(median(ViromeQC[ViromeQC$Type=="K",]$Reads_HQ/1000000, na.rm = T),1), color="K"),linewidth=1, linetype="dashed") +
  geom_vline(aes(xintercept=round(19098379/1000000,1), color="NC"), linewidth=1,linetype="dotted") +
  geom_vline(aes(xintercept=round(16196692/1000000,1), color="NC"), linewidth=1,linetype="dotted") +
  scale_fill_manual(breaks=c('M', 'K'), values=c("#00235B", "#E21818")) +
  scale_color_manual(breaks=c('M', 'K', 'NC'), values=c("#00235B", "#E21818", "black")) +
  my_theme +
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5))
dev.off()


# What is the distribution of duplicated reads in mothers and infants?
round(median(na.rm = T, N_duplicated$N_duplicates/1000000), 1)
round(sd(na.rm = T, N_duplicated$N_duplicates/1000000), 1)
pdf('05.PLOTS/03.GENERAL_STATS/N_duplicated.pdf', width=16/2.54, height=10/2.54)
ggplot() + 
  geom_histogram(data = N_duplicated[N_duplicated$Type=="M",], aes(N_duplicated[N_duplicated$Type=="M",]$N_duplicates/1000000, fill = "M"),  alpha = 0.2) + 
  geom_histogram(data = N_duplicated[N_duplicated$Type=="K",], aes(N_duplicated[N_duplicated$Type=="K",]$N_duplicates/1000000, fill = "K"), alpha = 0.2) +
  geom_density(data=N_duplicated, aes(x=N_duplicated$N_duplicates/1000000, y = (after_stat(count)/sum(after_stat(count)))*19000,color=Type), alpha=0.2) + 
  labs(x="Duplicated reads, M", y="N samples", fill="Type") +
  geom_vline(aes(xintercept=round(median(N_duplicated[N_duplicated$Type=="M",]$N_duplicates/1000000, na.rm = T),1), color="M"),linewidth=1, linetype="dashed") + 
  geom_vline(aes(xintercept=round(median(N_duplicated[N_duplicated$Type=="K",]$N_duplicates/1000000, na.rm = T),1), color="K"),linewidth=1, linetype="dashed") +
  geom_vline(aes(xintercept=round(16391316/1000000,1), color="NC"), linewidth=1,linetype="dotted") +
  geom_vline(aes(xintercept=round(5999324/1000000,1), color="NC"), linewidth=1,linetype="dotted") +
  scale_fill_manual(breaks=c('M', 'K'), values=c("#00235B", "#E21818")) +
  scale_color_manual(breaks=c('M', 'K', 'NC'), values=c("#00235B", "#E21818", "black")) +
  my_theme +
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5))
dev.off()


# What is the distribution of human reads in mothers and infants?
round(median(na.rm = T, N_human$N_human_reads/1000000), 1)
round(sd(na.rm = T, N_human$N_human_reads/1000000), 1)
pdf('05.PLOTS/03.GENERAL_STATS/N_human_reads.pdf', width=16/2.54, height=10/2.54)
ggplot() + 
  geom_histogram(data = N_human[N_human$Type=="M",], aes(N_human[N_human$Type=="M",]$N_human_reads/1000000, fill = "M"),  alpha = 0.2) + 
  geom_histogram(data = N_human[N_human$Type=="K",], aes(N_human[N_human$Type=="K",]$N_human_reads/1000000, fill = "K"), alpha = 0.2) +
  geom_density(data=N_human, aes(x=N_human$N_human_reads/1000000, y = (after_stat(count)/sum(after_stat(count)))*19000,color=Type), alpha=0.2) + 
  labs(x="Human reads, M", y="N samples", fill="Type") +
  geom_vline(aes(xintercept=round(median(N_human[N_human$Type=="M",]$N_human_reads/1000000, na.rm = T),1), color="M"),linewidth=1, linetype="dashed") + 
  geom_vline(aes(xintercept=round(median(N_human[N_human$Type=="K",]$N_human_reads/1000000, na.rm = T),1), color="K"),linewidth=1, linetype="dashed") +
  geom_vline(aes(xintercept=round(3629627/1000000,1), color="NC"), linewidth=1,linetype="dotted") +
  geom_vline(aes(xintercept=round(12219624/1000000,1), color="NC"), linewidth=1,linetype="dotted") +
  scale_fill_manual(breaks=c('M', 'K'), values=c("#00235B", "#E21818")) +
  scale_color_manual(breaks=c('M', 'K', 'NC'), values=c("#00235B", "#E21818", "black")) +
  my_theme +
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5))
dev.off()


# What is the distribution of the number of assembled contigs in mothers and infants?
round(median(na.rm = T, contigs_stat$contigs.....0.bp./1000), 1)
round(sd(na.rm = T, contigs_stat$contigs.....0.bp./1000), 1)
pdf('05.PLOTS/03.GENERAL_STATS/N_contigs.pdf', width=16/2.54, height=10/2.54)
ggplot() + 
  geom_histogram(data = contigs_stat[contigs_stat$Type=="M",], aes(contigs.....0.bp./1000, fill = "M"),  alpha = 0.2) + 
  geom_histogram(data = contigs_stat[contigs_stat$Type=="K",], aes(contigs.....0.bp./1000, fill = "K"), alpha = 0.2) +
  geom_density(data=contigs_stat, aes(x=contigs.....0.bp./1000, y = (after_stat(count)/sum(after_stat(count)))*20000,color=Type), alpha=0.2) + 
  labs(x="N contigs, thousands", y="N samples", fill="Type") +
  geom_vline(aes(xintercept=round(median(contigs_stat[contigs_stat$Type=="M",]$contigs.....0.bp., na.rm = T)/1000,0), color="M"), linetype="dashed") + 
  geom_vline(aes(xintercept=round(median(contigs_stat[contigs_stat$Type=="K",]$contigs.....0.bp., na.rm = T)/1000,0), color="K"), linetype="dashed") +
  geom_vline(aes(xintercept=round(343514/1000,0), color="NC"), linetype="dotted") +
  geom_vline(aes(xintercept=round(337449/1000,0), color="NC"), linetype="dotted") +
  scale_fill_manual(breaks=c('M', 'K'), values=c("#00235B", "#E21818")) +
  scale_color_manual(breaks=c('M', 'K', 'NC'), values=c("#00235B", "#E21818", "black")) +
  my_theme +
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5))
dev.off()

# quality of assembly additional: % of reads mapped to all assembled contigs
round(median(na.rm = T, VLP_metadata$all_perc_align), 1)

pdf('05.PLOTS/03.GENERAL_STATS/Read_alignment_to_all_contigs.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=all_perc_align, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="% of reads aligned to all contigs", x="") +
  my_theme2 +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()


# What is the distribution of the number of assembled contigs in mothers and infants?
round(median(na.rm = T, contigs_stat$contigs.....1000.bp./1000), 1)
round(sd(na.rm = T, contigs_stat$contigs.....1000.bp./1000), 1)
pdf('05.PLOTS/03.GENERAL_STATS/N_contigs_large.pdf', width=16/2.54, height=10/2.54)
ggplot() + 
  geom_histogram(data = contigs_stat[contigs_stat$Type=="M",], aes(contigs.....1000.bp./1000, fill = "M"),  alpha = 0.2) + 
  geom_histogram(data = contigs_stat[contigs_stat$Type=="K",], aes(contigs.....1000.bp./1000, fill = "K"), alpha = 0.2) +
  geom_density(data=contigs_stat, aes(x=contigs.....1000.bp./1000, y = (after_stat(count)/sum(after_stat(count)))*20000,color=Type), alpha=0.2) + 
  labs(x="N contigs > 1 kbp, thousands", y="N samples", fill="Type") +
  geom_vline(aes(xintercept=round(median(contigs_stat[contigs_stat$Type=="M",]$contigs.....1000.bp., na.rm = T)/1000,0), color="M"), linetype="dashed") + 
  geom_vline(aes(xintercept=round(median(contigs_stat[contigs_stat$Type=="K",]$contigs.....1000.bp., na.rm = T)/1000,0), color="K"), linetype="dashed") +
  geom_vline(aes(xintercept=round(24798/1000,0), color="NC"), linetype="dotted") +
  geom_vline(aes(xintercept=round(29221/1000,0), color="NC"), linetype="dotted") +
  scale_fill_manual(breaks=c('M', 'K'), values=c("#00235B", "#E21818")) +
  scale_color_manual(breaks=c('M', 'K', 'NC'), values=c("#00235B", "#E21818", "black")) +
  my_theme +
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5))
dev.off()

# % of reads mapped to contigs > 1 kbp
pdf('05.PLOTS/03.GENERAL_STATS/Read_alignment_to_own_contigs_1kbp.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=`1kbp_perc_align`, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="% of reads aligned to own contigs > 1kbp", x="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

median(VLP_metadata$own_vir_perc_align)
# % of reads mapped to virus contigs
pdf('05.PLOTS/03.GENERAL_STATS/Read_alignment_to_own_vir_contigs.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=`own_vir_perc_align`, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="% of reads aligned to own virus contigs", x="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the distribution of the virus enrichment score in Chiliadal data?
median(VLP_metadata$total.enrichmnet.score)
sd(VLP_metadata$total.enrichmnet.score)
ViromeQC_Chili <- melt(VLP_metadata[ , c("Sequencing_ID_VLP", "Type_VLP", 
                                          "SSU.rRNA.alignment.rate", "LSU.rRNA.alignment.rate", 
                                          "Bacterial_Markers.alignment.rate", "total.enrichmnet.score") ])
pdf('05.PLOTS/03.GENERAL_STATS/ViromeQC_stat.pdf', width=12/2.54, height=14/2.54)
ggplot(ViromeQC_Chili[!is.na(ViromeQC_Chili$Type_VLP),], aes(x=variable, y=value, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  facet_wrap(~variable, scales='free') + 
  labs(fill='Type', color='Type') +
  my_theme3 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

median(VLP_metadata$sc.enrichmnet.score)
median(VLP_metadata[VLP_metadata$Type_VLP=="M",]$sc.enrichmnet.score, na.rm = T)
median(VLP_metadata[VLP_metadata$Type_VLP=="K",]$sc.enrichmnet.score, na.rm = T)
sd(VLP_metadata$sc.enrichmnet.score)
ViromeQC_Chili_sc <- melt(VLP_metadata[ , c("Sequencing_ID_VLP", "Type_VLP", 
                                            "total.enrichmnet.score", "sc.enrichmnet.score") ])
# What is the distribution of enrichment score based on single copy bacterial markers only?
pdf('05.PLOTS/03.GENERAL_STATS/ViromeQC_stat_sc.pdf', width=12/2.54, height=14/2.54)
ggplot(ViromeQC_Chili_sc[!is.na(ViromeQC_Chili_sc$Type_VLP),], aes(x=variable, y=value, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  facet_wrap(~variable, scales='free') + 
  labs(fill='Type', color='Type') +
  ylim(0, 150) +
  my_theme3 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the distribution in CAV?
pdf('05.PLOTS/03.GENERAL_STATS/ViromeQC_CAV_pilot.pdf', width=12/2.54, height=14/2.54)
ggplot(ViromeQC_pilot[ViromeQC_pilot$Isolation_protocol=='VLP',], aes(x=variable, y=value, fill=preseq_prep)) + 
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  geom_sina(aes(color=preseq_prep), maxwidth=0.5, alpha=0.7) +
  facet_wrap(~variable, scales='free') + 
  labs(fill='Type', color='Type') +
  my_theme3 + 
  ggtitle("Pilot VLP") +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the distribution in CAV?
pdf('05.PLOTS/03.GENERAL_STATS/ViromeQC_CAV_pilot_MGS.pdf', width=12/2.54, height=14/2.54)
ggplot(ViromeQC_pilot[ViromeQC_pilot$Isolation_protocol=='MGS',], aes(x=variable, y=value, fill=preseq_prep)) + 
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  geom_sina(aes(color=preseq_prep), maxwidth=0.5, alpha=0.7) +
  facet_wrap(~variable, scales='free') + 
  labs(fill='Type', color='Type') +
  my_theme3 + 
  ggtitle("Pilot MGS") +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# TABLE OF ORIGIN: how many viruses are discovered by every tool?
summary_per_tool <- as.data.frame(matrix(NA, nrow=5, ncol=2))
summary_per_tool$V1 <- colnames(table_of_origin[,c(2:6)])
for (i in summary_per_tool$V1) {
  summary_per_tool[summary_per_tool$V1==i,]$V2 <- sum(table_of_origin[,i])
}

listInput <- list(CenoteTaker3=table_of_origin[table_of_origin$CenoteTaker3==1,]$V1, 
                  geNomad=table_of_origin[table_of_origin$geNomad==1,]$V1,
                  VIBRANT=table_of_origin[table_of_origin$VIBRANT==1,]$V1,
                  DeepVirFinder=table_of_origin[table_of_origin$DeepVirFinder==1,]$V1,
                  VirSorter2=table_of_origin[table_of_origin$VirSorter2==1,]$V1)

pdf('05.PLOTS/03.GENERAL_STATS/VD_redundant_tools_overlap.pdf', width=18/2.54, height=10/2.54)
upset(fromList(listInput), order.by = "freq", sets.bar.color = "#C00000", 
      number.angles = 20,
      sets.x.label = "N detected virus contigs", scale.sets = "identity",
      text.scale = c(1, 1, 1, 1, 1, 0.75))
dev.off()


# What is the distribution of the number of discovered viruses?
round(median(na.rm=T, VLP_metadata[VLP_metadata$Type_VLP=='M',]$N_disc_viruses))
round(median(na.rm=T, VLP_metadata[VLP_metadata$Type_VLP=='K',]$N_disc_viruses))
pdf('05.PLOTS/03.GENERAL_STATS/N_discovered_viruses.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=N_disc_viruses, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="N virus genomes and fragments", x="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# virusness of samples: % of virus contigs from contigs > 1kbp
round(median(na.rm = T, VLP_metadata$perc_viral_contigs), 1)
round(median(na.rm=T, VLP_metadata[VLP_metadata$Type_VLP=='M',]$perc_viral_contigs), 1)
round(median(na.rm=T, VLP_metadata[VLP_metadata$Type_VLP=='K',]$perc_viral_contigs), 1)
round(sd(na.rm = T, VLP_metadata$perc_viral_contigs), 1)

pdf('05.PLOTS/03.GENERAL_STATS/Contigs_virusness.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=perc_viral_contigs, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="% of virus contigs from all contigs > 1kbp", x="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# virusness of samples: % of virus total length from total length > 1kbp 

round(median(na.rm = T, VLP_metadata$perc_viral_length), 1)
round(median(na.rm=T, VLP_metadata[VLP_metadata$Type_VLP=='M',]$perc_viral_length), 1)
round(median(na.rm=T, VLP_metadata[VLP_metadata$Type_VLP=='K',]$perc_viral_length), 1)
round(sd(na.rm = T, VLP_metadata$perc_viral_length), 1)

pdf('05.PLOTS/03.GENERAL_STATS/Contigs_length_virusness.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=perc_viral_length, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="% of virus contigs length from all contigs > 1kbp", x="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# does virusness correlate with viromeqc enrichment score?
plot(VLP_metadata$N_disc_viruses, VLP_metadata$total.enrichmnet.score)

summary(lm(data = VLP_metadata, N_disc_viruses ~ sc.enrichmnet.score))

# how many contigs ps got extended?
median(VLP_metadata$N_extended_contigs/VLP_metadata$N_disc_viruses*100)
pdf('05.PLOTS/03.GENERAL_STATS/N_extended_contigs_ps.pdf', width=8/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=N_extended_contigs, fill=Type_VLP)) +
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.color = NA, width=0.5, alpha=0.5) + 
  labs(fill='Type', color='Type', y="N of COBRA-extended virus contigs", x="") + 
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

pdf('05.PLOTS/03.GENERAL_STATS/Perc_extended_contigs_ps.pdf', width=8/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=N_extended_contigs/N_disc_viruses*100, fill=Type_VLP)) +
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.color = NA, width=0.5, alpha=0.5) + 
  labs(fill='Type', color='Type', y="% of extended contigs from total virus contigs", x="") + 
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the increase in the total length of contigs ps?
median(VLP_metadata$Length_ext_viruses - VLP_metadata$Length_disc_viruses)
pdf('05.PLOTS/03.GENERAL_STATS/Total_length_increase_after_COBRA.pdf', width=8/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=(Length_ext_viruses - Length_disc_viruses), fill=Type_VLP)) +
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.color = NA, width=0.5, alpha=0.5) + 
  labs(fill='Type', color='Type', y="Total virus length increase after extension", x="") + 
  scale_y_log10() +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

pdf('05.PLOTS/03.GENERAL_STATS/perc_Total_length_increase_after_COBRA.pdf', width=8/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=(Length_ext_viruses - Length_disc_viruses)/Length_ext_viruses*100, fill=Type_VLP)) +
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.color = NA, width=0.5, alpha=0.5) + 
  labs(fill='Type', color='Type', y="% total virus length increase after extension", x="") + 
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the number of pruned contigs ps?
median(VLP_metadata$N_pruned_contigs/VLP_metadata$N_post_COB_viruses*100)
pdf('05.PLOTS/03.GENERAL_STATS/N_pruned_contigs_ps.pdf', width=8/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=N_pruned_contigs, fill=Type_VLP)) +
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.color = NA, width=0.5, alpha=0.5) + 
  labs(fill='Type', color='Type', y="N of CheckV-pruned virus contigs", x="") + 
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

pdf('05.PLOTS/03.GENERAL_STATS/Perc_pruned_contigs_ps.pdf', width=8/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=N_pruned_contigs/N_post_COB_viruses*100, fill=Type_VLP)) +
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.color = NA, width=0.5, alpha=0.5) + 
  labs(fill='Type', color='Type', y="% of pruned contigs from total virus contigs", x="") + 
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the decrease in the total length of contigs ps?
median(VLP_metadata$Length_ext_viruses - VLP_metadata$Length_pru_viruses)
pdf('05.PLOTS/03.GENERAL_STATS/Total_length_decrease_after_CHV.pdf', width=8/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=(Length_ext_viruses - Length_pru_viruses), fill=Type_VLP)) +
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.color = NA, width=0.5, alpha=0.5) + 
  labs(fill='Type', color='Type', y="Total virus length decrease after pruning", x="") + 
  scale_y_log10() +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

pdf('05.PLOTS/03.GENERAL_STATS/perc_Total_length_decrease_after_CHV.pdf', width=8/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=(Length_ext_viruses - Length_pru_viruses)/Length_ext_viruses*100, fill=Type_VLP)) +
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.color = NA, width=0.5, alpha=0.5) + 
  labs(fill='Type', color='Type', y="% total virus length decrease after pruning", x="") + 
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()


# What is the length distribution of virus contigs?
pdf('05.PLOTS/03.GENERAL_STATS/Virus_contigs_extended_pruned_quality_distribution.pdf', width=16/2.54, height=12/2.54)
ggplot() + 
  geom_histogram(data = extended_tof, aes(POST_CHV_length, color="All", fill="All"), alpha = 0.2, bins=60) + 
  geom_histogram(data = extended_tof[extended_tof$checkv_quality=="Not-determined",], aes(POST_CHV_length, color="Not-determined", fill="Not-determined"), alpha = 0.2, bins=60) +
  geom_histogram(data = extended_tof[extended_tof$checkv_quality=="Low-quality",], aes(POST_CHV_length, color="Low-quality", fill="Low-quality"), alpha = 0.2, bins=60) +
  geom_histogram(data = extended_tof[extended_tof$checkv_quality=="Medium-quality",], aes(POST_CHV_length, color="Medium-quality", fill="Medium-quality"), alpha = 0.2, bins=60) +
  geom_histogram(data = extended_tof[extended_tof$checkv_quality=="High-quality",], aes(POST_CHV_length, color="High-quality", fill="High-quality"), alpha = 0.2, bins=60) +
  geom_histogram(data = extended_tof[extended_tof$checkv_quality=="Complete",], aes(POST_CHV_length, color="Complete", fill="Complete"), alpha = 0.2, bins=60) +
  labs(x="Virus contig length, bp", y="Log10 N virus contigs or fragments", fill="Genome Quality", color="Genome Quality") +
  scale_color_manual(breaks=c("All", "Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"), 
                     values=c("#adadad", "#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")) + 
  scale_fill_manual(breaks=c("All", "Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"), 
                    values=c("#DDDDDD", "#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")) +
  scale_x_log10(breaks=c(100, 1000, 10000, 100000, 1000000), labels=c("100", "1,000", "10,000", "100,000", "1,000,000")) +
  scale_y_log10() +
  my_theme +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()
  
# How many complete, high-quality etc virus genomes do we have in the redundant DB?
quality_summary <- as.data.frame(table(extended_tof$checkv_quality))
quality_summary$Var1 <- factor(quality_summary$Var1, 
                               levels = c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"),
                               ordered = T)

qual_lab <- c("18,937", "27,835", "1,306,382", "33,204",  "1,116,671")

pdf('05.PLOTS/03.GENERAL_STATS/Virus_contigs_extended_pruned_quality.pdf', width=10/2.54, height=10/2.54)
ggplot(quality_summary, aes(Var1, Freq, fill=Var1)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label=qual_lab), vjust=-0.25, size=3) +
  labs(fill='CheckV quality', y="N virus contigs", x="") +
  my_theme +
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

extended_tof$checkv_quality <- factor(extended_tof$checkv_quality, 
                               levels = c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"),
                               ordered = T)

# % of reads mapped to extended virus contigs
median(VLP_metadata$own_vir_perc_align)
median(VLP_metadata$ext_vir_perc_align - VLP_metadata$own_vir_perc_align)
summary(VLP_metadata$ext_vir_perc_align - VLP_metadata$own_vir_perc_align)
pdf('05.PLOTS/03.GENERAL_STATS/Read_alignment_to_own_ext_vir_contigs.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=`ext_vir_perc_align`, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="% of reads aligned to extended virus contigs", x="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

pdf('05.PLOTS/03.GENERAL_STATS/Increase_in_alignment_after_COBRA.pdf', width=12/2.54, height=4/2.54)
ggplot(data=NULL, aes(x=(VLP_metadata$ext_vir_perc_align - VLP_metadata$own_vir_perc_align), y='')) + 
  geom_sina(maxwidth=0.5, alpha=0.7, color="#3B528BFF", maxwidth=0.5) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5, fill="#3B528BFF", width=0.5) +
  labs(x="% increase in alignment", y="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# % of reads mapped to extended pruned virus contigs
median(VLP_metadata$ext_vir_prun_perc_align)
pdf('05.PLOTS/03.GENERAL_STATS/Read_alignment_to_own_ext_prun_vir_contigs.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=`ext_vir_prun_perc_align`, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="% of reads aligned to extended pruned virus contigs", x="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()


# % of reads mapped to extended pruned virus contigs
median(VLP_metadata$ext_vir_perc_align - VLP_metadata$ext_vir_prun_perc_align)
pdf('05.PLOTS/03.GENERAL_STATS/Decrease_in_alignment_after_pruning.pdf', width=12/2.54, height=4/2.54)
ggplot(data=NULL, aes(x=(VLP_metadata$ext_vir_perc_align - VLP_metadata$ext_vir_prun_perc_align), y='')) + 
  geom_sina(maxwidth=0.5, alpha=0.7, color="#3B528BFF", maxwidth=0.5) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5, fill="#3B528BFF", width=0.5) +
  labs(x="% decrease in alignment (after pruning)", y="") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# Are extended contigs pruned more often compared to unextended?
pdf('05.PLOTS/03.GENERAL_STATS/Spearman_correlation_delta_mapping_extension_vs_pruning.pdf', width=10/2.54, height=12/2.54)
ggplot(VLP_metadata, aes(delta_ext, delta_pru, color=Type_VLP)) +
  geom_point(alpha=0.8) + 
  geom_smooth(inherit.aes = F, aes(delta_ext, delta_pru), method="lm") + 
  labs(color="Type", x="% increase in read mapping after COBRA extension", y="% decrease in read mapping after CheckV pruning") +
  annotate("text", x=12, y=-10, label="rho = -0.32\np-value=1.0e-27") +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

# What is the pruning distribution of delta
pdf('05.PLOTS/03.GENERAL_STATS/Spearman_correlation_delta_length_extension_vs_pruning.pdf', width=14/2.54, height=12/2.54)
ggplot(extended_tof[grep('E1_P1', extended_tof$New_CID),], aes(x=POST_CBR_length - Original_length,
                                                               y=POST_CHV_length - POST_CBR_length)) +
  geom_point(alpha=0.8) + 
  geom_smooth(method="lm") + 
  labs(x="Length difference between extended and original contigs", y="Length difference between pruned and extended contigs") +
  annotate("text", x=1e+06, y=-5e+05, label="rho = -0.49\np-value=2.0e-27") +
  my_theme2 + 
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

pdf('05.PLOTS/03.GENERAL_STATS/Mosaic_plot_ext_prun.pdf', width=16/2.54, height=18/2.54)
mosaicplot(ext_enrich, main = "Enrichment in pruning of extended contigs", color = T, sub = "p-value=4.9e-239")
dev.off()

# % of reads mapped to virus contigs dereplicated at 95% ANI over 85% of contig length
median(VLP_metadata$w_neg_der95)
pdf('05.PLOTS/03.GENERAL_STATS/Read_alignment_to_w_neg_der95.pdf', width=12/2.54, height=14/2.54)
ggplot(VLP_metadata[!is.na(VLP_metadata$Type_VLP),], aes(x=Type_VLP, y=`w_neg_der95`, fill=Type_VLP)) + 
  geom_sina(aes(color=Type_VLP), maxwidth=0.5, alpha=0.7) +
  geom_boxplot(outlier.colour = NA, width=0.5, alpha=0.5) +
  labs(fill='Type', color='Type', y="% of reads aligned to vOTUs (NGCTRL contigs+)", x="") +
  ylim(0,100) +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

median(VLP_metadata$w_neg_der95 - VLP_metadata$ext_vir_prun_perc_align)

# correlation between enrichment score and the number of reads aligning to dereplicated virus contigs
pdf('05.PLOTS/03.GENERAL_STATS/Spearman_corr_ViromeQC_total_perc_align_w_neg_der95.pdf', width=12/2.54, height=12/2.54)
ggplot(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,], aes(x=total.enrichmnet.score, y=w_neg_der95)) + 
  geom_point(aes(color=Type_VLP), alpha=0.8) + 
  geom_smooth(inherit.aes = F, aes(x=total.enrichmnet.score, y=w_neg_der95), method="lm") + 
  labs(color='Type', y="% of reads aligned to vOTUs (NGCTRL contigs+)", x="Total enrichment score (ViromeQC)") +
  annotate("text", x=5, y=60, label="rho = 0.12\np-value=1.5e-4") +
  scale_x_log10() +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

pdf('05.PLOTS/03.GENERAL_STATS/Spearman_corr_ViromeQC_sc_perc_align_w_neg_der95.pdf', width=12/2.54, height=12/2.54)
ggplot(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,], aes(x=sc.enrichmnet.score, y=w_neg_der95)) + 
  geom_point(aes(color=Type_VLP), alpha=0.8) + 
  geom_smooth(inherit.aes = F, aes(x=sc.enrichmnet.score, y=w_neg_der95), method="lm") + 
  labs(color='Type', y="% of reads aligned to vOTUs (NGCTRL contigs+)", x="Single-copy virus enrichment score (ViromeQC)") +
  annotate("text", x=30, y=60, label="rho = 0.64\np-value=5.5e-122") +
  scale_x_log10() +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

pdf('05.PLOTS/03.GENERAL_STATS/Spearman_corr_ViromeQC_SSU_perc_align_w_neg_der95.pdf', width=12/2.54, height=12/2.54)
ggplot(VLP_metadata[VLP_metadata$sc.enrichmnet.score<300,], aes(x=SSU.rRNA.alignment.rate, y=w_neg_der95)) + 
  geom_point(aes(color=Type_VLP), alpha=0.8) + 
  geom_smooth(inherit.aes = F, aes(x=SSU.rRNA.alignment.rate, y=w_neg_der95), method="lm") + 
  labs(color='Type', y="% of reads aligned to vOTUs (NGCTRL contigs+)", x="ALignment to SSU rRNA") +
  annotate("text", x=0.05, y=60, label="rho = -0.12\np-value=1.1e-4") +
  scale_x_log10() +
  my_theme2 + 
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1))
dev.off()

##############################
# OUTPUT
##############################
write.table(VLP_metadata, '06.CLEAN_DATA/Intermediate/VLP_metadata_temporary.txt', sep='\t', row.names = F)
