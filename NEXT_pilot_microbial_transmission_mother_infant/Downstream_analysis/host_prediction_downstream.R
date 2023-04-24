library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(ggpubr)
library(maditr)
library(tidyr)
library(utils)
library(lmerTest)

##################################### INITIAL LOADING AND DATA PREPARATION #####################################

meta_ini = read.table("C:\\Users\\Natal\\PycharmProjects\\iphop\\metadata\\VLP_metadata_with_phenos.txt", header=TRUE, sep='\t')

vlp_meta_ini = read.table("C:\\Users\\Natal\\PycharmProjects\\iphop\\metadata\\VLP_viral_contigs_metadata.txt", header=TRUE, sep='\t')
names(vlp_meta_ini)[names(vlp_meta_ini) == 'V1'] <- 'Virus'

assignments_ini = read.table("C:\\Users\\Natal\\PycharmProjects\\iphop\\final_data\\MERGED_Host_prediction_to_genus_m90.csv", header=TRUE, sep=',')
assignments_3000 <- assignments_ini %>% filter(Virus %in% vlp_meta_ini$Virus[vlp_meta_ini$length >= 3000])
assignments_3000[assignments_3000$Host.genus == 'd__Bacteria;p__Firmicutes_A;c__Clostridia;o__UMGS1883;f__;g__UMGS1540', "Host.genus"] <- 'd__Bacteria;p__Firmicutes_A;c__Clostridia;o__UMGS1883;f__UMGS1883;g__UMGS1540'
assignments_3000_fullclassif <- assignments_3000[grep("g__$", assignments_3000$Host.genus, invert = TRUE), ]

abundance_table_ini = read.table("C:\\Users\\Natal\\PycharmProjects\\iphop\\metadata\\RPKM_counts_VLP.txt", header=TRUE, sep='\t')
abundance_table_ini$Virus = row.names(abundance_table_ini)
abundance_table_ini <- abundance_table_ini %>% relocate(Virus)
rownames(abundance_table_ini) <- NULL
abundance_table_3000_fullclassif <- merge(x=abundance_table_ini, y=assignments_3000_fullclassif[, c("Virus", "Host.genus")], by="Virus", all.y = TRUE)

generalist_df <- as.data.frame(table(assignments_3000_fullclassif$Virus))
colnames(generalist_df) <- c("Virus", "genuses_infected")
generalist_df <- generalist_df %>% mutate(generalist = ifelse(genuses_infected > 1, 1, 0))

vlp_meta_ini_3000_fullclassif <- merge(x=generalist_df, y=vlp_meta_ini, by="Virus",all.x = TRUE)
abundance_table_3000_fullclassif <-  merge(x=abundance_table_3000_fullclassif, y=generalist_df[, c("Virus", "generalist")], by="Virus", all.x = TRUE)
abundance_table_3000_fullclassif <-  merge(x=abundance_table_3000_fullclassif, y=vlp_meta_ini[, c("Virus", "temperate")], by="Virus", all.x = TRUE)

##################################### TEST FOR BLAST & CRISPR ONLY HITS ##################################### 

assignments_generalist <- merge(x=generalist_df, y=assignments_3000_fullclassif, by="Virus", all.y = TRUE)
assignments_generalist <- assignments_generalist[assignments_generalist$generalist == 1, ]
assignments_generalist_RED <- assignments_generalist[grep(x=assignments_generalist$List.of.methods, pattern="blast|CRISPR"), ]

new_generalist_counts <- as.data.frame(table(as.character(assignments_generalist_RED$Virus)))
colnames(new_generalist_counts) <- c('Virus', 'genera_infected_RED')

assignments_generalist_RED <- merge(x=assignments_generalist_RED, y=new_generalist_counts, all.x= TRUE)
assignments_generalist_RED_95perc <- assignments_generalist_RED[assignments_generalist_RED$Confidence.score >= 95, ]

new_generalist_counts_95perc <- as.data.frame(table(as.character(assignments_generalist_RED_95perc$Virus)))
colnames(new_generalist_counts_95perc) <- c('Virus', 'genera_infected_RED_95perc')
assignments_generalist_RED_95perc <- merge(x=assignments_generalist_RED, y=new_generalist_counts_95perc, all.x= TRUE)


length(unique(assignments_generalist$Virus))
length(unique(assignments_generalist_RED$Virus))

#write.table(assignments_generalist, file = "C:\\Users\\Natal\\PycharmProjects\\iphop\\metadata\\Generalists_host_prediction.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
##################################### MATERNAL ABUNDANCE TABLE PREPARATION #####################################

abundance_table_3000_fullclassif_mom <- abundance_table_3000_fullclassif %>% 
  select(c('Virus', 'Host.genus', "generalist", "temperate"), meta_ini$Short_sample_ID[meta_ini$Type == 'Mother'])

abundance_table_3000_fullclassif_mom$vir_presence_cumulative <- rowSums(abundance_table_3000_fullclassif_mom[, grep("V$", colnames(abundance_table_3000_fullclassif_mom))], na.rm = TRUE)

abundance_table_3000_fullclassif_mom <- abundance_table_3000_fullclassif_mom %>% 
  mutate(vir_presence = ifelse(vir_presence_cumulative > 0, 1, 0)) %>%
  filter(vir_presence == 1) %>%
  mutate(generalist_presence_cumulative = ifelse(generalist == 1, vir_presence_cumulative, 0))

abundance_table_3000_fullclassif_mom_genus <- abundance_table_3000_fullclassif_mom %>%
  select(-"Virus")
abundance_table_3000_fullclassif_mom_genus <- aggregate(.~Host.genus, abundance_table_3000_fullclassif_mom_genus, sum)

abundance_table_3000_fullclassif_mom_genus_edited <- abundance_table_3000_fullclassif_mom_genus
abundance_table_3000_fullclassif_mom_genus_edited[, grep("V$", colnames(abundance_table_3000_fullclassif_mom_genus_edited))] <- sapply(abundance_table_3000_fullclassif_mom_genus_edited[, grep("V$", colnames(abundance_table_3000_fullclassif_mom_genus_edited))],  function(x) replace(x, x > 0, 1))
abundance_table_3000_fullclassif_mom_genus_edited$vir_prevalence <- rowSums(abundance_table_3000_fullclassif_mom_genus_edited[, grep("V$", colnames(abundance_table_3000_fullclassif_mom_genus_edited))], na.rm = TRUE)
abundance_table_3000_fullclassif_mom_genus_edited$genus <- gsub('.*g__', "", abundance_table_3000_fullclassif_mom_genus_edited$Host.genus)
abundance_table_3000_fullclassif_mom_genus_edited$phylum <- gsub('.*p__|;c__.*', "", abundance_table_3000_fullclassif_mom_genus_edited$Host.genus)
abundance_table_3000_fullclassif_mom_genus_edited$generalist_perc <- abundance_table_3000_fullclassif_mom_genus_edited$generalist / abundance_table_3000_fullclassif_mom_genus_edited$vir_presence * 100
abundance_table_3000_fullclassif_mom_genus_edited$temperate_perc <- abundance_table_3000_fullclassif_mom_genus_edited$temperate / abundance_table_3000_fullclassif_mom_genus_edited$vir_presence * 100

##################################### INFANT ABUNDANCE TABLE PREPARATION #####################################

abundance_table_3000_fullclassif_inf <- abundance_table_3000_fullclassif %>% 
  select(c('Virus', 'Host.genus', "generalist", "temperate"), meta_ini$Short_sample_ID[meta_ini$Type == 'Infant'])

abundance_table_3000_fullclassif_inf$vir_presence_cumulative <- rowSums(abundance_table_3000_fullclassif_inf[, grep("V$", colnames(abundance_table_3000_fullclassif_inf))], na.rm = TRUE)

abundance_table_3000_fullclassif_inf <- abundance_table_3000_fullclassif_inf %>% 
  mutate(vir_presence = ifelse(vir_presence_cumulative > 0, 1, 0)) %>%
  filter(vir_presence == 1) %>%
  mutate(generalist_presence_cumulative = ifelse(generalist == 1, vir_presence_cumulative, 0))

abundance_table_3000_fullclassif_inf_genus <- abundance_table_3000_fullclassif_inf %>%
  select(-"Virus")
abundance_table_3000_fullclassif_inf_genus <- aggregate(.~Host.genus, abundance_table_3000_fullclassif_inf_genus, sum)

abundance_table_3000_fullclassif_inf_genus_edited <- abundance_table_3000_fullclassif_inf_genus
abundance_table_3000_fullclassif_inf_genus_edited[, grep("V$", colnames(abundance_table_3000_fullclassif_inf_genus_edited))] <- sapply(abundance_table_3000_fullclassif_inf_genus_edited[, grep("V$", colnames(abundance_table_3000_fullclassif_inf_genus_edited))],  function(x) replace(x, x > 0, 1))
abundance_table_3000_fullclassif_inf_genus_edited$vir_prevalence <- rowSums(abundance_table_3000_fullclassif_inf_genus_edited[, grep("V$", colnames(abundance_table_3000_fullclassif_inf_genus_edited))], na.rm = TRUE)
abundance_table_3000_fullclassif_inf_genus_edited$genus <- gsub('.*g__', "", abundance_table_3000_fullclassif_inf_genus_edited$Host.genus)
abundance_table_3000_fullclassif_inf_genus_edited$phylum <- gsub('.*p__|;c__.*', "", abundance_table_3000_fullclassif_inf_genus_edited$Host.genus)
abundance_table_3000_fullclassif_inf_genus_edited$generalist_perc <- abundance_table_3000_fullclassif_inf_genus_edited$generalist / abundance_table_3000_fullclassif_inf_genus_edited$vir_presence * 100
abundance_table_3000_fullclassif_inf_genus_edited$temperate_perc <- abundance_table_3000_fullclassif_inf_genus_edited$temperate / abundance_table_3000_fullclassif_inf_genus_edited$vir_presence * 100


##################################### FIRST (PREVALENCE MOTHER/INFANT) GRAPH #####################################

## MOTHERS DATAFRAME FOR THE GRAPH
df_graph1_mom <- abundance_table_3000_fullclassif_mom_genus_edited[order(-abundance_table_3000_fullclassif_mom_genus_edited$vir_prevalence, -abundance_table_3000_fullclassif_mom_genus_edited$vir_presence), ]
df_graph1_mom <- df_graph1_mom[1:20, ]
df_graph1_mom <- df_graph1_mom[order(df_graph1_mom$vir_prevalence, df_graph1_mom$vir_presence), ]
df_graph1_mom$genus <- as.factor(df_graph1_mom$genus)
df_graph1_mom <- mutate(df_graph1_mom, genus = factor(genus, genus))

## INFANTS DATAFRAME FOR THE GRAPH
df_graph1_inf <- abundance_table_3000_fullclassif_inf_genus_edited[order(-abundance_table_3000_fullclassif_inf_genus_edited$vir_prevalence, -abundance_table_3000_fullclassif_inf_genus_edited$vir_presence), ]
df_graph1_inf <- df_graph1_inf[1:20, ]
df_graph1_inf <- df_graph1_inf[order(df_graph1_inf$vir_prevalence, df_graph1_inf$vir_presence), ]
df_graph1_inf$genus <- as.factor(df_graph1_inf$genus)
df_graph1_inf <- mutate(df_graph1_inf, genus = factor(genus, genus))


## MOTHER PANNEL GRAPHS
p1_mom1 <- ggplot(df_graph1_mom, aes(x=genus, y=vir_prevalence, fill=phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none")+
  labs(x="Maternal samples \n Host genus", fill = "Phylum")+
  expand_limits(y = 150) +
  scale_y_continuous(breaks=seq(from=0,to=150,by=50))+
  scale_fill_manual(values=c('#FF9C0E', '#1F77B4', '#ED3419'))

p1_mom2 <- ggplot(df_graph1_mom, aes(x=genus, y=vir_presence, fill=phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")+
  labs(fill = "Phylum")+
  expand_limits(y = 10000) +
  scale_y_log10()+
  scale_fill_manual(values=c('#FF9C0E', '#1F77B4', '#ED3419'))

p1_mom3 <- ggplot(df_graph1_mom, aes(x=genus, y=generalist_perc, fill=phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")+
  labs(y="% Generaists", x="Host genus", fill = "Phylum")+
  expand_limits(y = 100) +
  scale_y_continuous(breaks=seq(from=0,to=100,by=50),
                     labels=c("0","50","100"))+
  scale_fill_manual(values=c('#FF9C0E', '#1F77B4', '#ED3419'))

p1_mom4 <- ggplot(df_graph1_mom, aes(x=genus, y=temperate_perc, fill=phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")+
  labs(y="% Temperate", fill = "Phylum")+
  expand_limits(y = 100) +
  scale_y_continuous(breaks=seq(from=0,to=100,by=50),
                     labels=c("0","50","100"))+
  scale_fill_manual(values=c('#FF9C0E', '#1F77B4', '#ED3419'))


## INFANT PANNEL GRAPHS
p1_inf1 <- ggplot(df_graph1_inf, aes(x=genus, y=vir_prevalence, fill=phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(hjust=1))+
  labs(y="Detected in\nN samples", x="Infant samples \n Host genus", fill = "Phylum")+
  expand_limits(y = 150) +
  scale_y_continuous(breaks=seq(from=0,to=150,by=50),
                     labels=c("0","50","100","150"))+
  scale_fill_manual(values=c('#FF9C0E', '#1F77B4', '#C7C7C7', '#C61A09',  '#ED3419', '#FF6242', '#B446B3'))

p1_inf2 <- ggplot(df_graph1_inf, aes(x=genus, y=vir_presence, fill=phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(hjust=1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(y="Number of viruses\n assigned (log10)", x="Host genus", fill = "Phylum")+
  expand_limits(y = 10000) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values=c('#FF9C0E', '#1F77B4', '#C7C7C7', '#C61A09',  '#ED3419', '#FF6242', '#B446B3'))

p1_inf3 <- ggplot(df_graph1_inf, aes(x=genus, y=generalist_perc, fill=phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(hjust=1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(y="% Generaists", x="Host genus", fill = "Phylum")+
  expand_limits(y = 100) +
  scale_y_continuous(breaks=seq(from=0,to=100,by=50),
                     labels=c("0","50","100"))+
  scale_fill_manual(values=c('#FF9C0E', '#1F77B4', '#C7C7C7', '#C61A09',  '#ED3419', '#FF6242', '#B446B3'))

p1_inf4 <- ggplot(df_graph1_inf, aes(x=genus, y=temperate_perc, fill=phylum)) + 
  geom_bar(stat = "identity") +
  coord_flip()+
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(hjust=1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(y="% Temperate", x="Host genus", fill = "Phylum")+
  expand_limits(y = 100) +
  scale_y_continuous(breaks=seq(from=0,to=100,by=50),
                     labels=c("0","50","100"))+
  scale_fill_manual(values=c('#FF9C0E', '#1F77B4', '#C7C7C7', '#C61A09',  '#ED3419', '#FF6242', '#B446B3'))

##COMBINING INTO ONE FINAL GRAPH
plot1 <- p1_mom1 + p1_mom2 + p1_mom3 + p1_mom4 + p1_inf1 + p1_inf2 + p1_inf3 + p1_inf4 + plot_layout(ncol = 4, guides = "collect")

##################################### SECOND (PREVALENCE MOTHER/INFANT) GRAPH #####################################

df_graph2 <- as.data.frame(table(generalist_df$genuses_infected))
colnames(df_graph2) <- c("genuses_infected", "virus_count")
df_graph2 <- df_graph2 %>% mutate(gen_spec = ifelse(genuses_infected == 1, "specialist", "generalist"))


# df_graph2$genuses_infected <- as.character(df_graph2$genuses_infected)
# df_graph2[nrow(df_graph2) + 1,] <- c("Unk", nrow(abundance_table_ini) - nrow(generalist_df), "unknown")
df_graph2$gen_spec <- as.factor(df_graph2$gen_spec)
df_graph2$gen_spec <- factor(df_graph2$gen_spec,
                             levels = levels(df_graph2$gen_spec)[c(2, 1)]) #add 3 in c(2, 1) in case of coming back to unk


plot2 <- ggplot(df_graph2, aes(x=genuses_infected, y=virus_count, fill=gen_spec)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label=virus_count), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank())+
  labs(x="Number of predicted hosts (genus level)", y="Number of viral contigs (log10)", fill = "")+#
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values=c('navy', "violetred3")) # add '#C7C7C7' in case of coming back to unk
plot2
#################################### GENERALIST METATADA PREPARATION PER SAMPLE #####################################

generalist_per_sample_calc <- abundance_table_ini %>%
  select(-"Virus")
viral_abundance <- colSums(generalist_per_sample_calc)
viral_richness_check <- colSums(generalist_per_sample_calc != 0)

generalist_per_sample_calc <-  merge(x=abundance_table_ini, y=generalist_df[, c("Virus", "generalist")], by="Virus", all.y = TRUE)
generalist_per_sample_calc <- generalist_per_sample_calc %>%
  select(-"Virus")
viral_abundance_classified <- colSums(generalist_per_sample_calc)
viral_richness_classified <- colSums(generalist_per_sample_calc != 0)
generalist_per_sample_calc <- generalist_per_sample_calc[generalist_per_sample_calc$generalist == 1, ]
generalist_abundance <- colSums(generalist_per_sample_calc)
generalist_richness <- colSums(generalist_per_sample_calc != 0)

generalist_per_sample_calc_final <- as.data.frame(cbind(viral_abundance_classified, viral_richness_classified, generalist_abundance,generalist_richness))
generalist_per_sample_calc_final$generalist_perc <- generalist_per_sample_calc_final$generalist_richness / generalist_per_sample_calc_final$viral_richness_classified * 100
generalist_per_sample_calc_final$generalist_RA <- generalist_per_sample_calc_final$generalist_abundance / generalist_per_sample_calc_final$viral_abundance_classified * 100

generalist_per_sample_calc_final <- generalist_per_sample_calc_final[!(row.names(generalist_per_sample_calc_final) %in% c("generalist")), ]
generalist_per_sample_calc_final <- cbind(generalist_per_sample_calc_final, viral_abundance, viral_richness_check)
generalist_per_sample_calc_final$viruses_classified_RA <- generalist_per_sample_calc_final$viral_abundance_classified / generalist_per_sample_calc_final$viral_abundance * 100
generalist_per_sample_calc_final$viruses_classified_perc <- generalist_per_sample_calc_final$viral_richness_classified / generalist_per_sample_calc_final$viral_richness_check * 100

rownames(meta_ini) <- meta_ini$Short_sample_ID
meta_extended <- merge(x=meta_ini , y=generalist_per_sample_calc_final, by=0)
all(meta_extended$viral_richness == meta_extended$viral_richness_check)
rownames(meta_extended) <- meta_extended$Row.names
meta_extended <- meta_extended %>%
  select(-"Row.names")

meta_extended$Type <- as.factor(meta_extended$Type)
meta_extended$Type <- factor(meta_extended$Type,
                             levels = levels(meta_extended$Type)[c(2, 1)])

meta_extended$Timepoint <- as.factor(meta_extended$Timepoint)
meta_extended$Timepoint <- factor(meta_extended$Timepoint,
                             levels = levels(meta_extended$Timepoint)[c(7, 1, 2, 4, 5, 6, 3)])

meta_extended$NEXT_ID <- as.factor(meta_extended$NEXT_ID)
meta_extended$FAM_ID <- as.factor(meta_extended$FAM_ID)

##################################### THIRD (COMPARISON OF VIRAL/GENERALIST RICHNESS/ABUNDANCE PER SAMPLE) GRAPH #####################################

## Calculating models for the graph

p3_1_wilcox <- wilcox.test(meta_extended$viral_richness[meta_extended$Type == "Infant"], 
            meta_extended$viral_richness[meta_extended$Type == "Mother"], alternative = "two.sided")
p3_1_lm <- lmer(viral_richness ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended)


p3_2_wilcox <- wilcox.test(meta_extended$generalist_richness[meta_extended$Type == "Infant"], 
            meta_extended$generalist_richness[meta_extended$Type == "Mother"], alternative = "two.sided")
p3_2_lm <- lmer(generalist_richness ~ Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended)

p3_3_wilcox <- wilcox.test(meta_extended$generalist_perc[meta_extended$Type == "Infant"], 
                           meta_extended$generalist_perc[meta_extended$Type == "Mother"], alternative = "two.sided")
p3_3_lm <- lmer(generalist_perc ~  Type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended)

p3_4_wilcox <- wilcox.test(meta_extended$generalist_RA[meta_extended$Type == "Infant"], 
                           meta_extended$generalist_RA[meta_extended$Type == "Mother"], alternative = "two.sided")
p3_4_lm <- lmer(generalist_RA ~  Type+ DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended)

p3_1 <- ggplot(meta_extended, aes(x=Type, y=viral_richness, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
    legend.position="none",
    plot.title = element_text(size=11)) +
  annotate("text", x=1.8, y=20000, label= paste("p-value=", "2.2e-16 (wilcox);", "1.58e-07 (lmer)")) +
  labs(y="Phage richness")+
  xlab("")

p3_2 <- ggplot(meta_extended, aes(x=Type, y=generalist_richness, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=11)) +
  labs(y="Generalist richness")+
  annotate("text", x=1.8, y=850, label= paste("p-value=", "5.125e-16 (wilcox);", "3.04e-10 (lmer)")) +
  xlab("")

p3_3 <- ggplot(meta_extended, aes(x=Type, y=generalist_perc, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=11)) +
  expand_limits(y = 100) +
  labs(y="% generalists")+
  annotate("text", x=1.8, y=100, label= paste("p-value=", "9.994e-06 (wilcox);", "0.017437 (lmer)")) +
  xlab("")

p3_4 <- ggplot(meta_extended, aes(x=Type, y=generalist_RA, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=11)) +
  labs(y="Generalists relative abundance")+
  annotate("text", x=1.8, y=100, label= paste("p-value=", "0.1668 (wilcox);", "0.028673 (lmer)")) +
  xlab("")

plot3 <- p3_1 + p3_2 + p3_3 + p3_4 + plot_layout(ncol = 2, guides = "collect")
plot3

##################################### FOURTH (COMPARISON OF GENERALIST COUNTS PER SAMPLE PER TIMEPIONT) GRAPH #####################################


p4_mom1_lm <- lmer(generalist_richness ~ Timepoint_continuous + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Mother',])
p4_mom2_lm <- lmer(generalist_perc ~ Timepoint_continuous + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Mother',])
p4_mom3_lm <- lmer(generalist_RA ~ Timepoint_continuous + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Mother',])


p4_inf1_lm <- lmer(generalist_richness ~ Timepoint_continuous + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Infant',])
p4_inf2_lm <- lmer(generalist_perc ~ Timepoint_continuous + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Infant',])
p4_inf3_lm <- lmer(generalist_RA ~ Timepoint_continuous + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Infant',])

data_for_graph <- melt(meta_extended, id.vars = c("Short_sample_ID", "Timepoint", "Type") ,measure.vars = c("generalist_richness", "generalist_perc", "generalist_RA"))  

new.labs <- c("Generalist richness", "Generalist %", "Generalist RA")
names(new.labs) <- c("generalist_richness", "generalist_perc", "generalist_RA")

plot4 <- ggplot(data_for_graph, aes(x=Timepoint, y=value, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  facet_grid(variable~Type, scales = "free", labeller = labeller(variable = new.labs)) +
  scale_y_continuous() +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5))+
  labs(x="", y="", title="")
plot4

##################################### FIFTH (RICHNESS/ABUNDANCE PERCENT OF CLASSIFIED VIRUSES) GRAPH #####################################
df_graph5 <- melt(meta_extended[, c("Short_sample_ID", "viruses_classified_perc", "viruses_classified_RA")], id="Short_sample_ID")

plot5 <- ggplot(df_graph5, aes(x=variable, y=value, fill=variable)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('lightblue', 'lightblue')) +
  geom_jitter(aes(color=variable), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('lightblue4', 'lightblue4')) +
  scale_x_discrete(labels=c(viruses_classified_perc = "Richness % of assigned viruses", viruses_classified_RA = "Abundance % of assigned viruses")) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=16, hjust=0.5)) +
  labs(x="", y="%", title="")
  
plot5


##################################### MATERNAL/INFANT SPECIFIC VIRUSES TABLE PREPARATION #####################################

maternal_all_viruses <- abundance_table_3000_fullclassif_mom$Virus[!duplicated(abundance_table_3000_fullclassif_mom$Virus)]
infant_all_viruses <- abundance_table_3000_fullclassif_inf$Virus[!duplicated(abundance_table_3000_fullclassif_inf$Virus)]
unspecific_viruses <- intersect(maternal_all_viruses, infant_all_viruses)
maternal_specific_viruses <- setdiff(maternal_all_viruses, infant_all_viruses)
infant_specific_viruses <- setdiff(infant_all_viruses, maternal_all_viruses)

generalist_mom_specific <- merge(x=abundance_table_ini, y=generalist_df[, c("Virus", "generalist")], by="Virus", all.y = TRUE) %>%
  filter(Virus %in% maternal_specific_viruses) %>%
  select(-"Virus")

generalist_inf_specific <- merge(x=abundance_table_ini, y=generalist_df[, c("Virus", "generalist")], by="Virus", all.y = TRUE) %>%
  filter(Virus %in% infant_specific_viruses) %>%
  select(-"Virus")

generalist_unspecific <- merge(x=abundance_table_ini, y=generalist_df[, c("Virus", "generalist")], by="Virus", all.y = TRUE) %>%
  filter(Virus %in% unspecific_viruses) %>%
  select(-"Virus")


generalist_mom_specific_RA <- colSums(generalist_mom_specific[generalist_mom_specific$generalist == 1, ]) / colSums(generalist_mom_specific) * 100
generalist_mom_specific_perc <- colSums(generalist_mom_specific[generalist_mom_specific$generalist == 1, ] != 0) / colSums(generalist_mom_specific != 0) * 100

generalist_inf_specific_RA <- colSums(generalist_inf_specific[generalist_inf_specific$generalist == 1, ]) / colSums(generalist_inf_specific) * 100
generalist_inf_specific_perc <- colSums(generalist_inf_specific[generalist_inf_specific$generalist == 1, ] != 0) / colSums(generalist_inf_specific != 0) * 100

generalist_unspecific_RA <- colSums(generalist_unspecific[generalist_unspecific$generalist == 1, ]) / colSums(generalist_unspecific) * 100
generalist_unspecific_perc <- colSums(generalist_unspecific[generalist_unspecific$generalist == 1, ] != 0) / colSums(generalist_unspecific != 0) * 100

generalist_specific_final <- as.data.frame(cbind(generalist_mom_specific_RA, generalist_mom_specific_perc, generalist_inf_specific_RA, generalist_inf_specific_perc, generalist_unspecific_RA, generalist_unspecific_perc))
generalist_specific_final[is.na(generalist_specific_final)] <- 0

generalist_specific_final$generalist_specific_RA <- generalist_specific_final$generalist_mom_specific_RA + generalist_specific_final$generalist_inf_specific_RA
generalist_specific_final$generalist_specific_perc <- generalist_specific_final$generalist_mom_specific_perc + generalist_specific_final$generalist_inf_specific_perc

generalist_specific_final <- generalist_specific_final %>%
  select(-c(generalist_mom_specific_RA, generalist_mom_specific_perc, generalist_inf_specific_RA, generalist_inf_specific_perc))

generalist_specific_final <- generalist_specific_final[!(row.names(generalist_specific_final) %in% c("generalist")), ]

meta_extended <- merge(x=meta_extended, y=generalist_specific_final, by=0)
meta_extended <- meta_extended %>%
  select(-"Row.names")

##################################### SIXTH (COMPARISON OF GENERALISTS IN MOM/INF SPECIFIC VIRUSES PER SAMPLE) GRAPH #####################################

p6_1 <- ggplot(meta_extended, aes(x=Type, y=generalist_specific_perc, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=16, hjust=0.5)) +
  annotate("text", x=2.3, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_specific_perc[meta_extended$Type == "Infant"], 
                                                                    meta_extended$generalist_specific_perc[meta_extended$Type == "Mother"], alternative = "two.sided")$p.value, digits=3))) +
  labs(x="", y="Generalist percent", title="Specific viruses")

p6_2 <- ggplot(meta_extended, aes(x=Type, y=generalist_unspecific_perc, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=16, hjust=0.5)) +
  annotate("text", x=2.3, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_unspecific_perc[meta_extended$Type == "Infant"], 
                                                                    meta_extended$generalist_unspecific_perc[meta_extended$Type == "Mother"], alternative = "two.sided")$p.value, digits=3))) +
  labs(x="", y="", title="Unspecific viruses")

p6_3 <- ggplot(meta_extended, aes(x=Type, y=generalist_specific_RA, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=11)) +
  annotate("text", x=2.3, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_specific_RA[meta_extended$Type == "Infant"], 
                                                                    meta_extended$generalist_specific_RA[meta_extended$Type == "Mother"], alternative = "two.sided")$p.value, digits=3))) +
  labs(x="", y="Generalist relative abundance")

p6_4 <- ggplot(meta_extended, aes(x=Type, y=generalist_unspecific_RA, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=11)) +
  annotate("text", x=2.3, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_unspecific_RA[meta_extended$Type == "Infant"], 
                                                                    meta_extended$generalist_unspecific_RA[meta_extended$Type == "Mother"], alternative = "two.sided")$p.value, digits=3))) +
  labs(x="", y="")

plot6 <- p6_1 + p6_2 + p6_3 + p6_4 + plot_layout(ncol = 2, guides = "collect")


##################################### SEVENTH (COMPARISON OF GENERALIST RICHNESS/ABUNDANCE/PERCENT MOTHER VS BABIES PER SAMPLE ON SAME TIMEPOINS) GRAPH #####################################

p7_1 <- ggplot(meta_extended[meta_extended$Timepoint == "M1", ], aes(x=Type, y=generalist_richness, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  expand_limits(y = 900) +
  annotate("text", x=2.2, y=900, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_richness[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M1"], 
                                                                    meta_extended$generalist_richness[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M1"], alternative = "two.sided")$p.value, digits=3))) +
  labs(x="", y="Generalist richness", title="M1")

p7_2 <- ggplot(meta_extended[meta_extended$Timepoint == "M2", ], aes(x=Type, y=generalist_richness, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  annotate("text", x=2.2, y=900, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_richness[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M2"], 
                                                                    meta_extended$generalist_richness[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M2"], alternative = "two.sided")$p.value, digits=3))) +
  expand_limits(y = 900) +
  labs(x="", y="", title="M2")

p7_3 <- ggplot(meta_extended[meta_extended$Timepoint == "M3", ], aes(x=Type, y=generalist_richness, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  annotate("text", x=2.2, y=900, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_richness[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M3"], 
                                                                    meta_extended$generalist_richness[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M3"], alternative = "two.sided")$p.value, digits=3))) +
  expand_limits(y = 900) +
  labs(x="", y="", title="M3")



p7_4 <- ggplot(meta_extended[meta_extended$Timepoint == "M1", ], aes(x=Type, y=generalist_perc, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none") +
  annotate("text", x=2.2, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_perc[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M1"], 
                                                                    meta_extended$generalist_perc[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M1"], alternative = "two.sided")$p.value, digits=3))) +
  expand_limits(y = 100) +
  labs(x="", y="Generalist percent")

p7_5 <- ggplot(meta_extended[meta_extended$Timepoint == "M2", ], aes(x=Type, y=generalist_perc, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none") +
  annotate("text", x=2.2, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_perc[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M2"], 
                                                                    meta_extended$generalist_perc[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M2"], alternative = "two.sided")$p.value, digits=3))) +
  expand_limits(y = 100) +
  labs(x="", y="")

p7_6 <- ggplot(meta_extended[meta_extended$Timepoint == "M3", ], aes(x=Type, y=generalist_perc, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none") +
  annotate("text", x=2.2, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_perc[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M3"], 
                                                                    meta_extended$generalist_perc[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M3"], alternative = "two.sided")$p.value, digits=3))) +
  expand_limits(y = 100) +
  labs(x="", y="")



p7_7 <- ggplot(meta_extended[meta_extended$Timepoint == "M1", ], aes(x=Type, y=generalist_RA, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none") +
  annotate("text", x=2.2, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_RA[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M1"], 
                                                                    meta_extended$generalist_RA[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M1"], alternative = "two.sided")$p.value, digits=3))) +
  expand_limits(y = 100) +
  labs(x="", y="Generalists relative abundance")

p7_8 <- ggplot(meta_extended[meta_extended$Timepoint == "M2", ], aes(x=Type, y=generalist_RA, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none") +
  annotate("text", x=2.2, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_RA[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M2"], 
                                                                    meta_extended$generalist_RA[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M2"], alternative = "two.sided")$p.value, digits=3))) +
  expand_limits(y = 100) +
  labs(x="", y="")

p7_9 <- ggplot(meta_extended[meta_extended$Timepoint == "M3", ], aes(x=Type, y=generalist_RA, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#9966CC', '#DCD0FF')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#5E4153', '#B57EDC')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none") +
  annotate("text", x=2.2, y=100, label= paste("p-value=", signif(wilcox.test(meta_extended$generalist_RA[meta_extended$Type == "Infant" & meta_extended$Timepoint == "M3"], 
                                                                    meta_extended$generalist_RA[meta_extended$Type == "Mother" & meta_extended$Timepoint == "M3"], alternative = "two.sided")$p.value, digits=3))) +
  expand_limits(y = 100) +
  labs(x="", y="")


plot7 <- p7_1 + p7_2 + p7_3 + p7_4 + p7_5 + p7_6 + p7_7 + p7_8 + p7_9 + plot_layout(ncol = 3, guides = "collect")
plot7

##################################### EIGHTS (COMPARISON OF GENERALIST RICHNESS/ABUNDANCE/PERCENT BABIES EARLY VS LATE TIMEPOINS) GRAPH ####################################

meta_extended <- meta_extended %>% mutate(Timepoint_type = ifelse(Timepoint %in% c("M6", "M12"), "late", "early"))
meta_extended$Timepoint_type[meta_extended$Timepoint %in% c("P7", "B")] <- "pregnancy"
meta_extended$Timepoint_type <- as.factor(meta_extended$Timepoint_type)

p8_1_lm <- lmer(generalist_richness ~ Timepoint_type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Infant',])
p8_2_lm <- lmer(generalist_perc ~ Timepoint_type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Infant',])
p8_3_lm <- lmer(generalist_RA ~ Timepoint_type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Infant',])

p8_1 <- ggplot(meta_extended[meta_extended$Type == "Infant", ], aes(x=Timepoint_type, y=generalist_richness, fill=Timepoint_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d65f5f', '#dc7ec0')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#8c0800', '#a23582')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  annotate("text", x=1.5, y=900, label= paste("p-value=", "0.042 (wilcox);", "0.000197 (lm)")) +
  expand_limits(y = 900) +
  labs(x="", y="Generalist richness", title="Generalist richness")

p8_2 <- ggplot(meta_extended[meta_extended$Type == "Infant", ], aes(x=Timepoint_type, y=generalist_perc, fill=Timepoint_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d65f5f', '#dc7ec0')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#8c0800', '#a23582')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  expand_limits(y = 100) +
  annotate("text", x=1.5, y=100, label= paste("p-value=", "0.173 (wilcox);", "0.5731 (lm)")) +
  labs(x="", y="", title="Generalist percent")

p8_3 <- ggplot(meta_extended[meta_extended$Type == "Infant", ], aes(x=Timepoint_type, y=generalist_RA, fill=Timepoint_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d65f5f', '#dc7ec0')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#8c0800', '#a23582')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  annotate("text", x=1.5, y=100, label= paste("p-value=", "0.101 (wilcox);", "0.0173 (lm)")) +
  expand_limits(y = 100) +
  labs(x="", y="", title="Generalists relative abundance")


plot8 <- p8_1 + p8_2 + p8_3 + plot_layout(ncol = 3, guides = "collect")
plot8

##################################### NINGTH (COMPARISON OF GENERALIST RICHNESS/ABUNDANCE/PERCENT MOTHERS PREGNANCY VS EARLY TIMEPOINS) GRAPH ####################################

p9_1_lm <- lmer(generalist_richness ~ Timepoint_type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Mother',])
p9_2_lm <- lmer(generalist_perc ~ Timepoint_type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Mother',])
p9_3_lm <- lmer(generalist_RA ~ Timepoint_type + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = meta_extended[meta_extended$Type=='Mother',])

p9_1 <- ggplot(meta_extended[meta_extended$Type == "Mother", ], aes(x=Timepoint_type, y=generalist_richness, fill=Timepoint_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d65f5f', '#dc7ec0')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#8c0800', '#a23582')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  annotate("text", x=1.5, y=900, label= paste("p-value=0.1818 (lm)")) +
  expand_limits(y = 900) +
  labs(x="", y="Generalist richness", title="Generalist richness")

p9_2 <- ggplot(meta_extended[meta_extended$Type == "Mother", ], aes(x=Timepoint_type, y=generalist_perc, fill=Timepoint_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d65f5f', '#dc7ec0')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#8c0800', '#a23582')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  expand_limits(y = 100) +
  annotate("text", x=1.5, y=100, label= paste("p-value=0.0438 (lm)")) +
  labs(x="", y="", title="Generalist percent")

p9_3 <- ggplot(meta_extended[meta_extended$Type == "Mother", ], aes(x=Timepoint_type, y=generalist_RA, fill=Timepoint_type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c('#d65f5f', '#dc7ec0')) +
  geom_jitter(aes(color=Type), position=position_jitter(0.2), size=2, alpha= 1) +
  scale_color_manual(values=c('#8c0800', '#a23582')) +
  theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        legend.position="none",
        plot.title = element_text(size=14, hjust = 0.5)) +
  annotate("text", x=1.5, y=100, label= paste("p-value=0.0481 (lm)")) +
  expand_limits(y = 100) +
  labs(x="", y="", title="Generalists relative abundance")


plot9 <- p9_1 + p9_2 + p9_3 + plot_layout(ncol = 3, guides = "collect")
plot9


##################################### COMPARISON OF GENERALIST RICHNESS PERCENT IN MATERNAL VS INFANT SPECIFIC VIRUSES ####################################

perc_generalist_mom_specific <- sum(generalist_mom_specific$generalist) / length(generalist_mom_specific$generalist) * 100
perc_generalist_inf_specific <- sum(generalist_inf_specific$generalist) / length(generalist_inf_specific$generalist) * 100
perc_generalist_unspecific <- sum(generalist_unspecific$generalist) / length(generalist_unspecific$generalist) * 100

##################################### COMPARISON OF GENERALIST RICHNESS PERCENT IN MATERNAL VS INFANT SPECIFIC VIRUSES ####################################

assignment_classification_table <- assignments_3000_fullclassif[, c("Virus", "Host.genus")]
assignment_classification_table$phylum <- gsub('.*p__|;c__.*', "", assignment_classification_table$Host.genus)
assignment_classification_table$class <- gsub('.*c__|;o__.*', "", assignment_classification_table$Host.genus)
assignment_classification_table$order <- gsub('.*o__|;f__.*', "", assignment_classification_table$Host.genus)
assignment_classification_table$family <- gsub('.*f__|;g__.*', "", assignment_classification_table$Host.genus)
assignment_classification_table$genus <- gsub('.*g__', "", assignment_classification_table$Host.genus)

phylum_generalists <- unique(assignment_classification_table[, c("Virus", "phylum")])
class_generalists <- unique(assignment_classification_table[, c("Virus", "class")])
order_generalists <- unique(assignment_classification_table[, c("Virus", "order")])
family_generalists <- unique(assignment_classification_table[, c("Virus", "family")])

phylum_generalists <- as.data.frame(table(phylum_generalists$Virus))
class_generalists <- as.data.frame(table(class_generalists$Virus))
order_generalists <- as.data.frame(table(order_generalists$Virus))
family_generalists <- as.data.frame(table(family_generalists$Virus))

phylum_generalists <- as.data.frame(table(phylum_generalists$Freq))
class_generalists <- as.data.frame(table(class_generalists$Freq))
order_generalists <- as.data.frame(table(order_generalists$Freq))
family_generalists <- as.data.frame(table(family_generalists$Freq))

colnames(phylum_generalists) <- c("phylum_infected", "virus_number")
colnames(class_generalists) <- c("classes_infected", "virus_number")
colnames(order_generalists) <- c("orderds_infected", "virus_number")
colnames(family_generalists) <- c("familes_infected", "virus_number")

assignment_classification_table <-  merge(x=assignment_classification_table, y=generalist_df[, c("Virus", "generalist")], by="Virus", all.x = TRUE)
assignment_classification_table_generalists <- assignment_classification_table[assignment_classification_table$generalist == 1, ]
full_assingnment <- as.data.frame(table(assignment_classification_table_generalists$Host.genus))

##################################### LINKS EXPLORATION ####################################

test_df <- unique(test[, c("Virus", "Host.genus")])
df_new <- test_df %>%
  group_by(Virus) %>%
  summarise(Host.genus = list(sort(Host.genus))) 

datalist = list()

for (i in 1:nrow(df_new)){
  val = as.data.frame(t(combn(df_new$Host.genus[[i]], 2)))
  datalist[[i]] <- val
}

big_data = do.call(rbind, datalist)

big_data <- big_data %>%
  rowwise() %>%      # for each row
  mutate(links = paste(sort(c(V1, V2)), collapse = " - ")) %>%  # sort the teams alphabetically and then combine them separating with -
  ungroup()  

count_genera_pairs <- as.data.frame(table(big_data$links))
