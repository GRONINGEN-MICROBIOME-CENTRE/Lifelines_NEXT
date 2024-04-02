library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(UpSetR)



setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data")

df_maqsood <- as.data.frame(read_tsv("maqsood_stats.tsv"))

df_maqsood <- df_maqsood %>%
  mutate(cohort = "maqsood",
         input_files = "pe",
         clean_reads_comb = out_clmpf_paired + out_clmpf_unm1 + out_clmpf_unm2,
         clean_reads_se = NA,
         in_reads_se = NA,
         dedup_efficiency = ((in_clmpf_paired + in_clmpf_unm1 + in_clmpf_unm2) - (out_clmpf_paired + out_clmpf_unm1 + out_clmpf_unm2))/(in_clmpf_paired + in_clmpf_unm1 + in_clmpf_unm2)
         )

df_liang <- as.data.frame(read_tsv("liang_stats.tsv"))

df_liang <- df_liang %>%
  mutate(cohort = "liang",
         input_files = "pe",
         clean_reads_comb = out_clmpf_paired + out_clmpf_unm1 + out_clmpf_unm2,
         clean_reads_se = NA,
         in_reads_se = NA,
         dedup_efficiency = ((in_clmpf_paired + in_clmpf_unm1 + in_clmpf_unm2) - (out_clmpf_paired + out_clmpf_unm1 + out_clmpf_unm2))/(in_clmpf_paired + in_clmpf_unm1 + in_clmpf_unm2)
  )


df_walters_se <- as.data.frame(read_tsv("walters_read_stats_se.tsv"))

df_walters_se <- df_walters_se %>%
  mutate(cohort = "walters",
         input_files = "se",
         clean_reads_comb = out_clmpf_single,
         clean_reads_se = out_clmpf_single,
         in_reads_se = in_reads_comb,
         dedup_efficiency = (in_clmpf_single - out_clmpf_single)/in_clmpf_single
  )

df_walters_pe <- as.data.frame(read_tsv("walters_read_stats_pe.tsv"))

df_walters_pe <- df_walters_pe %>%
  mutate(cohort = "walters",
         in_reads_comb = in_reads_comb_single + in_reads_comb_paired,
         input_files = "pe + se",
         clean_reads_comb = out_clmpf_single + out_clmpf_paired + out_clmpf_unm1 + out_clmpf_unm2,
         clean_reads_se = out_clmpf_single,
         in_reads_se = in_reads_comb_single,
         dedup_efficiency = ((in_clmpf_single + in_clmpf_paired + in_clmpf_unm1 + in_clmpf_unm2) - (out_clmpf_single + out_clmpf_paired + out_clmpf_unm1 + out_clmpf_unm2))/(in_clmpf_single + in_clmpf_paired + in_clmpf_unm1 + in_clmpf_unm2)
  )


df_walters_assembly <- as.data.frame(read_tsv("walters_assembly_stats.tsv"))

df_garmaeva_qc1 <- as.data.frame(read_tsv("garmaeva_read_stats.tsv"))
df_garmaeva_qc2 <- as.data.frame(read_tsv("garmaeva_unmatched_count.tsv"))

df_garmaeva_qc <- merge(x=df_garmaeva_qc1, y=df_garmaeva_qc2, by="sample_id")

df_garmaeva_qc <- df_garmaeva_qc %>%
  mutate(cohort = "garmaeva",
         input_files = "pe",
         clean_reads_comb = out_paired + out_unmatched,
         clean_reads_se = NA,
         in_reads_se = NA,
         dedup_efficiency = NA
  )

df_garmaeva_as <- as.data.frame(read_tsv("garmaeva_assembly_stats.tsv"))

df_shah_qc <- as.data.frame(read_tsv("shah_qc_stats.tsv"))

df_shah_qc <- df_shah_qc %>%
  mutate(cohort = "shah",
         input_files = "pe + um",
         clean_reads_comb = in_reads_comb,
         clean_reads_se = NA,
         in_reads_se = NA,
         dedup_efficiency = NA
  )

df_shah_as <- as.data.frame(read_tsv("shah_assembly_stats.tsv"))

## Combining all the data into one dataframe

# Preparing Walters tables

df_walters_pe_as <- merge(x=df_walters_pe, y=df_walters_assembly, by="sample_id")
df_walters_se_as <- merge(x=df_walters_se, y=df_walters_assembly, by="sample_id")

# Preparing Garmaeva tables

df_garmaeva <- merge(x=df_garmaeva_qc, y=df_garmaeva_as, by="sample_id")

# Preparing Shah tables

df_shah <- merge(x=df_shah_qc, y=df_shah_as, by="sample_id")

# Preparing Shah data

# Combining everything in one table

columns_of_interest <- c("sample_id", "cohort", "input_files", "in_reads_comb", "in_reads_se", "clean_reads_comb", "clean_reads_se", "dedup_efficiency", "contigs_total", "contigs_1000", "N50")

df_maqsood_m <- df_maqsood[columns_of_interest]
df_liang_m <- df_liang[columns_of_interest]
df_walters_se_m <- df_walters_se_as[columns_of_interest]
df_walters_pe_m <- df_walters_pe_as[columns_of_interest]
df_garmaeva_m <- df_garmaeva[columns_of_interest]
df_shah_m <- df_shah[columns_of_interest]

qc_assembly_stats <- rbind(df_maqsood_m, df_liang_m, df_walters_se_m, df_walters_pe_m, df_garmaeva_m, df_shah_m)

qc_assembly_stats_melt <- melt(qc_assembly_stats[c(1, 2, 4, 6, 9, 10)], id = c("sample_id", "cohort")) 

variable.labs <- c("Raw reads", "Clean reads", "Contigs total", "Contigs longer 1000bp")
names(variable.labs) <- c("in_reads_comb", "clean_reads_comb", "contigs_total", "contigs_1000")


ggplot(qc_assembly_stats_melt[qc_assembly_stats_melt$variable %in% c("in_reads_comb", "clean_reads_comb"), ], aes(x=value, fill=cohort)) +
  geom_histogram(alpha = 0.2, binwidth = 0.05) +
  geom_density(aes(y = ..count../15), alpha = 0.2) +
  labs(x = "Number of reads per sample combined (log10)", y = "Sample number") +
  scale_x_log10() +
  facet_grid(variable ~ ., 
             labeller = labeller(variable = variable.labs)) +
  theme_bw() +
  theme(plot.title = element_text(size=20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))

ggplot(qc_assembly_stats_melt[qc_assembly_stats_melt$variable %in% c("contigs_total", "contigs_1000"), ], aes(x=value, fill=cohort)) +
  geom_histogram(alpha = 0.2, binwidth = 0.05) +
  geom_density(aes(y = ..count../15), alpha = 0.2) + 
  labs(x = "Number of contigs per sample (log10)", y = "Sample number") +
  scale_x_log10() +
  facet_grid(variable ~ ., 
             labeller = labeller(variable = variable.labs)) +
  theme_bw() +
  theme(plot.title = element_text(size=20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))

ggplot(qc_assembly_stats, aes(x=dedup_efficiency, fill=cohort)) +
  geom_histogram(alpha = 0.2, binwidth = 0.01) + 
  geom_density(aes(y = ..count../100), alpha = 0.2) +
  labs(x = "Deduplication efficiency", y = "Sample number") +
  theme_bw() +
  theme(plot.title = element_text(size=20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))


## Adding the number of discovered viruses per sample

to_garmaeva <- as.data.frame(read_tsv("table_of_origin_garmaeva"))
to_shah <- as.data.frame(read_tsv("table_of_origin_shah"))
to_walters <- as.data.frame(read_tsv("table_of_origin_walters"))
to_liang <- as.data.frame(read_tsv("table_of_origin_liang"))
to_maqsood <- as.data.frame(read_tsv("table_of_origin_maqsood"))
to_full <- rbind(to_garmaeva, to_shah, to_walters, to_liang, to_maqsood)


to_full_per_sample <- to_full
to_full_per_sample$total_viruses_discovered <- 1
to_full_per_sample$V1 <- sub("_NODE.*", "", to_full_per_sample$V1)
to_full_per_sample <- summarise_all(group_by(to_full_per_sample, V1), sum)

colnames(to_full_per_sample)[1] <- "sample_id"

qc_assembly_vd_stats <- merge(x=qc_assembly_stats, y=to_full_per_sample, by="sample_id", all = T)
qc_assembly_vd_stats[c("DeepVirFinder", "geNomad", "VIBRANT", "VirSorter2", 
                       "total_viruses_discovered")][is.na(qc_assembly_vd_stats[c("DeepVirFinder", "geNomad", "VIBRANT", "VirSorter2", 
                                                                                 "total_viruses_discovered")])] <- 0


summary_per_tool <- as.data.frame(matrix(NA, nrow=4, ncol=2))
summary_per_tool$V1 <- colnames(to_full[,c(2:5)])
for (i in summary_per_tool$V1) {
  summary_per_tool[summary_per_tool$V1==i,]$V2 <- sum(to_full[,i])
}



listInput <- list(geNomad=to_full[to_full$geNomad==1,]$V1,
                  VIBRANT=to_full[to_full$VIBRANT==1,]$V1,
                  DeepVirFinder=to_full[to_full$DeepVirFinder==1,]$V1,
                  VirSorter2=to_full[to_full$VirSorter2==1,]$V1)

pdf('C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\plots\\VD_redundant_tools_overlap.pdf', width=18/2.54, height=10/2.54)
upset(fromList(listInput), order.by = "freq", sets.bar.color = "#C00000", 
      number.angles = 20,
      sets.x.label = "N detected virus contigs", scale.sets = "identity",
      text.scale = c(1, 1, 1, 1, 1, 0.75))
dev.off()
