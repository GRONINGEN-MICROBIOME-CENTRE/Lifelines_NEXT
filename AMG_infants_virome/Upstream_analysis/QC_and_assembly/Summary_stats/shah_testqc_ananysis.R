library(readr)
library(ggplot2)
library(reshape2) 

setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data")

df_shah <- as.data.frame(read_tsv("shahtest_read_stats.tsv"))
df_shah$KTtotal_perc <- ((df_shah$KTri + df_shah$KTrbo) / df_shah$in_reads_comb) * 100  # KT refers to the adaptor trimming step
df_shah$KTremoved_perc <- (df_shah$KTtr / df_shah$in_reads_comb) * 100

df_shah$QT_perc <- ((df_shah$in_reads_comb - df_shah$KTtr) - ((df_shah$QTrpaired1*2) + df_shah$QTrunm1 + df_shah$QTrunm2)) / (df_shah$in_reads_comb - df_shah$KTtr) * 100  # KT refers to the quality trimming step

df_shah$HRr_perc <- (((df_shah$QTrpaired1*2) + df_shah$QTrunm1 + df_shah$QTrunm2) - ((df_shah$HRrpaired1*2) + df_shah$HRrunm1 + df_shah$HRrunm2)) / ((df_shah$QTrpaired1*2) + df_shah$QTrunm1 + df_shah$QTrunm2) * 100  # KT refers to the quality trimming step

df_shah$dedup_efficiency <- ((df_shah$in_clmpf_paired + df_shah$in_clmpf_unm1 + df_shah$in_clmpf_unm2) - (df_shah$out_clmpf_paired + df_shah$out_clmpf_unm1 + df_shah$out_clmpf_unm2))/(df_shah$in_clmpf_paired + df_shah$in_clmpf_unm1 + df_shah$in_clmpf_unm2) * 100

df_shah_ed <- df_shah[c(1, 18, 19, 20, 21, 22)]
df_shah_ed_melt <- melt(df_shah_ed, id = c("sample_id")) 

# New facet label names for variables
variable.labs <- c("Adapter trimmig: reads altered", "Adapter trimmig: reads removed", "Quality trimmed reads", "Human reads removed", "Dedupication efficiency")
names(variable.labs) <- c("KTtotal_perc", "KTremoved_perc", "QT_perc", "HRr_perc", "dedup_efficiency")

ggplot(df_shah_ed_melt, aes(x=value)) +
  geom_histogram(fill="lightblue", color="black", binwidth = 0.1) + 
  labs(x = "Percent reads %", y = "Frequency (total number of tested samples: 20)") +
  facet_grid(variable ~ ., 
             labeller = labeller(variable = variable.labs)) +
  theme_bw() +
  theme(plot.title = element_text(size=20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))
