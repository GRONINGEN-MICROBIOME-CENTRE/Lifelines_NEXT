#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
id <- as.integer(args[1])

Ks <- c(3, 6, 9, 12)

grid <- data.frame(K = Ks)

# critical fix:
this <- grid[id, , drop = FALSE]

cat("Job:", id, " → K =", this$K, "\n")
system(sprintf("Rscript run_mefisto_full.R %d", this$K))

