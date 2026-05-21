library(MOFA2)
library(clue)      # solve_LSAP
library(dplyr)
library(tidyr)
library(here)
library(purrr)
library(parallel)

Ks <- 2:12
bootstraps <- 1:15

model_dir <- here("80percent_bootstrap/trained_models/")

# helper function
match_factors_by_weights <- function(W_ref, W_cmp) {

  shared_features <- intersect(rownames(W_ref), rownames(W_cmp))
  if (length(shared_features) == 0) {
    stop("No shared features between models")
  }

  C <- abs(cor(
    W_ref[shared_features, ],
    W_cmp[shared_features, ],
    use = "pairwise.complete.obs" # pearson correlation is default
  ))

  assignment <- solve_LSAP(C, maximum = TRUE)

  tibble(
    factor_ref = seq_len(ncol(W_ref)),
    factor_cmp = assignment,
    correlation = C[cbind(seq_len(ncol(W_ref)), assignment)]
  )
}


#parallel version
n_cores <- 4

weight_stability_list <- mclapply(
  Ks,
  function(K) {

    message("Processing K = ", K)

    # ---- load all bootstraps for this K ----
    models <- lapply(bootstraps, function(b) {
      file <- file.path(
        model_dir,
        sprintf("MEFISTO_K%d_boot%d.hdf5", K, b)
      )
      m <- load_model(file, remove_inactive_factors = FALSE)
      W <- get_weights(m)[[1]]
      W
    })

    names(models) <- bootstraps

    # ---- all pairwise comparisons ----
    res <- list()
    idx <- 1

    for (i in 1:(length(bootstraps) - 1)) {
      for (j in (i + 1):length(bootstraps)) {

        b1 <- bootstraps[i]
        b2 <- bootstraps[j]

        matched <- match_factors_by_weights(
          models[[i]],
          models[[j]]
        )

        res[[idx]] <- matched %>%
          mutate(
            K = K,
            bootstrap1 = b1,
            bootstrap2 = b2
          )

        idx <- idx + 1
      }
    }

    bind_rows(res)
  },
  mc.cores = n_cores
)

#saveRDS(weight_stability_list, file = here("80percent_bootstrap/weight_stability_test.rds")) # as of 9 Feb 2026

weight_stability_list <- readRDS(file = here("80percent_bootstrap/weight_stability_test.rds"))


weight_stability_df <- bind_rows(weight_stability_list)



compute_Z_stability_for_K <- function(
    K,
    bootstraps,
    model_dir
) {
  
  message("Processing Z stability for K = ", K)
  
  # ---- load factor scores for all bootstraps ----
  Z_models <- lapply(bootstraps, function(b) {
    
    file <- file.path(
      model_dir,
      sprintf("MEFISTO_K%d_boot%d.hdf5", K, b)
    )
    
    if (!file.exists(file)) {
      stop("Missing file: ", file)
    }
    
    m <- load_model(file, remove_inactive_factors = FALSE)
    
    # assumes single view
    Z <- get_factors(m, factors = "all")[[1]]
    Z
  })
  
  names(Z_models) <- bootstraps
  
  # ---- helper: match factors by Z ----
  match_factors_by_scores <- function(Z1, Z2) {
    
    shared_samples <- intersect(rownames(Z1), rownames(Z2))
    
    if (length(shared_samples) < 5) {
      stop("Too few shared samples between bootstraps")
    }
    
    C <- abs(cor(
      Z1[shared_samples, , drop = FALSE],
      Z2[shared_samples, , drop = FALSE],
      use = "pairwise.complete.obs"
    ))
    
    assignment <- solve_LSAP(C, maximum = TRUE)
    
    tibble(
      factor_ref = seq_len(ncol(Z1)),
      factor_cmp = assignment,
      correlation_Z = C[cbind(seq_len(ncol(Z1)), assignment)]
    )
  }
  
  # ---- all pairwise comparisons ----
  res <- list()
  idx <- 1
  
  for (i in 1:(length(bootstraps) - 1)) {
    for (j in (i + 1):length(bootstraps)) {
      
      b1 <- bootstraps[i]
      b2 <- bootstraps[j]
      
      matched <- match_factors_by_scores(
        Z_models[[i]],
        Z_models[[j]]
      )
      
      res[[idx]] <- matched %>%
        mutate(
          K = K,
          bootstrap1 = b1,
          bootstrap2 = b2
        )
      
      idx <- idx + 1
    }
  }
  
  bind_rows(res)
}

library(parallel)

n_cores <- 5

score_stability_list <- mclapply(
  Ks,
  compute_Z_stability_for_K,
  bootstraps = bootstraps,
  model_dir = model_dir,
  mc.cores = n_cores
)

#saveRDS(score_stability_list, file = here("data_calculated/factor_value_corrs_80boot.rds")) # as of 9 Feb 2026

score_stability_list <- readRDS(file = here("data_calculated/factor_value_corrs_80boot.rds"))

score_stability_df <- dplyr::bind_rows(score_stability_list)

combined_stability <- weight_stability_df %>%
  left_join(
    score_stability_df,
    by = c(
      "K",
      "bootstrap1",
      "bootstrap2",
      "factor_ref",
      "factor_cmp"
    )
  )

library(dplyr)
library(tidyr)
library(ggplot2)

plot_df <- combined_stability %>%
  select(K, factor_ref, correlation, correlation_Z) %>%
  pivot_longer(
    cols = c(correlation, correlation_Z),
    names_to = "type",
    values_to = "corr"
  ) %>%
  mutate(
    type = recode(
      type,
      correlation   = "Factor weights (W)",
      correlation_Z = "Factor scores (Z)"
    )
  )


mean_df <- plot_df %>%
  group_by(K, factor_ref, type) %>%
  summarize(
    mean_corr = mean(corr, na.rm = TRUE),
    .groups = "drop"
  )

stability_plt <- ggplot(plot_df, aes(x = factor_ref, y = corr, color = type, shape = type)) +
  
  # individual bootstrap correlations
  geom_point(
    alpha = 0.6,
    size = 1,
    position = position_jitter(width = 0.15, height = 0)
  ) +
  
  # mean lines (one for W, one for Z)
  geom_line(
    data = mean_df,
    aes(x = factor_ref, y = mean_corr, group = type),
    linewidth = 1
  ) +
  
  # mean points with white fill
  geom_point(
    data = mean_df,
    aes(x = factor_ref, y = mean_corr, color = type),
    shape = 21,          # <-- key change
    size = 1,
    stroke = 1.2,
    fill = "white"
  ) +

  facet_wrap(~ K, scales = "free_x") +
  
  geom_hline(
    yintercept = 0.6,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  scale_x_continuous(
    breaks = 1:12,
    limits = c(1, 12)
  ) +
  
  labs(
    title = "MEFISTO factor stability across subsampling splits",
    x = "factor index",
    y = "correlation across bootstraps",
    color = "Estimated Parameter",
    shape = "Estimated Parameter"
  ) +
  
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5))

stability_plt

ggsave(
  filename = here("80percent_bootstrap/stability_results.pdf"),
  plot = stability_plt,
  width = 210,
  height = 148.5,
  units = "mm"
)

ggsave(
  filename = here("80percent_bootstrap/mefisto_factor_stability_top_half_A4.png"),
  plot = stability_plt,
  width = 210,
  height = 148.5,
  units = "mm"
)
