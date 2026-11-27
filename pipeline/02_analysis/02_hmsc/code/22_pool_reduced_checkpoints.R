# =============================================================================
# 22_pool_reduced_checkpoints.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine
#           at any time while the runs of 21_launch_reduced100.R progress)
# Author: Agostina Torres
#
# Pools the understory occurrence fits to reduced past data at checkpoints
# (number of reductions: 10, 20, 40, 60, 80, 100) and follows the stability of the
# pooled estimates as reductions accumulate.
#
# For each checkpoint reached by the finished fits, the first N fits in launch
# order (understory_reduced_datasets_100.csv: the 10 fits of run_oct09, then
# the new fits in random order) are pooled as in 14_pool_reduced_oct09.R: each
# chain is thinned so that the pooled posterior keeps about `target_samples`
# samples, every fit contributing equally. The pooled object can be passed to
# 04_hmsc_results.R for the full set of results and figures.
#
# Stability summaries per checkpoint (fixed effects, at mean elevation):
#   - species change in occurrence probability, present − past
#     (pnorm(b0 + b_time) − pnorm(b0) per posterior sample), with posterior
#     support for an increase or a decrease;
#   - trait effects on the time response (Gamma for time) with support.
# Convergence of each fit (PSRF of Beta and Gamma) is computed once and cached.
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/understory_reduced_datasets_100.csv (21)
#          fits in outputs/run_oct09/ and outputs/run_red100/
#          pipeline/02_analysis/02_hmsc/data_clean/understory_strict_hmsc_models_unfitted.rds
# Outputs: pipeline/02_analysis/02_hmsc/outputs/run_red100/
#            understory_strict_long_m1_time_elev_nsp_nos_redpoolNNN_pa.rds  (pooled fit per checkpoint)
#            red100_fit_convergence.csv          (PSRF per fit, cached)
#            red100_species_change_by_checkpoint.csv
#            red100_gamma_by_checkpoint.csv
#            red100_checkpoint_stability.csv     (agreement with the previous checkpoint)
#            red100_species_change_by_checkpoint.png
#            red100_gamma_by_checkpoint.png
# =============================================================================

library(here)
library(Hmsc)
library(coda)
library(tidyverse)

stratum     <- "understory"
pool        <- "strict"
run_mode    <- "long"
out_subdir  <- "run_red100"
base_pa     <- "m1_time_elev_nsp_nos_pa"
checkpoints <- c(10, 20, 40, 60, 80, 100)
target_samples <- 1000     # approximate size of each pooled posterior
support_thr    <- 0.90
psrf_target    <- 1.1

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
out_dir    <- file.path(module_dir, "outputs", out_subdir)
fit_path   <- function(id, subdir) file.path(module_dir, "outputs", subdir,
                                             paste0(stratum, "_", pool, "_", run_mode, "_", id, ".rds"))
out_file   <- function(name) file.path(out_dir, name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

map <- read.csv(file.path(module_dir, "data_clean", "understory_reduced_datasets_100.csv")) %>%
  arrange(launch_order) %>%
  mutate(path = fit_path(model_id, fit_subdir), done = file.exists(path))

# Fits are used in launch order; a gap (a fit still running) stops the sequence
n_ready <- if (all(map$done)) nrow(map) else which(!map$done)[1] - 1
reached <- checkpoints[checkpoints <= n_ready]
cat("Finished fits:", sum(map$done), "| consecutive in launch order:", n_ready,
    "| checkpoints reached:", paste(reached, collapse = ", "), "\n")
if (length(reached) == 0) stop("No checkpoint reached yet")

# -----------------------------------------------------------------------------
# 1. Convergence per fit (cached)
# -----------------------------------------------------------------------------
conv_file <- out_file("red100_fit_convergence.csv")
conv <- if (file.exists(conv_file)) read.csv(conv_file) else
  data.frame(model_id = character(), psrf_beta_max = numeric(), psrf_gamma_max = numeric())
todo <- setdiff(map$model_id[seq_len(n_ready)], conv$model_id)
for (id in todo) {
  m  <- readRDS(map$path[map$model_id == id])
  cc <- convertToCodaObject(m)
  conv <- rbind(conv, data.frame(
    model_id       = id,
    psrf_beta_max  = max(gelman.diag(cc$Beta,  multivariate = FALSE)$psrf[, 1]),
    psrf_gamma_max = max(gelman.diag(cc$Gamma, multivariate = FALSE)$psrf[, 1])))
  rm(m, cc); gc()
}
write.csv(conv, conv_file, row.names = FALSE)
bad <- conv$model_id[conv$psrf_beta_max > psrf_target | conv$psrf_gamma_max > psrf_target]
cat("Fits with Beta or Gamma PSRF >", psrf_target, ":", if (length(bad)) paste(bad, collapse = ", ") else "none", "\n")

# -----------------------------------------------------------------------------
# 2. Pooling and summaries per checkpoint
# -----------------------------------------------------------------------------
Y_obs <- local({
  u <- readRDS(file.path(module_dir, "data_clean", paste0(stratum, "_", pool, "_hmsc_models_unfitted.rds")))
  u$hmsc[[which(u$model_id == base_pa)]]$Y
})

pool_fits <- function(ids, paths) {
  keep_per_chain <- NULL
  post_all <- list(); Y_list <- list(); ref <- NULL
  for (k in seq_along(ids)) {
    m <- readRDS(paths[k])
    if (is.null(ref)) {
      ref <- m
      keep_per_chain <- max(2, floor(target_samples / (length(ids) * length(m$postList))))
    } else {
      stopifnot(identical(m$spNames, ref$spNames), identical(m$covNames, ref$covNames),
                identical(m$trNames, ref$trNames), identical(m$samples, ref$samples))
    }
    idx <- round(seq(m$samples / keep_per_chain, m$samples, length.out = keep_per_chain))
    post_all <- c(post_all, lapply(m$postList, function(ch) ch[idx]))
    Y_list[[ids[k]]] <- m$Y
    rm(m); gc()
  }
  pooled <- ref
  pooled$postList <- post_all
  pooled$samples  <- keep_per_chain
  pooled$thin     <- ref$thin * (ref$samples / keep_per_chain)
  pooled$Y_reduced_mean <- Reduce(`+`, Y_list) / length(Y_list)
  pooled$Y_observed     <- Y_obs
  pooled$pooled_from    <- ids
  pooled
}

summarise_draws_cols <- function(s) {
  tibble(mean = colMeans(s), lower_90 = apply(s, 2, quantile, 0.05),
         upper_90 = apply(s, 2, quantile, 0.95), p_positive = colMeans(s > 0))
}

checkpoint_summary <- function(m, n) {
  post <- poolMcmcChains(m$postList)
  B <- simplify2array(lapply(post, `[[`, "Beta"))      # covariates x species x samples
  G <- simplify2array(lapply(post, `[[`, "Gamma"))     # covariates x traits x samples
  i0 <- which(m$covNames == "(Intercept)"); it <- which(m$covNames == "timepresent")
  # Change in occurrence probability at mean elevation (elev_sc = 0)
  d <- t(pnorm(B[i0, , ] + B[it, , ]) - pnorm(B[i0, , ]))
  sp <- summarise_draws_cols(d) %>% mutate(species = m$spNames, .before = 1)
  g <- t(G[it, , ]); colnames(g) <- m$trNames
  g <- g[, colnames(g) != "(Intercept)", drop = FALSE]
  gm <- summarise_draws_cols(g) %>% mutate(trait = colnames(g), .before = 1)
  list(species = mutate(sp, n_reductions = n, n_samples = dim(B)[3]),
       gamma   = mutate(gm, n_reductions = n, n_samples = dim(B)[3]))
}

res_sp <- list(); res_g <- list()
for (n in reached) {
  pooled_id <- sprintf("m1_time_elev_nsp_nos_redpool%03d_pa", n)
  pooled_path <- fit_path(pooled_id, out_subdir)
  if (file.exists(pooled_path)) {
    pooled <- readRDS(pooled_path)
  } else {
    pooled <- pool_fits(map$model_id[seq_len(n)], map$path[seq_len(n)])
    saveRDS(pooled, pooled_path)
  }
  s <- checkpoint_summary(pooled, n)
  res_sp[[as.character(n)]] <- s$species
  res_g[[as.character(n)]]  <- s$gamma
  cat("Checkpoint", n, "reductions:", s$species$n_samples[1], "pooled samples ->", basename(pooled_path), "\n")
  rm(pooled); gc()
}

support_dir <- function(p) case_when(p >= support_thr ~ "increase", 1 - p >= support_thr ~ "decrease",
                                     TRUE ~ "not supported")
sp_all <- bind_rows(res_sp) %>% mutate(direction = support_dir(p_positive))
g_all  <- bind_rows(res_g)  %>% mutate(support = pmax(p_positive, 1 - p_positive),
                                       direction = support_dir(p_positive))
write_csv(sp_all, out_file("red100_species_change_by_checkpoint.csv"))
write_csv(g_all,  out_file("red100_gamma_by_checkpoint.csv"))

# Agreement of each checkpoint with the previous one: correlation of species
# change relative to the mean across species, and number of species whose
# supported direction changed
stab <- map_dfr(seq_along(reached)[-1], function(k) {
  a <- filter(sp_all, n_reductions == reached[k - 1]); b <- filter(sp_all, n_reductions == reached[k])
  ab <- inner_join(a, b, by = "species", suffix = c("_prev", "_now"))
  tibble(n_previous = reached[k - 1], n_reductions = reached[k],
         r_species_change = cor(ab$mean_prev - mean(ab$mean_prev), ab$mean_now - mean(ab$mean_now)),
         n_direction_changed = sum(ab$direction_prev != ab$direction_now))
})
counts <- sp_all %>% count(n_reductions, direction) %>% pivot_wider(names_from = direction, values_from = n, values_fill = 0)
cat("\nSupported species changes per checkpoint:\n"); print(counts)
if (nrow(stab) > 0) { cat("\nAgreement with the previous checkpoint:\n"); print(stab) }
write_csv(left_join(counts, stab, by = "n_reductions"), out_file("red100_checkpoint_stability.csv"))

# -----------------------------------------------------------------------------
# 3. Figures
# -----------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()))
col_n <- setNames(colorRampPalette(c("#F2A65A", "#9A3F12"))(length(reached)), reached)

sp_order <- sp_all %>% filter(n_reductions == max(reached)) %>% arrange(mean) %>% pull(species)
p_sp <- sp_all %>%
  mutate(species = factor(str_replace_all(species, "_", " "), levels = str_replace_all(sp_order, "_", " ")),
         n = factor(n_reductions)) %>%
  ggplot(aes(mean, species, colour = n)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lower_90, xmax = upper_90), position = position_dodge(width = 0.7), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.7), size = 1.2) +
  scale_colour_manual(values = col_n, name = "Reductions") +
  labs(x = "Change in occurrence probability (present − past), mean elevation", y = NULL) +
  theme(axis.text.y = element_text(face = "italic", size = 7))
ggsave(out_file("red100_species_change_by_checkpoint.png"), p_sp,
       width = 8, height = 0.16 * length(sp_order) + 2, dpi = 300, bg = "white")

p_g <- g_all %>%
  mutate(n = factor(n_reductions)) %>%
  ggplot(aes(mean, trait, colour = n)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lower_90, xmax = upper_90), position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  scale_colour_manual(values = col_n, name = "Reductions") +
  labs(x = "Trait effect on the occurrence response to time (Γ)", y = NULL)
ggsave(out_file("red100_gamma_by_checkpoint.png"), p_g, width = 7, height = 4, dpi = 300, bg = "white")

cat("\nOutputs written to", out_dir, "\n")
