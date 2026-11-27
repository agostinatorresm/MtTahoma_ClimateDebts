# =============================================================================
# 28_pool_compare_sparse12.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (after 27_launch_sparse12_reduced.R)
# Author: Agostina Torres
#
# Pools the 10 understory occurrence fits without the 12 sparsely sampled
# species (run_sparse12) and compares them with the main occurrence model fitted
# to the same 10 reduced past datasets (run_oct09). Both are pooled in the same
# way as the main model (each chain thinned so that the pooled posterior keeps
# about `target_samples` samples, every fit contributing equally).
#
# Comparison (fixed effects, at mean elevation):
#   - species change in occurrence probability, present − past, for the species
#     present in both models: correlation of posterior means (centred on the
#     mean across species, as the common change is not interpreted) and number
#     of species whose supported direction differs;
#   - trait effects on the time response (Gamma for time): posterior mean, 90%
#     interval and support in each model.
#
# Inputs : outputs/run_sparse12/ (10 medium fits; understory_sparse12_reduced_datasets.csv)
#          outputs/run_oct09/   (10 long fits; understory_reduced_datasets_oct09.csv)
# Outputs: outputs/run_sparse12/
#            understory_strict_medium_m1_time_elev_nsp_nos_sparse12_pool010_pa.rds (pooled fit)
#            sparse12_species_change_comparison.csv
#            sparse12_gamma_comparison.csv
#            sparse12_comparison_summary.csv
#            sparse12_comparison.png
# =============================================================================

library(here)
library(Hmsc)
library(tidyverse)
library(patchwork)
library(hrbrthemes)

stratum <- "understory"; pool <- "strict"
target_samples <- 1000
support_thr    <- 0.90

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
out_dir    <- file.path(module_dir, "outputs", "run_sparse12")
fit_path   <- function(id, subdir, mode) file.path(module_dir, "outputs", subdir, paste0(stratum, "_", pool, "_", mode, "_", id, ".rds"))

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.text = element_text(colour = "black", size = 11),
          axis.title.x = element_text(colour = "black", size = 12, hjust = 0.5),
          axis.title.y = element_text(colour = "black", size = 12, hjust = 0.5),
          plot.title = element_text(colour = "black", size = 13, face = "bold"))
}

# -----------------------------------------------------------------------------
# 1. Pooling (as in 22_pool_reduced_checkpoints.R)
# -----------------------------------------------------------------------------
pool_fits <- function(paths) {
  post_all <- list(); ref <- NULL; keep <- NULL
  for (k in seq_along(paths)) {
    m <- readRDS(paths[k])
    if (is.null(ref)) { ref <- m; keep <- max(2, floor(target_samples / (length(paths) * length(m$postList)))) }
    else stopifnot(identical(m$spNames, ref$spNames), identical(m$covNames, ref$covNames),
                   identical(m$trNames, ref$trNames), identical(m$samples, ref$samples))
    idx <- round(seq(m$samples / keep, m$samples, length.out = keep))
    post_all <- c(post_all, lapply(m$postList, function(ch) ch[idx]))
    rm(m); invisible(gc())
  }
  ref$postList <- post_all; ref$thin <- ref$thin * (ref$samples / keep); ref$samples <- keep
  ref
}

sparse_map <- read.csv(file.path(module_dir, "data_clean", "understory_sparse12_reduced_datasets.csv"))
main_map   <- read.csv(file.path(module_dir, "data_clean", "understory_reduced_datasets_oct09.csv"))
sparse_paths <- fit_path(sparse_map$model_id, "run_sparse12", "medium")
main_paths   <- fit_path(main_map$model_id, "run_oct09", "long")
stopifnot("Not all sparse12 fits finished" = all(file.exists(sparse_paths)),
          "Not all run_oct09 fits found"   = all(file.exists(main_paths)))

pooled_sparse_file <- file.path(out_dir, "understory_strict_medium_m1_time_elev_nsp_nos_sparse12_pool010_pa.rds")
m_sparse <- if (file.exists(pooled_sparse_file)) readRDS(pooled_sparse_file) else { p <- pool_fits(sparse_paths); saveRDS(p, pooled_sparse_file); p }
m_main   <- pool_fits(main_paths)

# -----------------------------------------------------------------------------
# 2. Species change and trait effects
# -----------------------------------------------------------------------------
summ <- function(s) tibble(mean = colMeans(s), lower_90 = apply(s, 2, quantile, 0.05),
                           upper_90 = apply(s, 2, quantile, 0.95), p_positive = colMeans(s > 0))
direction <- function(p) case_when(p >= support_thr ~ "increase", 1 - p >= support_thr ~ "decrease", TRUE ~ "not supported")
estimates <- function(m, model) {
  post <- poolMcmcChains(m$postList)
  B <- simplify2array(lapply(post, `[[`, "Beta")); G <- simplify2array(lapply(post, `[[`, "Gamma"))
  i0 <- which(m$covNames == "(Intercept)"); it <- which(m$covNames == "timepresent")
  d <- t(pnorm(B[i0, , ] + B[it, , ]) - pnorm(B[i0, , ]))
  g <- t(G[it, , ]); colnames(g) <- m$trNames; g <- g[, colnames(g) != "(Intercept)", drop = FALSE]
  list(species = summ(d) %>% mutate(species = m$spNames, model = model, .before = 1),
       gamma   = summ(g) %>% mutate(trait = colnames(g), model = model, .before = 1))
}
e_main <- estimates(m_main, "main"); e_sparse <- estimates(m_sparse, "sparse12")

sp_cmp <- inner_join(e_main$species, e_sparse$species, by = "species", suffix = c("_main", "_sparse12")) %>%
  mutate(direction_main = direction(p_positive_main), direction_sparse12 = direction(p_positive_sparse12))
g_cmp <- bind_rows(e_main$gamma, e_sparse$gamma) %>%
  mutate(support = pmax(p_positive, 1 - p_positive), direction = direction(p_positive))
write_csv(sp_cmp, file.path(out_dir, "sparse12_species_change_comparison.csv"))
write_csv(g_cmp, file.path(out_dir, "sparse12_gamma_comparison.csv"))

centre <- function(x) x - mean(x)
summary_tab <- tibble(
  n_species_compared = nrow(sp_cmp),
  r_species_change   = cor(centre(sp_cmp$mean_main), centre(sp_cmp$mean_sparse12)),
  n_direction_differs = sum(sp_cmp$direction_main != sp_cmp$direction_sparse12),
  n_opposite_supported = sum((sp_cmp$direction_main == "increase" & sp_cmp$direction_sparse12 == "decrease") |
                             (sp_cmp$direction_main == "decrease" & sp_cmp$direction_sparse12 == "increase")))
write_csv(summary_tab, file.path(out_dir, "sparse12_comparison_summary.csv"))
cat("\nSpecies change, main vs without 12 sparse species:\n"); print(summary_tab)
cat("\nSpecies whose supported direction differs:\n")
print(filter(sp_cmp, direction_main != direction_sparse12) %>%
        transmute(species, mean_main = round(mean_main, 3), direction_main,
                  mean_sparse12 = round(mean_sparse12, 3), direction_sparse12), n = Inf)
cat("\nTrait effects on the time response (Gamma):\n")
print(g_cmp %>% transmute(trait, model, mean = signif(mean, 3), lower_90 = signif(lower_90, 3),
                          upper_90 = signif(upper_90, 3), support = round(support, 2), direction) %>%
        arrange(trait, model), n = Inf)

# -----------------------------------------------------------------------------
# 3. Figure
# -----------------------------------------------------------------------------
pA <- ggplot(sp_cmp, aes(mean_main, mean_sparse12)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey50") +
  geom_errorbar(aes(ymin = lower_90_sparse12, ymax = upper_90_sparse12), colour = "grey75", width = 0) +
  geom_errorbarh(aes(xmin = lower_90_main, xmax = upper_90_main), colour = "grey75", height = 0) +
  geom_point(size = 1.8, colour = "grey20") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3.8,
           label = sprintf("r = %.2f (n = %d)", summary_tab$r_species_change, summary_tab$n_species_compared)) +
  labs(title = "A. Species change in occurrence", x = "Main model", y = "Without 12 sparse species") + theme_mora()
pB <- g_cmp %>%
  mutate(model = factor(recode(model, main = "Main model", sparse12 = "Without 12 sparse species"),
                        levels = c("Main model", "Without 12 sparse species"))) %>%
  ggplot(aes(mean, trait, colour = model)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey50") +
  geom_linerange(aes(xmin = lower_90, xmax = upper_90), position = position_dodge(width = 0.5), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.5), size = 2.2) +
  scale_colour_manual(values = c("Main model" = "grey20", "Without 12 sparse species" = "#1A759F"), name = NULL) +
  labs(title = "B. Trait effects on the time response", x = "Γ (time), 90% interval", y = NULL) +
  theme_mora() + theme(legend.position = "bottom")
fig_file <- file.path(out_dir, "sparse12_comparison.png")
ggsave(fig_file, pA + pB + plot_layout(widths = c(1, 1.1)), width = 12, height = 5.5, dpi = 300, bg = "white")
cat("\nFigure written to", fig_file, "\n")
