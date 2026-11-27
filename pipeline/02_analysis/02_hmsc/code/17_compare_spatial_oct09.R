# =============================================================================
# 17_compare_spatial_oct09.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine
#           after 03_check_convergence.R for run_oct09)
# Author: Agostina Torres
#
# Understory: models with a spatially explicit plot-level effect (latent
# factors following a Gaussian process on plot coordinates; medium runs)
# compared with the main structure (non-spatial plot level; long runs), all on
# the observed data:
#   m1 spatial  (~ time x elevation) vs m1 non-spatial: whether spatial
#               structure in the residual plot effect changes the time response
#   m0 spatial  (~ time)             vs m1 non-spatial: whether a spatial plot
#               effect without elevation gives the same time response
# for occurrence and paired relative cover. Compared through Gamma for time
# and Beta for time per species; the spatial scale of each latent factor
# (Alpha) is summarised from the convergence output of 03.
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/run_oct09/understory_strict_<run_mode>_<model_id>.rds
#          pipeline/02_analysis/02_hmsc/outputs/understory_strict_medium_run_oct09_conv_alpha_overview.csv (03)
#          pipeline/02_analysis/02_hmsc/outputs/understory_strict_<run_mode>_run_oct09_conv_summary.csv (03)
# Outputs: pipeline/02_analysis/02_hmsc/outputs/run_oct09_compare/spatial_*.csv
#          pipeline/02_analysis/02_hmsc/figures/run_oct09_compare/spatial_structures.png
#          pipeline/04_robustness/figures/figS5_5_hmsc_spatial_structures.png
# =============================================================================

library(here)
library(tidyverse)
library(Hmsc)
library(hrbrthemes)
library(ggpubr)

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
res_dir    <- file.path(module_dir, "outputs")
fit_dir    <- file.path(res_dir, "run_oct09")
out_dir    <- file.path(res_dir, "run_oct09_compare")
fig_dir    <- file.path(module_dir, "figures", "run_oct09_compare")
walk(c(out_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)
sig_level  <- 0.90
sup_dir <- here("pipeline", "04_robustness", "figures")   # supplementary and diagnostic figures
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line  = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 11),
          axis.title.y = element_text(colour = "black", size = 11),
          legend.text  = element_text(colour = "black"),
          strip.text   = element_text(colour = "black", face = "bold", size = 10),
          plot.title   = element_text(colour = "black", size = 11, face = "bold"),
          plot.margin  = margin(12, 10, 10, 10))
}
trait_labels <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity",
                  SLA = "SLA", seed_mass = "Seed mass", woodinesswoody = "Woody (vs non-woody)")
trait_lab <- function(x) factor(coalesce(trait_labels[x], x), levels = rev(trait_labels))
col_struct <- c(m1_nonspatial = "black", m1_spatial = "#1A759F", m0_spatial = "#E07B39")
lab_struct <- c(m1_nonspatial = "Time x elevation, non-spatial (main, long)",
                m1_spatial    = "Time x elevation, spatial (medium)",
                m0_spatial    = "Time only, spatial (medium)")

structures <- tribble(
  ~response,    ~structure,      ~run_mode, ~model_id,
  "occurrence", "m1_nonspatial", "long",    "m1_time_elev_nsp_nos_pa",
  "occurrence", "m1_spatial",    "medium",  "m1_time_elev_nos_pa",
  "occurrence", "m0_spatial",    "medium",  "m0_time_nos_pa",
  "rel_cover",  "m1_nonspatial", "long",    "m1_time_elev_nsp_nos_rel",
  "rel_cover",  "m1_spatial",    "medium",  "m1_time_elev_nos_rel",
  "rel_cover",  "m0_spatial",    "medium",  "m0_time_nos_rel"
)

summ <- function(s) {
  tibble(mean = colMeans(s), q05 = apply(s, 2, quantile, 0.05), q95 = apply(s, 2, quantile, 0.95),
         support = colMeans(s > 0))
}

res <- pmap(structures, function(response, structure, run_mode, model_id) {
  m <- readRDS(file.path(fit_dir, paste0("understory_strict_", run_mode, "_", model_id, ".rds")))
  post <- poolMcmcChains(m$postList)
  i_t  <- which(m$covNames == "timepresent")
  B <- t(sapply(post, function(p) p$Beta[i_t, ]));  colnames(B) <- m$spNames
  G <- t(sapply(post, function(p) p$Gamma[i_t, ])); colnames(G) <- m$trNames
  G <- G[, colnames(G) != "(Intercept)", drop = FALSE]
  rm(m, post); gc()
  list(gamma = summ(G) %>% mutate(trait = colnames(G), response = response, structure = structure),
       beta  = tibble(species = colnames(B), beta_time = colMeans(B), response = response, structure = structure))
})
gamma <- bind_rows(map(res, "gamma"))
beta  <- bind_rows(map(res, "beta"))

beta_cor <- beta %>%
  pivot_wider(names_from = structure, values_from = beta_time) %>%
  group_by(response) %>%
  summarise(r_m1spatial_vs_main = round(cor(m1_spatial, m1_nonspatial), 3),
            r_m0spatial_vs_main = round(cor(m0_spatial, m1_nonspatial), 3), .groups = "drop")
cat("Beta(time) per species, correlation with the main structure:\n"); print(beta_cor)
cat("\nGamma for time:\n")
print(gamma %>% select(response, structure, trait, mean, q05, q95, support) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf)

# Spatial scale (Alpha) of the latent factors, and convergence of each model
alpha_file <- file.path(res_dir, "understory_strict_medium_run_oct09_conv_alpha_overview.csv")
if (file.exists(alpha_file)) {
  alpha <- read_csv(alpha_file, show_col_types = FALSE) %>%
    filter(model_id %in% structures$model_id)
  cat("\nSpatial scale per latent factor (proportion of samples at 0; range of chain medians, km):\n")
  print(alpha, n = Inf)
  write_csv(alpha, file.path(out_dir, "spatial_alpha.csv"))
}
conv <- map_dfr(c("long", "medium"), function(rm_) {
  f <- file.path(res_dir, paste0("understory_strict_", rm_, "_run_oct09_conv_summary.csv"))
  if (file.exists(f)) read_csv(f, show_col_types = FALSE) else NULL
}) %>%
  filter(model_id %in% structures$model_id, group %in% c("Beta", "Gamma")) %>%
  select(model_id, group, psrf_max, prop_psrf_high)
cat("\nConvergence (Beta, Gamma):\n"); print(conv, n = Inf)

write_csv(gamma, file.path(out_dir, "spatial_gamma.csv"))
write_csv(beta,  file.path(out_dir, "spatial_beta_time.csv"))
write_csv(beta_cor, file.path(out_dir, "spatial_beta_correlation.csv"))

panel <- function(resp, title) {
  gamma %>% filter(response == resp) %>%
    mutate(structure = factor(structure, levels = names(lab_struct))) %>%
    ggplot(aes(mean, trait_lab(trait), colour = structure)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = q05, xmax = q95), width = 0, linewidth = 0.8,
                  position = position_dodge(width = 0.6), orientation = "y") +
    geom_point(size = 2.2, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = col_struct, labels = lab_struct, name = NULL) +
    labs(x = "Trait effect on the response to time (Γ)", y = NULL, title = title) +
    theme_mora() + theme(legend.position = "bottom", legend.direction = "vertical")
}
fig <- ggarrange(panel("occurrence", "Occurrence (observed data)"),
                 panel("rel_cover", "Relative cover (paired)"),
                 ncol = 2, labels = c("A", "B"), common.legend = TRUE, legend = "bottom")
ggsave(file.path(fig_dir, "spatial_structures.png"), fig, width = 12, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(sup_dir, "figS5_5_hmsc_spatial_structures.png"), fig, width = 12, height = 5.5, dpi = 300, bg = "white")
cat("\nFigure written to", fig_dir, "\n")
