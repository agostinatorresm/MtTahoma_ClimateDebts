# =============================================================================
# 24_trait_effects_elevation.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run where the final fits are
#           stored; no refitting)
# Author: Agostina Torres
#
# Trait effects on the species' change between past and present along the
# elevational gradient, from the final HMSC fits (~ time * elevation).
#
# Species responses to time have two parts, the time effect at mean elevation
# (covariate timepresent) and its change with elevation (timepresent:elev_sc),
# and traits explain both through Gamma. The effect of a trait on the change
# between past and present at a given elevation is therefore
#     Gamma[timepresent, trait] + Gamma[timepresent:elev_sc, trait] * elev_sc,
# computed for each posterior sample along the sampled elevation range. Because
# the model is linear in elevation, each trait effect changes linearly along the
# gradient; the posterior support of Gamma[timepresent:elev_sc, trait] states
# whether that change is clear. Effects are on the link scale (probit for
# occurrence, log10 relative abundance for abundance) per SD of the trait.
#
# Inputs : final fits (paths in section 0, as in 23_species_change_centered.R)
#          pipeline/02_analysis/02_hmsc/data_clean/<stratum>_strict_hmsc_inputs.rds
#            (plot elevation, to express elevation in metres)
# Outputs: pipeline/02_analysis/02_hmsc/outputs/centered/trait_effects_elevation.csv
#          pipeline/02_analysis/02_hmsc/outputs/centered/trait_elevation_interaction.csv
# Figure : 05c_trait_effects_elevation_figure.R
# =============================================================================

library(here)
library(Hmsc)
library(tidyverse)

# -----------------------------------------------------------------------------
# 0. Settings: final fits
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
fits <- tribble(
  ~stratum,     ~response,    ~path,
  "overstory",  "occurrence", file.path(module_dir, "outputs", "overstory_strict_medium_m1_time_elev_nf2_nsp_nos_pa.rds"),
  "overstory",  "abundance",  file.path(module_dir, "outputs", "overstory_strict_medium_m1_time_elev_nf2_nsp_nos_rel.rds"),
  "understory", "occurrence", file.path(module_dir, "outputs", "run_oct09", "understory_strict_long_m1_time_elev_nsp_nos_redpool_pa.rds"),
  "understory", "abundance",  file.path(module_dir, "outputs", "run_oct09", "understory_strict_long_m1_time_elev_nsp_nos_rel.rds")
)
n_grid  <- 50
out_dir <- file.path(module_dir, "outputs", "centered")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
missing <- fits$path[!file.exists(fits$path)]
if (length(missing) > 0) stop("Fits not found (edit section 0):\n", paste(missing, collapse = "\n"))

# -----------------------------------------------------------------------------
# 1. Trait effects along elevation
# -----------------------------------------------------------------------------
summ <- function(x) c(mean = mean(x), q05 = unname(quantile(x, 0.05)), q25 = unname(quantile(x, 0.25)),
                      q75 = unname(quantile(x, 0.75)), q95 = unname(quantile(x, 0.95)), p_positive = mean(x > 0))

along <- list(); inter <- list()
for (i in seq_len(nrow(fits))) {
  st <- fits$stratum[i]; resp <- fits$response[i]
  m <- readRDS(fits$path[i])
  inp <- readRDS(file.path(module_dir, "data_clean", paste0(st, "_strict_hmsc_inputs.rds")))
  elev_mu <- mean(inp$plot_covariates$elevation); elev_sd <- sd(inp$plot_covariates$elevation)
  
  it  <- which(m$covNames == "timepresent")
  itx <- which(m$covNames == "timepresent:elev_sc")
  stopifnot("Model has no time x elevation term" = length(it) == 1, length(itx) == 1)
  
  G <- simplify2array(lapply(poolMcmcChains(m$postList), `[[`, "Gamma"))   # covariates x traits x samples
  traits <- setdiff(m$trNames, "(Intercept)")
  grid <- seq(min(m$XData$elev_sc), max(m$XData$elev_sc), length.out = n_grid)
  
  for (tr in traits) {
    j <- which(m$trNames == tr)
    g0 <- G[it, j, ]; g1 <- G[itx, j, ]
    eff <- sapply(grid, function(e) g0 + g1 * e)        # samples x grid
    along[[paste(st, resp, tr)]] <- as_tibble(t(apply(eff, 2, summ))) %>%
      mutate(stratum = st, response = resp, trait = tr, elev_sc = grid,
             elevation = elev_mu + grid * elev_sd, .before = 1)
    inter[[paste(st, resp, tr)]] <- as_tibble(t(summ(g1))) %>%
      mutate(stratum = st, response = resp, trait = tr, .before = 1)
  }
  cat(st, resp, ": traits", paste(traits, collapse = ", "), "\n")
  rm(m, G); gc()
}

along <- bind_rows(along)
inter <- bind_rows(inter) %>% mutate(support = pmax(p_positive, 1 - p_positive),
                                     direction = case_when(p_positive >= 0.90 ~ "increases with elevation",
                                                           p_positive <= 0.10 ~ "decreases with elevation",
                                                           TRUE ~ "not supported"))
write_csv(along, file.path(out_dir, "trait_effects_elevation.csv"))
write_csv(inter, file.path(out_dir, "trait_elevation_interaction.csv"))
cat("\nChange of trait effects with elevation (Gamma for time x elevation):\n")
print(inter %>% select(stratum, response, trait, mean, q05, q95, support, direction), n = Inf)
