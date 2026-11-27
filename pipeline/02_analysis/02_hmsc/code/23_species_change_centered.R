# =============================================================================
# 23_species_change_centered.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run where the final fits are
#           stored; no refitting)
# Author: Agostina Torres
#
# Species change between past and present relative to the mean change across
# species ("centered" change), from the final HMSC fits.
#
# Changes shared by all species are not interpretable in the understory: the
# reduction of the past data removes the common decline in occurrence, and
# relative cover increases for most species because fewer species have cover
# in the present. Centering removes this common component and shows which
# species increased or decreased relative to the others. It is done within each
# posterior sample (species change minus the mean change across species in that
# sample), so the credible intervals and the posterior support of the centered
# change carry the uncertainty of the common component.
#
# Species change (as in 04_hmsc_results.R):
#   occurrence: mean over all plots of P(present) − P(past), full model
#               (fixed effects and plot-level random effect);
#   abundance : mean over the plots where the species occurs of the change in
#               log10 relative abundance (the plot-level random effect is the
#               same at both times and cancels out).
#
# Inputs : final fits (paths in section 0)
# The figure is made by 05b_species_change_centered_figure.R from the table.
# Outputs: pipeline/02_analysis/02_hmsc/outputs/centered/species_change_centered.csv
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
support_thr <- 0.90
out_dir <- file.path(module_dir, "outputs", "centered")
walk(out_dir, dir.create, recursive = TRUE, showWarnings = FALSE)
missing <- fits$path[!file.exists(fits$path)]
if (length(missing) > 0) stop("Fits not found (edit section 0):\n", paste(missing, collapse = "\n"))

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------
# Plot covariates (one row per plot) and design matrices for past and present
plot_design <- function(m) {
  pl <- as.character(m$studyDesign$plot)
  xd <- m$XData %>% mutate(plot = pl) %>% distinct(plot, .keep_all = TRUE)
  mk <- function(t) {
    nd <- xd %>% mutate(time = factor(t, levels = levels(m$XData$time)))
    X <- model.matrix(m$XFormula, nd)
    stopifnot(identical(colnames(X), m$covNames))
    X
  }
  list(plot = xd$plot, X_past = mk("past"), X_present = mk("present"))
}

# Plot-level random effect (Eta x Lambda) for one posterior sample, rows in the
# order of `plots`
plot_effect <- function(m, s, plots) {
  r <- which(m$rLNames == "plot")
  lev <- levels(m$dfPi[, r])
  E <- s$Eta[[r]] %*% s$Lambda[[r]]
  E[match(plots, lev), , drop = FALSE]
}

# Species change per posterior sample (samples x species)
change_draws <- function(m, response, occupied = NULL) {
  dz <- plot_design(m)
  post <- poolMcmcChains(m$postList)
  out <- t(sapply(post, function(s) {
    if (response == "occurrence") {
      re <- plot_effect(m, s, dz$plot)
      colMeans(pnorm(dz$X_present %*% s$Beta + re) - pnorm(dz$X_past %*% s$Beta + re))
    } else {
      d <- (dz$X_present - dz$X_past) %*% s$Beta
      occ <- occupied[dz$plot, m$spNames, drop = FALSE]
      colSums(d * occ) / colSums(occ)
    }
  }))
  colnames(out) <- m$spNames
  out
}

summarise_change <- function(d) {
  tibble(species = colnames(d), mean = colMeans(d),
         q05 = apply(d, 2, quantile, 0.05), q25 = apply(d, 2, quantile, 0.25),
         q75 = apply(d, 2, quantile, 0.75), q95 = apply(d, 2, quantile, 0.95),
         support = colMeans(d > 0)) %>%
    mutate(direction = case_when(support >= support_thr ~ "increase",
                                 support <= 1 - support_thr ~ "decrease",
                                 TRUE ~ "not supported"))
}

# -----------------------------------------------------------------------------
# 2. Raw and centered species change
# -----------------------------------------------------------------------------
res <- list()
for (st in unique(fits$stratum)) {
  m_pa  <- readRDS(fits$path[fits$stratum == st & fits$response == "occurrence"])
  m_abu <- readRDS(fits$path[fits$stratum == st & fits$response == "abundance"])
  # Plots where each species occurs at either time (observed data)
  Y_occ <- if (!is.null(m_pa$Y_observed)) m_pa$Y_observed else m_pa$Y
  occupied <- rowsum(Y_occ, group = as.character(m_pa$studyDesign$plot)) > 0
  
  for (resp in c("occurrence", "abundance")) {
    m <- if (resp == "occurrence") m_pa else m_abu
    d <- change_draws(m, resp, occupied)
    d_c <- d - rowMeans(d)                      # centered within each sample
    res[[paste(st, resp)]] <- bind_rows(
      summarise_change(d)   %>% mutate(version = "raw"),
      summarise_change(d_c) %>% mutate(version = "centered")) %>%
      mutate(stratum = st, response = resp, .before = 1)
    cat(st, resp, ": mean change across species =", signif(mean(rowMeans(d)), 3), "\n")
  }
  rm(m_pa, m_abu); gc()
}
change <- bind_rows(res)
write_csv(change, file.path(out_dir, "species_change_centered.csv"))
cat("\nSupported centered changes:\n")
print(change %>% filter(version == "centered") %>% count(stratum, response, direction) %>%
        pivot_wider(names_from = direction, values_from = n, values_fill = 0))

cat("\nTable written to", out_dir, "(figure: 05b_species_change_centered_figure.R)\n")
