# =============================================================================
# 07_simple_glm_comparison.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine,
#           where the fitted models (.rds) are stored, after 04_hmsc_results.R
#           for the final models fitted to the observed data)
# Author: Agostina Torres
#
# Species-level change between past and present from single-species models,
# for comparison with the HMSC estimates of the final models (time x
# elevation, non-spatial plot-level effect, no sample-level effect), fitted to
# the observed data. For each species and layer:
#   occurrence : probit GLM, presence ~ time * elev_sc
#   abundance  : linear model, log10 abundance ~ time * elev_sc, fitted on the
#                samples where the species is present
# The fixed effects are those of the HMSC models; the models have no plot
# random effect and none of the parts that HMSC shares among species (traits,
# phylogeny, residual species associations). They are the same single-species
# models used for the reduced past datasets (04_robustness/06). Species are
# modelled one at a time. Data, species, samples, elevation scaling and
# abundance response are those of the HMSC models (taken from the fitted
# model objects). Because the past and present samples of a plot are treated
# as independent, the intervals are not directly comparable with those of
# HMSC; the comparison concerns the mean change and its direction.
#
# Change between past and present as in 04_hmsc_results.R:
#   occurrence : per plot, P(present) - P(past) at the plot's elevation,
#                averaged over all plots (HMSC: with the plot-level random
#                effect, delta_species_plotRE.csv)
#   abundance  : change in log10 abundance at the plot's elevation, averaged
#                over the plots where the species occurs (HMSC:
#                delta_species.csv; the plot effect cancels in the difference)
# Uncertainty: 1000 draws of the coefficients from their estimated sampling
# distribution (multivariate normal). A change is "supported" when at least
# sig_level of the draws are on one side of zero, as for HMSC.
#
# Occurrence models with very large coefficients (quasi-complete separation,
# e.g. no occurrences at one time) or that do not converge give unstable
# estimates; they are flagged.
#
# Inputs : 02_hmsc/outputs/**/<stratum>_strict_<run_mode>_<model>_{pa,<response>}.rds
#            (searched in outputs/ and its subfolders)
#          02_hmsc/outputs/<stratum>_strict_<run_mode>_<model>_<response>_res_delta_species{,_plotRE}.csv (04)
# Outputs: 02_hmsc/outputs/<stratum>_strict_glm_time_elev_<response>_delta_species.csv
#          02_hmsc/outputs/glm_vs_hmsc_<tag>_<response>_agreement.csv
#          02_hmsc/figures/<fig_subdir>/ms1__glm_species_change_<response>_<tag>.png
#          02_hmsc/figures/<fig_subdir>/ms1__hmsc_vs_glm_species_change_<response>_<tag>.png
# =============================================================================

library(here)
library(tidyverse)
library(MASS, exclude = "select")
library(hrbrthemes)
library(ggpubr)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
pool       <- "strict"
# Final models (observed data) and their run lengths
run_mode   <- c(overstory = "medium", understory = "long")
hmsc_model <- c(overstory = "m1_time_elev_nf2_nsp_nos", understory = "m1_time_elev_nsp_nos")
tag        <- "m1_final"            # added to output names
fig_subdir <- "m1_final"            # subfolder of figures/
response   <- "rel"                 # abundance response (as in the HMSC models)
strata     <- c("overstory", "understory")
n_draws    <- 1000
sig_level  <- 0.90
seed       <- 1

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
out_dir    <- file.path(module_dir, "outputs")
fig_dir    <- file.path(module_dir, "figures", fig_subdir)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Plot style (as in 05_hmsc_figures_combined.R)
if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid       = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x      = element_line(colour = "black", linewidth = 0.6),
          axis.line.y      = element_line(colour = "black", linewidth = 0.6),
          axis.ticks       = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          axis.text        = element_text(colour = "black"),
          axis.title.x     = element_text(colour = "black", size = 11),
          axis.title.y     = element_text(colour = "black", size = 11),
          legend.text      = element_text(colour = "black"),
          legend.title     = element_text(colour = "black"),
          strip.text       = element_text(colour = "black", face = "bold", size = 10),
          plot.title       = element_text(colour = "black", size = 11, face = "bold"),
          plot.subtitle    = element_text(colour = "black"),
          plot.caption     = element_text(colour = "black"),
          plot.margin      = margin(12, 10, 10, 10))
}
col_by_response <- c(occurrence = "#E07B39", log10_abundance = "#9A3F12")
col_model       <- c(HMSC = "#9A3F12", GLM = "grey45")
fold_labels <- function(x) paste0("×", signif(10^x, 2))
sp_label    <- function(x) str_replace_all(x, "_", " ")
strata_labels <- c(overstory = "Overstory", understory = "Understory")
abundance_labels <- c(abu = "abundance", rel = "relative abundance",
                      abucnt = "stem density", relcnt = "relative density")

summarise_draws <- function(d) {
  tibble(mean = colMeans(d),
         q05 = apply(d, 2, quantile, 0.05), q25 = apply(d, 2, quantile, 0.25),
         q75 = apply(d, 2, quantile, 0.75), q95 = apply(d, 2, quantile, 0.95),
         support = colMeans(d > 0)) %>%
    mutate(direction = case_when(support >= sig_level     ~ "increase",
                                 support <= 1 - sig_level ~ "decrease",
                                 TRUE                     ~ "not supported"))
}

# -----------------------------------------------------------------------------
# 1. Single-species models per layer
# -----------------------------------------------------------------------------
# Design rows for past and present at given elevations, in the column order of
# the coefficients
design_rows <- function(fit, elev) {
  nd <- function(tm) data.frame(time = factor(tm, levels = c("past", "present")), elev_sc = elev)
  f  <- formula(fit)[-2]
  list(past    = model.matrix(f, nd(rep("past", length(elev))))[, names(coef(fit)), drop = FALSE],
       present = model.matrix(f, nd(rep("present", length(elev))))[, names(coef(fit)), drop = FALSE])
}

# Draws of the change in occurrence probability from a probit GLM, averaged
# over all plots at their elevation
occ_draws <- function(fit, elev) {
  b <- mvrnorm(n_draws, coef(fit), vcov(fit))
  X <- design_rows(fit, elev)
  rowMeans(pnorm(b %*% t(X$present)) - pnorm(b %*% t(X$past)))
}

# Locate a fitted model in outputs/ or its subfolders
find_fit <- function(file) {
  hits <- list.files(out_dir, pattern = paste0("^", file, "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("Fitted model not found under outputs/: ", file)
  if (length(hits) > 1) message("Several copies of ", file, "; using ", hits[1])
  hits[1]
}

set.seed(seed)
species_elev_list <- list()
glm_results <- map_dfr(strata, function(st) {
  read_fit <- function(id) readRDS(find_fit(paste0(st, "_", pool, "_", run_mode[[st]], "_", id, ".rds")))
  m_pa  <- read_fit(paste0(hmsc_model[[st]], "_pa"))
  m_abu <- read_fit(paste0(hmsc_model[[st]], "_", response))
  
  # Samples (plot x survey)
  samples <- tibble(plot = factor(as.character(m_pa$studyDesign$plot)),
                    time = factor(as.character(m_pa$XData$time), levels = c("past", "present")),
                    elev_sc = m_pa$XData$elev_sc)
  plot_elev <- samples %>% distinct(plot, elev_sc) %>% arrange(plot)
  plots <- as.character(plot_elev$plot)
  # Plots where each species occurs in either period (for averaging abundance change)
  occ_plot <- rowsum(m_pa$Y, group = as.character(samples$plot))[plots, , drop = FALSE] > 0
  # Mean (scaled) elevation of the plots where each species occurs, for
  # ordering species in the figures
  species_elev_list[[st]] <<- tibble(stratum = st, species = colnames(occ_plot),
                                     mean_elevation = colSums(occ_plot * plot_elev$elev_sc) / colSums(occ_plot))
  
  # Occurrence: probit GLM per species
  occ <- map_dfr(m_pa$spNames, function(sp) {
    d <- mutate(samples, y = m_pa$Y[, sp])
    warn <- NULL
    fit <- withCallingHandlers(
      tryCatch(glm(y ~ time * elev_sc, family = binomial(link = "probit"), data = d),
               error = function(e) NULL),
      warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") })
    if (is.null(fit)) return(NULL)
    draws <- occ_draws(fit, plot_elev$elev_sc)
    summarise_draws(matrix(draws, ncol = 1)) %>%
      mutate(species = sp, response = "occurrence",
             flagged = !is.null(warn) || !fit$converged || any(abs(coef(fit)) > 8), .before = 1)
  })
  
  # Log10 abundance when present: linear model per species
  abu <- map_dfr(m_abu$spNames, function(sp) {
    d <- mutate(samples, y = m_abu$Y[, sp]) %>% filter(!is.na(y)) %>% droplevels()
    if (n_distinct(d$time) < 2 || nrow(d) < 6) return(NULL)
    fit <- tryCatch(lm(y ~ time * elev_sc, data = d), error = function(e) NULL)
    if (is.null(fit) || anyNA(coef(fit))) return(NULL)
    b <- mvrnorm(n_draws, coef(fit), vcov(fit))
    # Change in log10 abundance at the elevation of the occupied plots, averaged
    e <- plot_elev$elev_sc[occ_plot[, sp]]
    X <- design_rows(fit, e)
    dX <- colMeans(X$present - X$past)
    summarise_draws(matrix(drop(b %*% dX), ncol = 1)) %>%
      mutate(species = sp, response = "log10_abundance", flagged = FALSE, .before = 1)
  })
  
  res <- bind_rows(occ, abu) %>% mutate(stratum = st, .before = 1)
  write_csv(res, file.path(out_dir, paste0(st, "_", pool, "_glm_time_elev_", response, "_delta_species.csv")))
  cat(st, ": occurrence models fitted:", nrow(occ), "of", length(m_pa$spNames),
      "(flagged:", sum(occ$flagged), ") | abundance models fitted:", nrow(abu), "of",
      length(m_abu$spNames), "\n")
  res
})

# Species mean elevation (scaled; plots where the species occurs in either
# period), for ordering as in the HMSC species-change figure
species_elev <- bind_rows(species_elev_list)

# HMSC estimates (04) for the comparison: occurrence change with the plot
# random effect; log10 abundance change (plot effect cancels)
hmsc_results <- map_dfr(strata, function(st) {
  f <- function(x) file.path(out_dir, paste0(st, "_", pool, "_", run_mode[[st]], "_", hmsc_model[[st]], "_",
                                             response, "_res_", x))
  f_occ <- f("delta_species_plotRE.csv"); f_abu <- f("delta_species.csv")
  if (!file.exists(f_occ) || !file.exists(f_abu)) {
    message("HMSC results not found (run 04): ", basename(f_occ)); return(NULL)
  }
  bind_rows(read_csv(f_occ, show_col_types = FALSE),
            read_csv(f_abu, show_col_types = FALSE) %>% filter(response == "log10_abundance")) %>%
    mutate(stratum = st, .before = 1)
})

all_res <- bind_rows(mutate(glm_results, model = "GLM"),
                     mutate(hmsc_results, model = "HMSC")) %>%
  filter(response %in% c("occurrence", "log10_abundance")) %>%
  left_join(species_elev, by = c("stratum", "species")) %>%
  mutate(supported = factor(direction != "not supported", levels = c(TRUE, FALSE),
                            labels = c(paste0("Supported (≥ ", sig_level, ")"), "Not supported")))

# Agreement between GLM and HMSC
agreement <- all_res %>%
  select(stratum, species, response, model, mean, direction) %>%
  pivot_wider(names_from = model, values_from = c(mean, direction)) %>%
  filter(!is.na(mean_GLM), !is.na(mean_HMSC))
cat("\nAgreement GLM vs HMSC (correlation of species changes; same supported direction):\n")
print(agreement %>% group_by(stratum, response) %>%
        summarise(n = n(), r = round(cor(mean_GLM, mean_HMSC), 3),
                  same_direction = sum(direction_GLM == direction_HMSC), .groups = "drop"))
write_csv(agreement, file.path(out_dir, paste0("glm_vs_hmsc_", tag, "_", response, "_agreement.csv")))

# -----------------------------------------------------------------------------
# 2. Figures
# -----------------------------------------------------------------------------
abundance_label <- abundance_labels[[response]]
row_heights <- c(overstory = 0.75, understory = 1.6)

change_panel <- function(dat, s, resp, xlab, x_labels = waiver(), show_species = TRUE, compare = FALSE) {
  d   <- dat %>% filter(stratum == s, response == resp)
  ord <- dat %>% filter(stratum == s, response == "occurrence") %>%
    distinct(species, mean_elevation) %>% arrange(mean_elevation) %>% pull(species)
  d <- d %>% mutate(species = factor(sp_label(species), levels = sp_label(ord)))
  pos <- if (compare) position_dodge(width = 0.6) else "identity"
  p <- ggplot(d, aes(y = species, x = mean, alpha = supported, group = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50")
  if (compare) {
    p <- p +
      geom_linerange(aes(xmin = q05, xmax = q95, colour = model), linewidth = 0.5, position = pos) +
      geom_linerange(aes(xmin = q25, xmax = q75, colour = model), linewidth = 1.2, position = pos) +
      geom_point(aes(colour = model), size = 1.6, position = pos) +
      scale_colour_manual(values = col_model, name = NULL)
  } else {
    p <- p +
      geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.6, colour = col_by_response[[resp]]) +
      geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 1.4, colour = col_by_response[[resp]]) +
      geom_point(size = 2, colour = col_by_response[[resp]])
  }
  # Abundance: x axis limited to the range of the HMSC intervals (with a
  # margin), so that the very wide intervals of single-species models fitted
  # to few cover values do not compress the panel; intervals beyond the
  # limits are cut at the panel edge (no data removed)
  if (resp == "log10_abundance" && any(d$model == "HMSC")) {
    r <- range(c(d$q05[d$model == "HMSC"], d$q95[d$model == "HMSC"]), na.rm = TRUE)
    p <- p + coord_cartesian(xlim = r + c(-0.3, 0.3) * diff(r))
  }
  p + scale_alpha_manual(values = c(1, 0.3), name = NULL, drop = FALSE) +
    scale_x_continuous(labels = x_labels) +
    scale_y_discrete(drop = FALSE) +
    labs(x = xlab, y = if (show_species) "Species (low to high mean elevation)" else NULL,
         title = if (show_species) strata_labels[[s]] else NULL) +
    theme_mora() +
    theme(axis.text.y = if (show_species) element_text(face = "italic", colour = "black", size = if (s == "overstory") 10 else 8)
          else element_blank(),
          legend.position = "top")
}

make_figure <- function(dat, compare, name) {
  rows <- map(strata, function(s) {
    ggarrange(change_panel(dat, s, "occurrence", "Change in occurrence probability", compare = compare),
              change_panel(dat, s, "log10_abundance", paste("Change in", abundance_label, "when present"),
                           x_labels = fold_labels, show_species = FALSE, compare = compare),
              ncol = 2, widths = c(1.5, 1), legend = "none", align = "h")
  })
  fig <- ggarrange(plotlist = rows, ncol = 1, heights = row_heights[strata], labels = c("A", "B"),
                   common.legend = TRUE, legend = "top",
                   legend.grob = get_legend(change_panel(dat, strata[1], "occurrence", "", compare = compare)))
  ggsave(file.path(fig_dir, paste0("ms1__", name, "_", response, "_", tag, ".png")), fig, device = "png", dpi = 300,
         width = 6.5, height = sum(row_heights[strata]) * 4.2 + 0.6, scale = 1.8, bg = "white", limitsize = FALSE)
}

make_figure(filter(all_res, model == "GLM"), compare = FALSE, name = "glm_species_change")
if (nrow(hmsc_results) > 0) make_figure(all_res, compare = TRUE, name = "hmsc_vs_glm_species_change")
cat("\nFigures written to", fig_dir, "\n")
