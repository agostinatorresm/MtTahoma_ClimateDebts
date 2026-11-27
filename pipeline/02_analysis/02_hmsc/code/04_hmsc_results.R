# =============================================================================
# 04_hmsc_results.R
# Author: Agostina Torres
#
# Results of the HMSC occurrence and abundance models of one layer, analysed
# together as the two parts of the hurdle model. Intended for the final,
# converged run; on shorter runs the estimates are preliminary.
#
# Species-level change between surveys is expressed on the response scale as
# the predicted difference present - past:
#   - delta occurrence: change in occurrence probability, averaged over all
#     plots (occurrence model)
#   - delta abundance: change in abundance conditional on presence, for the
#     abundance response chosen in `abundance_response` (relative basal area
#     or relative cover by default; see 01_prepare_hmsc_inputs.R), averaged
#     over the plots where the species occurs in either survey
# Predictions use the fixed effects for each posterior sample, with random
# effects at zero (a typical plot), at each plot's own elevation.
#
# A posterior quantity is considered supported when at least `sig_level` of its
# posterior distribution lies on one side of zero (posterior support >= 0.90 or
# <= 0.10).
#
# Figures:
#   1. Posterior of delta occurrence and delta abundance per species
#   2. Delta occurrence and delta abundance against each trait, with the trait
#      effects on the time response (Gamma) alongside
#   3. Delta occurrence and delta abundance along the elevational gradient, per
#      species, and the predicted rate of change of the community-weighted
#      thermal and moisture affinity (CTI, CMI); CTI and CMI in the past and
#      present are shown in a separate figure
#   Abundance changes are displayed on the log10 scale (fold change); changes in
#   original units are kept in the output tables.
#   4. Residual species associations (Omega) per random level
#   5. Phylogenetic signal in species responses (rho)
#   6. Variance partitioning and explanatory power
#
# Usage  : Rscript 04_hmsc_results.R <stratum> <abundance_response> <run_mode>
#          [fit_subdir] [model_base] [occurrence model id]
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_<model_id>.rds
#          pipeline/02_analysis/02_hmsc/data_clean/<stratum>_<pool>_hmsc_inputs.rds
# Outputs: pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_<model>_res_*.csv
#          pipeline/02_analysis/02_hmsc/figures/<stratum>_<pool>_<run_mode>_<model>_res_*.png
# =============================================================================

library(here)
library(tidyverse)
library(Hmsc)
library(hrbrthemes)
library(viridis)
library(ggpubr)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
stratum  <- "overstory"   # "overstory" | "understory"
pool     <- "strict"
run_mode <- "medium"

# Abundance response analysed with the occurrence model:
#   overstory : "rel" (relative basal area), "abu" (basal area % of plot area),
#               "relcnt" (relative density), "abucnt" (stem density)
#   understory: "rel" (relative cover), "abu" (cover)
abundance_response <- "rel"

# Folder of the fitted models inside 02_hmsc/outputs ("" = outputs itself).
# Results and figures are always written to outputs/ and figures/.
fit_subdir <- ""

# Layer and abundance response from the command line, if given, e.g.
# Rscript 04_hmsc_results.R overstory relcnt
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) >= 1) stratum            <- cmd_args[1]
if (length(cmd_args) >= 2) abundance_response <- cmd_args[2]
if (length(cmd_args) >= 3) run_mode           <- cmd_args[3]
if (length(cmd_args) >= 4) fit_subdir         <- cmd_args[4]

# Occurrence and abundance models analysed together (same fixed effects), per
# layer: main structure (time x elevation, non-spatial plot level)
model_base <- c(overstory = "m1_time_elev_nf2_nsp", understory = "m1_time_elev_nsp")
# Another fixed-effect structure (e.g. "m0_time_nsp", time only) from the command line
if (length(cmd_args) >= 5) model_base[[stratum]] <- cmd_args[5]
model_pa   <- paste0(model_base[[stratum]], "_pa")
model_abu  <- paste0(model_base[[stratum]], "_", abundance_response)
# Another occurrence model analysed with the same abundance model, from the
# command line (e.g. the pooled posterior of the occurrence models fitted to
# reduced past datasets, "m1_time_elev_nsp_nos_redpool_pa"; see
# 14_pool_reduced_oct09.R)
if (length(cmd_args) >= 6) model_pa <- cmd_args[6]

sig_level <- 0.90   # posterior support threshold (one side of zero)

# Flagged species (from 01c_prelaunch_checks.R, <pool>_species_flags.csv).
# All species are in the models; flags only change how species are shown in
# the species-level figures:
#   "highlight" : flagged species marked with a dagger and open symbols
#   "omit"      : flagged species left out of the species-level figures
#   "none"      : flags ignored
# Community-level quantities (Gamma, CTI/CMI, variance partitioning) always
# use all species.
# Flags marked in the figures. Data-limitation flags that are common and mild
# (quasi-separation: coefficient partly informed by the prior; patchy: covered
# by the sensitivity analysis in Supplementary Methods S3) are listed in the
# species table only; add them here to mark them in the figures as well.
flag_display <- "highlight"
flags_used   <- c("identification", "non_climatic", "stand_dynamics",
                  "few_abundance_records")
n_grid    <- 30     # points along the elevational gradient
seed      <- 1

# Mean interval between surveys (decades), to express predicted change in
# CTI and CMI as rates comparable with the thermophilization analyses
survey_interval_decades <- 4.05

# Colours shared with the thermophilization figures (temperature / precipitation)
col_temp <- "#F26B5B"
col_prec <- "#1A759F"

# Colours for supported increases / decreases and unsupported changes
# Supported changes / effects (posterior support >= sig_level, either sign) in
# black, unsupported in light grey; the sign is read from the x-axis
col_sig <- c(increase = "black", decrease = "black", "not supported" = "grey75")

# -----------------------------------------------------------------------------
# Plot style (shared across project figures)
# -----------------------------------------------------------------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
}

# In non-interactive runs (Rscript) printed plots go to a temporary png device,
# which uses the Windows fonts registered above
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

# Panel labels placed in the top-left corner, clear of axis titles
arrange_panels <- function(..., labels = NULL) {
  ggarrange(..., labels = labels, font.label = list(size = 12, face = "bold"),
            hjust = 0, vjust = 1, label.x = 0.005, label.y = 0.995)
}

# Axis labels for changes on the log10 scale, shown as fold change
fold_labels <- function(x) paste0("×", signif(10^x, 2))

save_fig <- function(plot, file, width, height) {
  ggsave(file, plot, device = "png", dpi = 300, width = width, height = height,
         scale = 1.8, bg = "transparent", limitsize = FALSE)
}

# -----------------------------------------------------------------------------
# Paths and inputs
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
run_tag    <- paste(stratum, pool, sep = "_")
model_dir  <- file.path(module_dir, "outputs")
# Figures go to figures/ or to a subfolder set with options(mora.fig_subdir = "...")
fig_dir    <- file.path(module_dir, "figures", getOption("mora.fig_subdir", ""))
walk(c(model_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

# Output names: <stratum>_<pool>_<run_mode>_<model>_res_<name>, where <model>
# is the shared model name plus the abundance response (e.g.
# m1_time_elev_nsp_rel)
model_name <- paste0(str_remove(model_pa, "_pa$"), "_", abundance_response)
res_file <- function(name) file.path(model_dir, paste0(run_tag, "_", run_mode, "_", model_name, "_res_", name))
res_fig  <- function(name) file.path(fig_dir,   paste0(run_tag, "_", run_mode, "_", model_name, "_res_", name))

fit_dir <- if (nzchar(fit_subdir)) file.path(model_dir, fit_subdir) else model_dir
read_model <- function(id) {
  f <- file.path(fit_dir, paste0(run_tag, "_", run_mode, "_", id, ".rds"))
  if (!file.exists(f)) stop("Fitted model not found: ", f)
  readRDS(f)
}
m_pa  <- read_model(model_pa)
m_abu <- read_model(model_abu)

inp <- readRDS(file.path(module_dir, "data_clean", paste0(run_tag, "_hmsc_inputs.rds")))

abundance_labels <- list(
  overstory  = c(abu = "basal area (% of plot area)", rel = "relative basal area",
                 abucnt = "stem density (stems per m2)", relcnt = "relative density"),
  understory = c(abu = "cover (%)", rel = "relative cover",
                 abuall = "cover (%), all presences", relall = "relative cover, all presences")
)
abundance_label <- abundance_labels[[stratum]][[abundance_response]]
# Data are taken from the fitted model objects (species, samples, traits,
# scaled elevation); the inputs file supplies raw trait values and raw
# elevation only.
traits_raw <- inp$TrData_raw
missing_tr <- setdiff(m_pa$spNames, traits_raw$species)
if (length(missing_tr) > 0) {
  warning("Species of the fitted model without raw traits in the current inputs: ",
          paste(missing_tr, collapse = ", "))
}
sample_info <- tibble(plot_id = as.character(m_pa$studyDesign$plot),
                      elev_sc = if ("elev_sc" %in% names(m_pa$XData)) m_pa$XData$elev_sc else NA_real_) %>%
  left_join(inp$plot_covariates %>% select(plot_id, elevation), by = "plot_id")
has_elev <- "elev_sc" %in% m_pa$covNames
if (has_elev) {
  # Scaling of elevation used in the fitted model (elevation = mean + sd * elev_sc)
  elev_fit  <- coef(lm(elevation ~ elev_sc, data = distinct(sample_info)))
  elev_mean <- unname(elev_fit[1]); elev_sd <- unname(elev_fit[2])
} else {
  elev_mean <- mean(sample_info$elevation, na.rm = TRUE); elev_sd <- 1
}

trait_labels <- c("MAT_mean" = "Thermal affinity (°C)", "MAP_mean" = "Moisture affinity (mm)",
                  "seed_mass" = "Seed mass (mg)", "wood_density" = "Wood density (g/cm3)",
                  "SLA" = "SLA (mm2/mg)", "woodiness" = "Woodiness",
                  "woodinesswoody" = "Woody (vs non-woody)")
relabel  <- function(x, labels) ifelse(x %in% names(labels), labels[x], x)
sp_label <- function(x) str_replace_all(x, "_", " ")

# Species flags per response: a flag applies to occurrence when it affects
# "pa" or "all", and to abundance when it affects the abundance response or
# "all"
flags_file <- file.path(module_dir, "outputs", paste0(pool, "_species_flags.csv"))
species_flags <- if (flag_display != "none" && file.exists(flags_file)) {
  read_csv(flags_file, show_col_types = FALSE) %>%
    filter(stratum == !!stratum, flag %in% flags_used)
} else {
  tibble(stratum = character(), species = character(), flag = character(), affects = character())
}
if (flag_display != "none" && !file.exists(flags_file)) {
  warning("Species flags file not found (run 01c_prelaunch_checks.R): ", flags_file)
}
flagged_for <- function(resp) {
  key <- if (resp == "occurrence") "pa" else abundance_response
  species_flags %>%
    filter(str_detect(affects, paste0("\\b(", key, "|all)\\b"))) %>%
    distinct(species) %>%
    pull(species)
}
flagged_occ <- flagged_for("occurrence")
flagged_abu <- flagged_for("abundance")
flagged_any <- unique(species_flags$species)
flagged_any
sp_label_flag <- function(x) {
  lab <- sp_label(x)
  if (flag_display == "highlight") lab <- if_else(x %in% flagged_any, paste0(lab, " \u2020"), lab)
  lab
}
# Adds `flagged` (per response) and removes flagged species when omitted
apply_flags <- function(d) {
  d <- d %>%
    mutate(flagged = if_else(response == "occurrence",
                             species %in% flagged_occ, species %in% flagged_abu))
  if (flag_display == "omit") d <- filter(d, !flagged)
  d
}
shape_flag <- scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
# Posterior samples of a parameter as an array (rows x columns x samples)
post_array <- function(m, par) {
  simplify2array(lapply(poolMcmcChains(m$postList), `[[`, par))
}

# Design matrix for new covariate values, with the model's factor levels
make_X <- function(m, time, elev_sc) {
  xd <- data.frame(time = factor(time, levels = levels(m$XData$time)), elev_sc = elev_sc)
  X  <- model.matrix(m$XFormula, xd)
  stopifnot(identical(colnames(X), m$covNames))
  X
}

# Summary of posterior samples (rows = samples, columns = units)
summarise_post <- function(s) {
  tibble(mean    = colMeans(s),
         q05     = apply(s, 2, quantile, 0.05),
         q25     = apply(s, 2, quantile, 0.25),
         q75     = apply(s, 2, quantile, 0.75),
         q95     = apply(s, 2, quantile, 0.95),
         support = colMeans(s > 0)) %>%
    mutate(direction = factor(case_when(support >= sig_level     ~ "increase",
                                        support <= 1 - sig_level ~ "decrease",
                                        TRUE                     ~ "not supported"),
                              levels = c("increase", "decrease", "not supported")))
}

# Plot-level elevation and species occupancy (plots where present in either survey)
plot_elev <- sample_info %>%
  distinct(plot_id, .keep_all = TRUE) %>%
  mutate(elev_sc = coalesce(elev_sc, 0))
# Occupancy from the observed data: a pooled occurrence model fitted to reduced
# past datasets carries the observed matrix as Y_observed
Y_occ    <- if (!is.null(m_pa$Y_observed)) m_pa$Y_observed else m_pa$Y
occ_plot <- rowsum(Y_occ, group = as.character(m_pa$studyDesign$plot)) > 0
occ_plot <- occ_plot[plot_elev$plot_id, , drop = FALSE]

X_past    <- make_X(m_pa, rep("past", nrow(plot_elev)),    plot_elev$elev_sc)
X_present <- make_X(m_pa, rep("present", nrow(plot_elev)), plot_elev$elev_sc)

B_pa  <- post_array(m_pa,  "Beta")    # covariates x species x samples
B_abu <- post_array(m_abu, "Beta")
n_samples <- dim(B_pa)[3]

# -----------------------------------------------------------------------------
# 1. Species-level change: delta occurrence and delta abundance
# -----------------------------------------------------------------------------
# Delta occurrence: mean over all plots of P(present) - P(past)
delta_occ_s <- t(sapply(seq_len(n_samples), function(s) {
  colMeans(pnorm(X_present %*% B_pa[, , s]) - pnorm(X_past %*% B_pa[, , s]))
}))
colnames(delta_occ_s) <- m_pa$spNames

# Delta abundance: mean over occupied plots of 10^mu(present) - 10^mu(past),
# with mu the predicted log10 abundance conditional on presence; and the
# corresponding change on the log10 scale (multiplicative change)
occ_abu <- occ_plot[, m_abu$spNames, drop = FALSE]
delta_abu_s <- t(sapply(seq_len(n_samples), function(s) {
  d <- 10^(X_present %*% B_abu[, , s]) - 10^(X_past %*% B_abu[, , s])
  colSums(d * occ_abu) / colSums(occ_abu)
}))
delta_log_abu_s <- t(sapply(seq_len(n_samples), function(s) {
  d <- (X_present - X_past) %*% B_abu[, , s]
  colSums(d * occ_abu) / colSums(occ_abu)
}))
colnames(delta_abu_s) <- colnames(delta_log_abu_s) <- m_abu$spNames

delta_species <- bind_rows(
  summarise_post(delta_occ_s)     %>% mutate(species = m_pa$spNames,  response = "occurrence"),
  summarise_post(delta_abu_s)     %>% mutate(species = m_abu$spNames, response = "abundance"),
  summarise_post(delta_log_abu_s) %>% mutate(species = m_abu$spNames, response = "log10_abundance")
) %>%
  relocate(species, response)
delta_species %>% count(response, direction)

write_csv(delta_species, res_file("delta_species.csv"))

# Delta occurrence including the plot-level random effect: each plot keeps its
# own latent effect (Eta x Lambda, shared by the two surveys) and the change is
# averaged over plots. Without the random effect the change refers to a
# typical plot (random effects at zero); with it, to the plots actually
# sampled. The sample-level random effect is set to zero in both. Changes in
# log10 abundance are unaffected (additive effects cancel in the difference).
i_plot_rl <- which(names(m_pa$ranLevels) == "plot")
if (length(i_plot_rl) == 1) {
  post_pa     <- poolMcmcChains(m_pa$postList)
  plot_levels <- levels(m_pa$dfPi[, i_plot_rl])
  idx_plot    <- match(plot_elev$plot_id, plot_levels)
  stopifnot("plots not matched to the plot random level" = !anyNA(idx_plot))
  delta_occ_re_s <- t(sapply(seq_len(n_samples), function(s) {
    lam <- post_pa[[s]]$Lambda[[i_plot_rl]]
    if (length(dim(lam)) == 3) lam <- lam[, , 1]
    lam <- matrix(lam, ncol = length(m_pa$spNames))
    L <- post_pa[[s]]$Eta[[i_plot_rl]][idx_plot, , drop = FALSE] %*% lam
    colMeans(pnorm(X_present %*% B_pa[, , s] + L) - pnorm(X_past %*% B_pa[, , s] + L))
  }))
  colnames(delta_occ_re_s) <- m_pa$spNames
  delta_species_re <- summarise_post(delta_occ_re_s) %>%
    mutate(species = m_pa$spNames, response = "occurrence", .before = 1)
  write_csv(delta_species_re, res_file("delta_species_plotRE.csv"))
  cat("Delta occurrence with the plot random effect: mean over species",
      round(mean(delta_species_re$mean), 3), "(without:",
      round(mean(filter(delta_species, response == "occurrence")$mean), 3), ")\n")
}

# Species ordered by delta occurrence
sp_order <- delta_species %>%
  filter(response == "occurrence") %>%
  arrange(mean) %>%
  pull(species)

# Abundance change is displayed on the log10 scale (fold change), because
# changes in original units differ by orders of magnitude among species and
# are dominated by species with large, uncertain abundances
abundance_axis <- paste0("Change in ", abundance_label, " when present (fold)")

# With flag_display = "omit", a species flagged for one response only keeps
# its row (empty in the panel of the flagged response), so that the two
# panels stay aligned
delta_plot_data <- delta_species %>%
  apply_flags()
sp_order_shown <- intersect(sp_order, delta_plot_data$species)
delta_plot_data <- delta_plot_data %>%
  mutate(species = factor(sp_label_flag(species), levels = sp_label_flag(sp_order_shown)))

delta_panel <- function(resp, xlab, x_labels = waiver(), show_species = TRUE) {
  delta_plot_data %>%
    filter(response == resp) %>%
    ggplot(aes(y = species, x = mean, colour = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.4, show.legend = TRUE) +
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 1.2, show.legend = TRUE) +
    geom_point(aes(shape = flagged), size = 1.8, fill = "white", stroke = 0.8, show.legend = TRUE) +
    scale_colour_manual(values = col_sig, guide = "none") +
    shape_flag +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(labels = x_labels) +
    labs(x = xlab, y = NULL) +
    theme_mora() +
    theme(axis.text.y = if (show_species) element_text(face = "italic") else element_blank(),
          legend.position = "top")
}
p_delta_occ <- delta_panel("occurrence", "Change in occurrence probability")
p_delta_occ
p_delta_abu <- delta_panel("log10_abundance", abundance_axis, x_labels = fold_labels,
                           show_species = FALSE)
p_delta_abu

p_delta <- arrange_panels(p_delta_occ, p_delta_abu, ncol = 2, widths = c(1.45, 1),
                          labels = c("A", "B"), common.legend = TRUE, legend = "top",
                          align = "h")
p_delta

save_fig(p_delta, res_fig("delta_species.png"), width = 7, height = 0.16 * length(sp_order_shown) + 2)

# -----------------------------------------------------------------------------
# 2. Change against traits, with trait effects on the time response (Gamma)
# -----------------------------------------------------------------------------
# Gamma: effect of each trait on the species' time response (Beta for time,
# i.e. the change at mean elevation on the link scale)
gamma_time <- function(m, response) {
  G <- post_array(m, "Gamma")                 # covariates x traits x samples
  g <- t(G[which(m$covNames == "timepresent"), , ])
  colnames(g) <- m$trNames
  summarise_post(g[, colnames(g) != "(Intercept)", drop = FALSE]) %>%
    mutate(trait = setdiff(m$trNames, "(Intercept)"), response = response, .before = 1)
}
gamma_tbl <- bind_rows(gamma_time(m_pa, "occurrence"), gamma_time(m_abu, "abundance"))
gamma_tbl

write_csv(gamma_tbl, res_file("gamma_time.csv"))

model_traits <- all.vars(m_pa$TrFormula)
delta_traits <- delta_species %>%
  filter(response %in% c("occurrence", "log10_abundance")) %>%
  apply_flags() %>%
  left_join(traits_raw %>% select(species, all_of(model_traits)), by = "species")

# One panel per trait; the y-axis title is shown on the first panel of a row
trait_panel <- function(resp, trait, show_ylab = TRUE) {
  # Species without a value for this trait in the current inputs are not shown
  d <- delta_traits %>% filter(response == resp, !is.na(.data[[trait]]))
  ylab <- if (!show_ylab) NULL else if (resp == "occurrence") "Change in occurrence probability" else abundance_axis
  p <- ggplot(d, aes(x = .data[[trait]], y = mean, colour = direction)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = col_sig, guide = "none") +
    labs(x = relabel(trait, trait_labels), y = ylab) +
    theme_mora()
  if (resp == "log10_abundance") p <- p + scale_y_continuous(labels = fold_labels)
  if (is.numeric(d[[trait]])) {
    p <- p + geom_errorbar(aes(ymin = q05, ymax = q95), width = 0, alpha = 0.5) +
      geom_point(aes(shape = flagged), size = 2, stroke = 0.8, show.legend = TRUE) + shape_flag
    if (trait %in% inp$log_traits) p <- p + scale_x_log10()
  } else {
    p <- p + geom_boxplot(aes(group = .data[[trait]]), colour = "grey40",
                          outlier.shape = NA, width = 0.5) +
      geom_jitter(aes(shape = flagged), width = 0.12, size = 2, stroke = 0.8, show.legend = TRUE) + shape_flag
  }
  p
}

gamma_panel <- function(resp) {
  gamma_tbl %>%
    filter(response == resp) %>%
    mutate(trait = factor(relabel(trait, trait_labels),
                          levels = rev(relabel(setdiff(m_pa$trNames, "(Intercept)"), trait_labels)))) %>%
    ggplot(aes(y = trait, x = mean, colour = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.5) +
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 1.3) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = col_sig, guide = "none") +
    labs(x = "Gamma (trait effect on time response)", y = NULL) +
    theme_mora()
}

trait_row <- function(resp_delta, resp_gamma) {
  panels <- c(imap(model_traits, ~ trait_panel(resp_delta, .x, show_ylab = .y == 1)),
              list(gamma_panel(resp_gamma)))
  ggarrange(plotlist = panels, nrow = 1, legend = "none",
            widths = c(rep(1, length(model_traits)), 1.3))
}
row_occ <- trait_row("occurrence", "occurrence")
row_occ
row_abu <- trait_row("log10_abundance", "abundance")
row_abu

fig_traits <- arrange_panels(row_occ, row_abu, ncol = 1, labels = c("A", "B"),
                             common.legend = TRUE, legend = "bottom",
                             legend.grob = get_legend(trait_panel("occurrence", model_traits[1]) +
                                                        theme(legend.position = "bottom")))
fig_traits

save_fig(fig_traits, res_fig("delta_traits_gamma.png"),
         width = 2.2 * (length(model_traits) + 1.3), height = 5.5)

# -----------------------------------------------------------------------------
# 3. Change along the elevational gradient
# -----------------------------------------------------------------------------
if (has_elev) {
  elev_grid <- seq(min(plot_elev$elev_sc), max(plot_elev$elev_sc), length.out = n_grid)
  Xg_past    <- make_X(m_pa, rep("past", n_grid),    elev_grid)
  Xg_present <- make_X(m_pa, rep("present", n_grid), elev_grid)

  # Species elevational range (plots where present in either survey); the
  # curves are restricted to this range
  sp_range <- map_dfr(colnames(occ_plot), function(sp) {
    e <- plot_elev$elev_sc[occ_plot[, sp]]
    tibble(species = sp, elev_min = min(e), elev_max = max(e))
  })

  delta_grad <- function(B, Xp, Xa, transform, response) {
    arr <- simplify2array(lapply(seq_len(dim(B)[3]), function(s) {
      transform(Xp %*% B[, , s]) - transform(Xa %*% B[, , s])
    }))                                        # grid x species x samples
    map_dfr(seq_len(dim(arr)[2]), function(j) {
      summarise_post(t(arr[, j, ])) %>%
        mutate(species = dimnames(B)[[2]][j], elev_sc = elev_grid, .before = 1)
    }) %>%
      mutate(response = response)
  }
  dimnames(B_pa)[[2]]  <- m_pa$spNames
  dimnames(B_abu)[[2]] <- m_abu$spNames

  delta_elev <- bind_rows(
    delta_grad(B_pa,  Xg_present, Xg_past, pnorm,    "occurrence"),
    delta_grad(B_abu, Xg_present, Xg_past, identity, "log10_abundance")
  ) %>%
    left_join(sp_range, by = "species") %>%
    filter(elev_sc >= elev_min, elev_sc <= elev_max) %>%
    mutate(elevation = elev_sc * elev_sd + elev_mean) %>%
    apply_flags() %>%
    left_join(select(traits_raw, species, MAT_mean), by = "species")

  write_csv(delta_elev, res_file("delta_elevation_species.csv"))

  elev_panel <- function(resp, ylab) {
    delta_elev %>%
      filter(response == resp) %>%
      ggplot(aes(x = elevation, y = mean, group = species, colour = MAT_mean)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_line(aes(linetype = flagged), linewidth = 0.6, alpha = 0.8) +
      scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
      scale_colour_viridis_c(option = "viridis", name = "Thermal\naffinity (°C)") +
      labs(x = "Elevation (m)", y = ylab) +
      theme_mora()
  }
  p_elev_occ <- elev_panel("occurrence", "Change in occurrence probability")
  p_elev_occ
  p_elev_abu <- elev_panel("log10_abundance", abundance_axis) +
    scale_y_continuous(labels = fold_labels)
  p_elev_abu

  # Predicted community-weighted affinities (CTI, CMI) from the occurrence
  # model: sum_j(p_j * affinity_j) / sum_j(p_j), per posterior sample
  # Species of the fitted model with a raw affinity in the inputs
  aff    <- tibble(species = m_pa$spNames) %>%
    left_join(select(traits_raw, species, MAT_mean, MAP_mean), by = "species")
  keep_j <- which(!is.na(aff$MAT_mean) & !is.na(aff$MAP_mean))
  if (length(keep_j) < nrow(aff)) {
    message("CTI/CMI from the occurrence model without (no affinity in current inputs): ",
            paste(aff$species[-keep_j], collapse = ", "))
  }
  aff <- aff[keep_j, ]
  cwm <- map_dfr(seq_len(n_samples), function(s) {
    p_a <- pnorm(Xg_past %*% B_pa[, keep_j, s])
    p_p <- pnorm(Xg_present %*% B_pa[, keep_j, s])
    tibble(sample = s, elev_sc = elev_grid,
           CTI_past = as.vector(p_a %*% aff$MAT_mean) / rowSums(p_a),
           CTI_present = as.vector(p_p %*% aff$MAT_mean) / rowSums(p_p),
           CMI_past = as.vector(p_a %*% aff$MAP_mean) / rowSums(p_a),
           CMI_present = as.vector(p_p %*% aff$MAP_mean) / rowSums(p_p))
  }) %>%
    pivot_longer(-c(sample, elev_sc), names_to = c("index", "time"), names_sep = "_")

  cwm_summary <- cwm %>%
    group_by(index, time, elev_sc) %>%
    summarise(mean = mean(value), q05 = quantile(value, 0.05),
              q95 = quantile(value, 0.95), .groups = "drop") %>%
    mutate(elevation = elev_sc * elev_sd + elev_mean)
  # Change (present - past) per posterior sample, also expressed per decade
  cwm_change <- cwm %>%
    pivot_wider(names_from = time, values_from = value) %>%
    mutate(change = present - past) %>%
    group_by(index, elev_sc) %>%
    summarise(mean = mean(change), q05 = quantile(change, 0.05),
              q95 = quantile(change, 0.95), support = mean(change > 0), .groups = "drop") %>%
    mutate(elevation = elev_sc * elev_sd + elev_mean,
           across(c(mean, q05, q95), ~ .x / survey_interval_decades, .names = "rate_{.col}"))
  cwm_change

  write_csv(cwm_summary, res_file("cwm_elevation.csv"))
  write_csv(cwm_change, res_file("cwm_change_elevation.csv"))

  cwm_panel <- function(ix, colour, ylab) {
    cwm_summary %>%
      filter(index == ix) %>%
      ggplot(aes(x = elevation, y = mean, colour = time, fill = time)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2, colour = NA) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c(past = alpha(colour, 0.45), present = colour),
                          labels = c(past = "~1975", present = "~2015"), name = NULL) +
      scale_fill_manual(values = c(past = colour, present = colour),
                        labels = c(past = "~1975", present = "~2015"), name = NULL) +
      labs(x = "Elevation (m)", y = ylab) +
      theme_mora() +
      theme(legend.position = "top")
  }
  p_cti <- cwm_panel("CTI", col_temp, "Predicted CTI (°C)")
  p_cti
  p_cmi <- cwm_panel("CMI", col_prec, "Predicted CMI (mm)")
  p_cmi

  # Rate of change of the predicted community-weighted affinities along the
  # gradient (posterior mean and 90% interval)
  rate_panel <- function(ix, colour, ylab) {
    cwm_change %>%
      filter(index == ix) %>%
      ggplot(aes(x = elevation, y = rate_mean)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_ribbon(aes(ymin = rate_q05, ymax = rate_q95), fill = colour, alpha = 0.25) +
      geom_line(colour = colour, linewidth = 1) +
      labs(x = "Elevation (m)", y = ylab) +
      theme_mora()
  }
  p_cti_rate <- rate_panel("CTI", col_temp, "Predicted CTI change (°C per decade)")
  p_cti_rate
  p_cmi_rate <- rate_panel("CMI", col_prec, "Predicted CMI change (mm per decade)")
  p_cmi_rate

  fig_elev <- arrange_panels(p_elev_occ, p_elev_abu, p_cti_rate, p_cmi_rate,
                             ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
  fig_elev

  fig_cwm <- arrange_panels(p_cti, p_cmi, ncol = 2, labels = c("A", "B"))
  fig_cwm

  save_fig(fig_elev, res_fig("delta_elevation.png"), width = 7, height = 6)
  save_fig(fig_cwm, res_fig("cwm_past_present.png"), width = 7, height = 3.2)
}

# -----------------------------------------------------------------------------
# 4. Residual species associations (Omega)
# -----------------------------------------------------------------------------
association_panel <- function(m, response) {
  assoc <- computeAssociations(m)
  imap(assoc, function(a, k) {
    level <- names(m$ranLevels)[k]
    cor_m <- a$mean
    keep  <- a$support >= sig_level | a$support <= 1 - sig_level
    ord   <- hclust(as.dist(1 - cor_m))$order
    sp    <- sp_label(rownames(cor_m))[ord]
    df <- as_tibble(cor_m, rownames = "sp1") %>%
      pivot_longer(-sp1, names_to = "sp2", values_to = "correlation") %>%
      mutate(supported = as.vector(t(keep)),
             value = if_else(supported & sp1 != sp2, correlation, NA_real_),
             sp1 = factor(sp_label(sp1), levels = sp),
             sp2 = factor(sp_label(sp2), levels = rev(sp)),
             level = level, response = response)
    p <- ggplot(df, aes(x = sp1, y = sp2, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_viridis_c(option = "viridis", limits = c(-1, 1), na.value = "grey95",
                           name = "Residual\ncorrelation") +
      labs(x = NULL, y = NULL,
           title = paste0(str_to_sentence(response), ": ", level, " level"),
           subtitle = paste0("Coloured: posterior support >= ", sig_level)) +
      theme_mora() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = 7),
            axis.text.y = element_text(face = "italic", size = 7),
            axis.line = element_blank())
    list(plot = p, data = df)
  })
}
assoc_pa  <- association_panel(m_pa,  "occurrence")
assoc_abu <- association_panel(m_abu, "abundance")

assoc_data <- bind_rows(map(c(assoc_pa, assoc_abu), "data"))
write_csv(assoc_data, res_file("residual_associations.csv"))

fig_assoc <- ggarrange(plotlist = map(c(assoc_pa, assoc_abu), "plot"),
                       ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
fig_assoc

save_fig(fig_assoc, res_fig("residual_associations.png"),
         width = 0.12 * length(m_pa$spNames) * 2 + 3, height = 0.12 * length(m_pa$spNames) * 2 + 2)

# -----------------------------------------------------------------------------
# 5. Phylogenetic signal (rho)
# -----------------------------------------------------------------------------
# rho = proportion of the residual variation in species responses (beyond
# traits) that is structured by phylogeny (0 = none, 1 = Brownian motion)
rho_post <- bind_rows(
  tibble(response = "occurrence", rho = sapply(poolMcmcChains(m_pa$postList),  `[[`, "rho")),
  tibble(response = "abundance",  rho = sapply(poolMcmcChains(m_abu$postList), `[[`, "rho"))
)
rho_summary <- rho_post %>%
  group_by(response) %>%
  summarise(mean = round(mean(rho), 3), median = round(median(rho), 3),
            q05 = round(quantile(rho, 0.05), 3), q95 = round(quantile(rho, 0.95), 3),
            prob_rho_above_0 = round(mean(rho > 0), 3), .groups = "drop")
rho_summary

p_rho <- ggplot(rho_post, aes(x = rho, fill = response)) +
  geom_histogram(binwidth = 0.02, boundary = 0, colour = "white", alpha = 0.8,
                 position = "identity") +
  scale_fill_viridis_d(option = "viridis", begin = 0.2, end = 0.8, name = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Phylogenetic signal (rho)", y = "Posterior samples") +
  theme_mora() +
  theme(legend.position = "top")
p_rho

write_csv(rho_summary, res_file("phylogenetic_signal.csv"))
save_fig(p_rho, res_fig("phylogenetic_signal.png"), width = 4, height = 3)

# -----------------------------------------------------------------------------
# 6. Variance partitioning and explanatory power
# -----------------------------------------------------------------------------
vp_tables <- imap(list(occurrence = m_pa, abundance = m_abu), function(m, response) {
  group <- ifelse(m$covNames %in% c("timepresent", "timepresent:elev_sc"), 2, 1)
  names <- c("Elevation", "Time")
  if (!has_elev) { group <- rep(1, length(m$covNames)); names <- "Time" }
  vp <- computeVariancePartitioning(m, group = group, groupnames = names)
  set.seed(seed)
  pred <- computePredictedValues(m, expected = TRUE)
  fit  <- evaluateModelFit(hM = m, predY = pred)
  list(
    pred = pred,
    vp = as_tibble(vp$vals, rownames = "component") %>%
      pivot_longer(-component, names_to = "species", values_to = "proportion") %>%
      mutate(response = response),
    r2_traits = tibble(response = response, covariate = names(vp$R2T$Beta),
                       r2_traits = round(vp$R2T$Beta, 3)),
    fit = as_tibble(fit) %>% mutate(species = m$spNames, response = response, .before = 1)
  )
})
vp_long   <- bind_rows(map(vp_tables, "vp"))
r2_traits <- bind_rows(map(vp_tables, "r2_traits"))
fit_tbl   <- bind_rows(map(vp_tables, "fit"))

vp_long %>% group_by(response, component) %>%
  summarise(mean_proportion = round(mean(proportion), 3), .groups = "drop")
r2_traits
fit_tbl %>% group_by(response) %>%
  summarise(across(where(is.numeric), ~ round(mean(.x, na.rm = TRUE), 3)))

p_vp <- vp_long %>%
  mutate(species = factor(sp_label(species), levels = sp_label(rev(sp_order)))) %>%
  ggplot(aes(x = species, y = proportion, fill = component)) +
  geom_col(width = 0.85) +
  facet_wrap(~ response, ncol = 1) +
  scale_fill_viridis_d(option = "viridis", name = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of explained variance") +
  theme_mora() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic", size = 7),
        legend.position = "top")
p_vp

write_csv(vp_long, res_file("variance_partitioning.csv"))
write_csv(r2_traits, res_file("variance_explained_by_traits.csv"))
write_csv(fit_tbl, res_file("model_fit.csv"))
save_fig(p_vp, res_fig("variance_partitioning.png"),
         width = 0.15 * length(m_pa$spNames) + 3, height = 5)

# -----------------------------------------------------------------------------
# 7. Posterior predictive check: observed vs predicted, past and present
# -----------------------------------------------------------------------------
# Predicted values for every sample (plot x survey), with all random effects,
# from the posterior (computePredictedValues, expected values: occurrence
# probability; log10 abundance conditional on presence). Per species and
# survey:
#   occurrence : observed prevalence vs predicted mean occurrence probability
#   abundance  : observed mean log10 abundance vs predicted mean, over the
#                samples where the species is present
# and the change between surveys (present - past), observed vs predicted.
time_s <- as.character(m_pa$XData$time)
# For a pooled occurrence model fitted to reduced past datasets, the observed
# values are the mean over the reduced datasets (Y_reduced_mean)
ppc_one <- function(m, pred, response) {
  Y <- if (!is.null(m$Y_reduced_mean)) m$Y_reduced_mean else m$Y
  map_dfr(c("past", "present"), function(tt) {
    rows <- which(time_s == tt)
    obs  <- colMeans(Y[rows, , drop = FALSE], na.rm = TRUE)
    # mean prediction over the same cells as the observed mean
    mask <- !is.na(Y[rows, , drop = FALSE])
    post <- apply(pred[rows, , , drop = FALSE], 3, function(p) colSums(p * mask) / colSums(mask))
    tibble(species = m$spNames, response = response, time = tt, observed = obs,
           pred_mean = rowMeans(post),
           pred_q05 = apply(post, 1, quantile, 0.05),
           pred_q95 = apply(post, 1, quantile, 0.95))
  })
}
ppc <- bind_rows(ppc_one(m_pa,  vp_tables$occurrence$pred, "occurrence"),
                 ppc_one(m_abu, vp_tables$abundance$pred,  "log10_abundance"))
ppc_delta <- ppc %>%
  select(species, response, time, observed, pred_mean) %>%
  pivot_wider(names_from = time, values_from = c(observed, pred_mean)) %>%
  transmute(species, response,
            observed_change  = observed_present - observed_past,
            predicted_change = pred_mean_present - pred_mean_past)
write_csv(ppc, res_file("posterior_predictive_check.csv"))
write_csv(ppc_delta, res_file("posterior_predictive_check_change.csv"))
cat("Posterior predictive check, correlation observed vs predicted:\n")
print(ppc %>% group_by(response, time) %>%
        summarise(r = round(cor(observed, pred_mean, use = "complete.obs"), 3), .groups = "drop"))
print(ppc_delta %>% group_by(response) %>%
        summarise(r_change = round(cor(observed_change, predicted_change, use = "complete.obs"), 3)))

col_time <- c(past = "#F4A259", present = "#B5501C")
resp_lab <- c(occurrence = "Occurrence (prevalence)", log10_abundance = paste0("Abundance when present (log10 ", abundance_label, ")"))
p_ppc_level <- ppc %>%
  mutate(response = factor(resp_lab[response], levels = resp_lab)) %>%
  ggplot(aes(x = observed, y = pred_mean, colour = time)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = pred_q05, ymax = pred_q95), width = 0, alpha = 0.5) +
  geom_point(size = 1.8) +
  facet_wrap(~ response, scales = "free") +
  scale_colour_manual(values = col_time, name = NULL) +
  labs(x = "Observed", y = "Predicted (posterior mean, 90% interval)") +
  theme_mora() + theme(legend.position = "top")
p_ppc_change <- ppc_delta %>%
  mutate(response = factor(resp_lab[response], levels = resp_lab)) %>%
  ggplot(aes(x = observed_change, y = predicted_change)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.8, colour = "grey20") +
  facet_wrap(~ response, scales = "free") +
  labs(x = "Observed change (present - past)", y = "Predicted change") +
  theme_mora()
fig_ppc <- arrange_panels(p_ppc_level, p_ppc_change, ncol = 1, labels = c("A", "B"))
save_fig(fig_ppc, res_fig("posterior_predictive_check.png"), width = 6.5, height = 6.5)

# -----------------------------------------------------------------------------
# 8. Species-level change from the full model (all random effects)
# -----------------------------------------------------------------------------
# Change between surveys from the predictions of the full model
# (computePredictedValues, expected values, fixed effects plus the plot- and
# sample-level random effects), for each posterior sample: per plot, predicted
# present minus predicted past, averaged over all plots (occurrence
# probability) or over the plots where the species occurs in either survey
# (log10 abundance when present). Unlike the changes in section 1, which refer
# to a typical plot (random effects at zero), these refer to the plots actually
# sampled and include the plot x survey (sample-level) variation captured by
# the latent factors.
sample_plot <- as.character(m_pa$studyDesign$plot)
plots_both  <- intersect(sample_plot[time_s == "past"], sample_plot[time_s == "present"])
row_past    <- match(plots_both, sample_plot[time_s == "past"])
row_present <- match(plots_both, sample_plot[time_s == "present"])
idx_past    <- which(time_s == "past")[row_past]
idx_present <- which(time_s == "present")[row_present]
occ_full    <- occ_plot[match(plots_both, plot_elev$plot_id), , drop = FALSE]

full_delta <- function(m, pred, response, occupied = NULL) {
  d <- pred[idx_present, , , drop = FALSE] - pred[idx_past, , , drop = FALSE]   # plots x species x samples
  s <- if (is.null(occupied)) {
    apply(d, 3, colMeans)
  } else {
    w <- occupied[, m$spNames, drop = FALSE]
    apply(d, 3, function(x) colSums(x * w) / colSums(w))
  }
  summarise_post(t(s)) %>% mutate(species = m$spNames, response = response, .before = 1)
}
delta_species_full <- bind_rows(
  full_delta(m_pa,  vp_tables$occurrence$pred, "occurrence"),
  full_delta(m_abu, vp_tables$abundance$pred,  "log10_abundance", occupied = occ_full)
)
write_csv(delta_species_full, res_file("delta_species_full.csv"))
cat("Species change from the full model (all random effects):", length(plots_both), "plots with both surveys\n")
print(delta_species_full %>% count(response, direction))
