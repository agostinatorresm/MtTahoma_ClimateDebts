# =============================================================================
# 05_hmsc_figures_combined.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the local computer, after
#           04_hmsc_results.R; copy the <stratum>_..._res_*.csv tables from the
#           computing machine to pipeline/02_analysis/02_hmsc/outputs/ first)
# Author: Agostina Torres
#
# Usage (final models):
#   Rscript 05_hmsc_figures_combined.R medium,long m1_time_elev_nf2_nsp_nos_rel \
#           m1_time_elev_nsp_nos_redpool_rel m1_final
#
# Combined overstory + understory figures from the result tables written by
# 04_hmsc_results.R (one run of 04 per layer, same run mode and abundance
# response). Nothing is recomputed from the models.
#
# Figures:
#   1. Species change: change in occurrence probability and fold change in
#      abundance when present, per species, ordered by mean elevation;
#      random effects at zero (1) and full model (1b; 1c occurrence only).
#   2. Change along the elevational gradient: per-species change in
#      occurrence and abundance, and the predicted rate of change of the
#      community thermal index (CTI) from the occurrence model.
#   3. Trait effects on the time response (Gamma), occurrence and abundance
#      (main text; larger text, layer names written vertically with the tree
#      and herb icons). This figure is opened at the end of the script.
#   4. Residual species associations (Omega), sample and plot levels; one
#      figure per level with overstory/understory x occurrence/abundance.
#   5. Variance partitioning (mean proportion per component) and
#      phylogenetic signal (rho).
#   6. Occurrence change against abundance change per species, coloured by
#      thermal affinity.
#   7. Posterior predictive check: observed vs predicted prevalence and
#      abundance when present per species, past and present surveys.
#
# A change is "supported" when at least 90% of its posterior lies on one side
# of zero (as in 04).
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_<model>_res_*.csv
#          pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/
#            data_clean/species_climatic_affinity.csv
# Outputs: pipeline/02_analysis/02_hmsc/figures/<tag>/ms1__fig*_<pool>_<run_mode>_<tag>.png
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(patchwork)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
pool     <- "strict"
run_mode <- "medium"

# Model name used by 04 for each layer: <fixed-effect structure>_<abundance response>
models <- c(overstory  = "m0_time_nf2_nsp_rel",
            understory = "m0_time_nsp_rel")
tag <- "m0_rel"       # added to the figure names
fig_prefix <- "ms1__"  # prefix of the figure files

# Run mode, models and tag from the command line, if given:
# Rscript 05_hmsc_figures_combined.R <run_mode> <overstory model> <understory model> <tag>
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) >= 1) run_mode <- cmd_args[1]
# Run mode per layer: one value for both, or "overstory,understory" (e.g.
# "medium,long" when the final overstory and understory runs differ in length)
run_modes <- strsplit(run_mode, ",")[[1]]
run_modes <- setNames(rep_len(run_modes, 2), c("overstory", "understory"))
run_mode  <- paste(unique(run_modes), collapse = "-")   # used in the figure names
if (length(cmd_args) >= 3) models   <- c(overstory = cmd_args[2], understory = cmd_args[3])
if (length(cmd_args) >= 4) tag      <- cmd_args[4]
sig_level <- 0.90

abundance_labels <- c(abu = "abundance", rel = "relative abundance",
                      abucnt = "stem density", relcnt = "relative density")
abundance_label <- abundance_labels[[str_extract(models[["overstory"]], "[a-z]+$")]]

strata_labels <- c(overstory = "Overstory", understory = "Understory")
# Supported changes / effects (posterior support >= sig_level, either sign) in
# black, unsupported in light grey; the sign is read from the x-axis
col_sig <- c(increase = "black", decrease = "black", "not supported" = "grey75")
col_temp <- "#F26B5B"   # coral: temperature / CTI (as in the thermophilization figures)

# -----------------------------------------------------------------------------
# Style
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
fold_labels <- function(x) paste0("×", signif(10^x, 2))
sp_label    <- function(x) str_replace_all(x, "_", " ")
save_fig <- function(plot, name, width, height) {
  f <- file.path(fig_dir, paste0(fig_prefix, name, "_", pool, "_", run_mode, "_", tag, ".png"))
  ggsave(f, plot, device = "png", dpi = 300, width = width, height = height,
         scale = 1.6, bg = "white", limitsize = FALSE)
  invisible(f)
}

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
res_dir    <- file.path(module_dir, "outputs")
# Figures go to a subfolder named after the tag (e.g. figures/m1_final/), or to
# the subfolder set with options(mora.fig_subdir = "...")
fig_dir    <- file.path(module_dir, "figures", getOption("mora.fig_subdir", tag))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_res <- function(name) {
  imap_dfr(models, function(model, s) {
    f <- file.path(res_dir, paste0(s, "_", pool, "_", run_modes[[s]], "_", model, "_res_", name))
    if (!file.exists(f)) {
      message("Not found (run 04 for ", s, "): ", basename(f))
      return(NULL)
    }
    read_csv(f, show_col_types = FALSE) %>% mutate(stratum = s, .before = 1)
  })
}

affinity <- read_csv(here("pipeline", "01_data_prep", "05_species_distribution",
                          "02_integrate_gbif_chelsa", "data_clean", "species_climatic_affinity.csv"),
                     show_col_types = FALSE) %>%
  transmute(species = str_replace_all(species, " ", "_"), MAT_mean, MAP_mean)

# Mean elevation of the plots where each species occurs (either survey), used to
# order species in the species-change figure
cleaning_dir <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
plot_info  <- read_csv(cleaning_dir("outputs", "plot_info.csv"), show_col_types = FALSE)
plot_key   <- intersect(c("plot_id", "plot", "Plot", "PlotID"), names(plot_info))[1]
elev_col   <- intersect(c("elevation", "elevation_m", "elev", "elev_m"), names(plot_info))[1]
plot_elev <- plot_info %>%
  transmute(plot_id = as.character(.data[[plot_key]]), elevation = as.numeric(.data[[elev_col]])) %>%
  filter(!is.na(elevation)) %>%
  distinct(plot_id, .keep_all = TRUE)
species_elev <- map_dfr(names(models), function(s) {
  read_csv(cleaning_dir("outputs", paste0(s, "_long.csv")), show_col_types = FALSE) %>%
    filter(occurrence > 0) %>%
    mutate(plot_id = as.character(plot_id)) %>%
    left_join(plot_elev, by = "plot_id") %>%
    group_by(species) %>%
    summarise(mean_elevation = mean(elevation, na.rm = TRUE), .groups = "drop") %>%
    mutate(species = str_replace_all(species, " ", "_"), stratum = s)
})

delta     <- read_res("delta_species.csv") %>% left_join(affinity, by = "species") %>%
  left_join(species_elev, by = c("stratum", "species"))
gamma     <- read_res("gamma_time.csv")
delta_el  <- read_res("delta_elevation_species.csv")
cwm_ch    <- read_res("cwm_change_elevation.csv")
assoc     <- read_res("residual_associations.csv")
vp        <- read_res("variance_partitioning.csv")
rho       <- read_res("phylogenetic_signal.csv")

stopifnot("delta_species tables not found for any stratum" = nrow(delta) > 0)
direction_levels <- c("increase", "decrease", "not supported")
delta <- delta %>% mutate(direction = factor(direction, levels = direction_levels),
                          stratum_lab = factor(strata_labels[stratum], levels = strata_labels))

# -----------------------------------------------------------------------------
# 1. Species change
# -----------------------------------------------------------------------------
# Occurrence and abundance in two shades of orange, used in all figures;
# supported changes / effects (posterior support >= sig_level) in full colour,
# unsupported ones transparent
col_resp <- c(Occurrence = "#E07B39", Abundance = "#9A3F12")
col_by_response <- c(occurrence = col_resp[["Occurrence"]], log10_abundance = col_resp[["Abundance"]])

# Species ordered by the mean elevation of the plots where they occur (lowest
# at the bottom, highest at the top), the same order in the occurrence and
# abundance panels of a layer
delta_panel <- function(s, resp, xlab, x_labels = waiver(), show_species = TRUE, dat = delta) {
  d <- dat %>% filter(stratum == s, response == resp)
  ord <- dat %>% filter(stratum == s, response == "occurrence") %>%
    arrange(mean_elevation) %>% pull(species)
  d <- d %>% mutate(species   = factor(sp_label(species), levels = sp_label(ord)),
                    supported = factor(direction != "not supported", levels = c(TRUE, FALSE),
                                       labels = c(paste0("Supported (\u2265 ", sig_level, ")"), "Not supported")))
  ggplot(d, aes(y = species, x = mean, alpha = supported)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.6, colour = col_by_response[[resp]]) +
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 1.4, colour = col_by_response[[resp]]) +
    geom_point(size = 2, colour = col_by_response[[resp]]) +
    scale_alpha_manual(values = c(1, 0.3), name = NULL, drop = FALSE) +
    scale_x_continuous(labels = x_labels) +
    scale_y_discrete(drop = FALSE) +
    labs(x = xlab, y = if (show_species) "Species (low to high mean elevation)" else NULL,
         title = if (show_species) strata_labels[[s]] else NULL) +
    theme_mora() +
    theme(axis.text.y = if (show_species) element_text(face = "italic", colour = "black", size = if (s == "overstory") 10 else 8) else element_blank(),
          legend.position = "top")
}
n_sp <- delta %>% filter(response == "occurrence") %>% count(stratum) %>% deframe()
row_for <- function(s, dat = delta) {
  ggarrange(delta_panel(s, "occurrence", "Change in occurrence probability", dat = dat),
            delta_panel(s, "log10_abundance", paste("Change in", abundance_label, "when present"),
                        x_labels = fold_labels, show_species = FALSE, dat = dat),
            ncol = 2, widths = c(1.5, 1), legend = "none", align = "h")
}
strata_present <- intersect(names(models), unique(delta$stratum))
# Relative height of each layer's row: the overstory row (12 species) is given
# more space per species than the understory row, so both are readable
row_heights <- c(overstory = 0.75, understory = 1.6)
fig1 <- ggarrange(plotlist = map(strata_present, row_for), ncol = 1,
                  heights = row_heights[strata_present], labels = LETTERS[seq_along(strata_present)],
                  common.legend = TRUE, legend = "top",
                  legend.grob = get_legend(delta_panel(strata_present[1], "occurrence", "") +
                                             theme(legend.position = "top")))
save_fig(fig1, "fig1_species_change", width = 6.5, height = sum(row_heights[strata_present]) * 4.2 + 0.6)

# Same figure with the change predicted by the full model (fixed effects plus
# plot- and sample-level random effects), averaged over the sampled plots
delta_full <- read_res("delta_species_full.csv")
if (nrow(delta_full) > 0) {
  delta_full <- delta_full %>%
    mutate(direction = factor(direction, levels = direction_levels)) %>%
    left_join(species_elev, by = c("stratum", "species"))
  fig1_full <- ggarrange(plotlist = map(strata_present, ~ row_for(.x, dat = delta_full)), ncol = 1,
                         heights = row_heights[strata_present], labels = LETTERS[seq_along(strata_present)],
                         common.legend = TRUE, legend = "top",
                         legend.grob = get_legend(delta_panel(strata_present[1], "occurrence", "", dat = delta_full) +
                                                    theme(legend.position = "top")))
  save_fig(fig1_full, "fig1b_species_change_full_model", width = 6.5,
           height = sum(row_heights[strata_present]) * 4.2 + 0.6)
  # Occurrence only
  fig1_full_occ <- ggarrange(plotlist = map(strata_present, ~ delta_panel(.x, "occurrence",
                                                                          "Change in occurrence probability",
                                                                          dat = delta_full)),
                             ncol = 1, heights = row_heights[strata_present],
                             labels = LETTERS[seq_along(strata_present)], align = "v",
                             common.legend = TRUE, legend = "top")
  save_fig(fig1_full_occ, "fig1c_occurrence_change_full_model", width = 4.2,
           height = sum(row_heights[strata_present]) * 4.2 + 0.6)
}

# -----------------------------------------------------------------------------
# 2. Change along the elevational gradient
# -----------------------------------------------------------------------------
if (nrow(delta_el) > 0) {
  elev_panel <- function(resp, ylab, y_labels = waiver()) {
    delta_el %>%
      filter(response == resp) %>%
      mutate(stratum_lab = factor(strata_labels[stratum], levels = strata_labels)) %>%
      ggplot(aes(x = elevation, y = mean, group = species, colour = MAT_mean)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      facet_wrap(~ stratum_lab, nrow = 1) +
      scale_colour_viridis_c(option = "plasma", name = "Thermal\naffinity (°C)") +
      scale_y_continuous(labels = y_labels) +
      labs(x = "Elevation (m)", y = ylab) +
      theme_mora()
  }
  p_el_occ <- elev_panel("occurrence", "Change in occurrence probability")
  p_el_abu <- elev_panel("log10_abundance", paste("Change in", abundance_label, "(fold)"), fold_labels)
  p_cti <- cwm_ch %>%
    filter(index == "CTI") %>%
    mutate(stratum_lab = factor(strata_labels[stratum], levels = strata_labels)) %>%
    ggplot(aes(x = elevation, y = rate_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = rate_q05, ymax = rate_q95), fill = col_temp, alpha = 0.25) +
    geom_line(colour = col_temp, linewidth = 0.9) +
    facet_wrap(~ stratum_lab, nrow = 1) +
    labs(x = "Elevation (m)", y = "Predicted CTI change\n(°C per decade)") +
    theme_mora()
  fig2 <- ggarrange(p_el_occ, p_el_abu, p_cti, ncol = 1, labels = c("A", "B", "C"),
                    common.legend = TRUE, legend = "right")
  save_fig(fig2, "fig2_elevation", width = 6.5, height = 7.5)
}

# -----------------------------------------------------------------------------
# 3. Trait effects on the time response (Gamma)
# -----------------------------------------------------------------------------
trait_labels <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity",
                  seed_mass = "Seed mass", wood_density = "Wood density",
                  SLA = "SLA", woodinesswoody = "Woody (vs non-woody)")
# Main-text figure, with larger text than the other figures. One row per layer
# (overstory on top, understory below), occurrence and abundance side by side;
# the x axis of each response covers both layers, so the two rows are directly
# comparable. The layer name is written vertically on the right of each row,
# together with the tree (overstory) or herb (understory) icon of Figures 1
# and 3b; one legend, under the figure.
gamma_fig <- gamma %>%
  mutate(trait = recode(trait, !!!trait_labels),
         trait = factor(trait, levels = rev(unique(trait_labels))),
         response = factor(str_to_sentence(response), levels = c("Occurrence", "Abundance")),
         supported = factor(direction != "not supported", levels = c(TRUE, FALSE),
                            labels = c(paste0("Supported (\u2265 ", sig_level, ")"), "Not supported")))
gamma_range <- gamma_fig %>%
  group_by(response) %>%
  reframe(x = range(c(q05, q95)))

# Direction labels over the top row: Gamma > 0, species with higher trait values
# fared better; Gamma < 0, species with lower trait values fared better
direction_labels <- function(n_traits) {
  y_arrow <- n_traits + 0.75; y_text <- n_traits + 1.05
  list(annotate("segment", x = -0.02, xend = -0.16, y = y_arrow, yend = y_arrow, colour = "grey40",
                linewidth = 0.5, arrow = arrow(length = unit(2, "mm"), type = "closed")),
       annotate("segment", x = 0.02, xend = 0.16, y = y_arrow, yend = y_arrow, colour = "grey40",
                linewidth = 0.5, arrow = arrow(length = unit(2, "mm"), type = "closed")),
       annotate("text", x = -0.02, y = y_text, label = "lower\ntrait values\nfare better", hjust = 1, vjust = 0,
                size = 5.4, fontface = "italic", lineheight = 0.9, colour = "grey30"),
       annotate("text", x = 0.02, y = y_text, label = "higher\ntrait values\nfare better", hjust = 0, vjust = 0,
                size = 5.4, fontface = "italic", lineheight = 0.9, colour = "grey30"))
}

gamma_row <- function(s) {
  d   <- gamma_fig %>% filter(stratum == s) %>% mutate(trait = droplevels(trait))
  top <- s == "overstory"
  n_traits <- nlevels(d$trait)
  ggplot(d, aes(y = trait, x = mean, colour = response, alpha = supported)) +
    geom_blank(data = mutate(gamma_range, trait = factor(levels(d$trait)[1], levels = levels(d$trait))),
               aes(x = x, y = trait), inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.6) +
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 1.5) +
    geom_point(size = 2.6) +
    (if (top) direction_labels(n_traits) else NULL) +
    scale_y_discrete(expand = expansion(add = c(0.6, if (top) 2.6 else 0.6))) +
    coord_cartesian(clip = "off") +   # direction labels may extend past the panel edge
    facet_wrap(~ response, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = col_resp, guide = "none") +
    scale_alpha_manual(values = c(1, 0.25), name = NULL, drop = FALSE,
                       guide = if (top) "none" else "legend") +   # one legend, from the understory row
    # three axis labels per panel: zero and one step on each side
    scale_x_continuous(breaks = function(lim) { step <- if (diff(lim) > 1.2) 0.5 else 0.2; c(-step, 0, step) }) +
    labs(x = if (top) NULL else "Trait effect on change between past and present (\u0393)", y = NULL) +
    theme_mora() +
    theme(axis.text    = element_text(colour = "black", size = 18),
          axis.text.y  = element_text(colour = "black", size = 21),
          axis.title.x = element_text(colour = "black", size = 20),
          legend.text  = element_text(colour = "black", size = 19),
          strip.text   = if (top) element_text(colour = "black", face = "bold", size = 21, hjust = 0) else element_blank(),
          axis.text.x  = if (top) element_blank() else element_text(colour = "black", size = 18),
          axis.ticks.x = if (top) element_blank() else element_line(colour = "black", linewidth = 0.5),
          axis.line.x  = if (top) element_blank() else element_line(colour = "black", linewidth = 0.6))
}

layer_icons <- c(overstory = "\U0001F332", understory = "\U0001F33F")
# Layer name and icon written together vertically, as one label
layer_label <- function(s) {
  wrap_elements(panel = grid::textGrob(paste(layer_icons[[s]], strata_labels[[s]]), rot = -90,
                                       gp = grid::gpar(fontfamily = "Arial Narrow", fontface = "bold", fontsize = 24)))
}

n_traits <- gamma_fig %>% distinct(stratum, trait) %>% count(stratum) %>% deframe()
fig3 <- wrap_plots(gamma_row("overstory"), layer_label("overstory"),
                   gamma_row("understory"), layer_label("understory"),
                   ncol = 2, widths = c(1, 0.06),
                   heights = n_traits[c("overstory", "understory")]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig3_file <- save_fig(fig3, "fig3_traits_gamma", width = 6.8, height = 5.2)

# -----------------------------------------------------------------------------
# 4. Residual associations
# -----------------------------------------------------------------------------
if (nrow(assoc) > 0) {
  assoc_panel <- function(s, resp, lev) {
    d <- assoc %>% filter(stratum == s, response == resp, level == lev)
    lv <- unique(d$sp1)
    d %>%
      mutate(sp1 = factor(sp1, levels = lv), sp2 = factor(sp2, levels = rev(lv))) %>%
      ggplot(aes(x = sp1, y = sp2, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.15) +
      scale_fill_gradient2(low = viridis(1, begin = 0.15), mid = "white", high = col_temp,
                           limits = c(-1, 1), na.value = "grey96", name = "Residual\ncorrelation") +
      labs(x = NULL, y = NULL, title = paste0(strata_labels[[s]], ", ", resp, " (", lev, ")")) +
      theme_mora() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic",
                                       size = if (length(lv) > 30) 4.5 else 7),
            axis.text.y = element_text(face = "italic", size = if (length(lv) > 30) 4.5 else 7),
            axis.line = element_blank())
  }
  combos <- assoc %>% distinct(stratum, response, level) %>%
    arrange(match(stratum, names(models)), response, level)
  fig4 <- ggarrange(plotlist = pmap(combos, ~ assoc_panel(..1, ..2, ..3)),
                    ncol = 2, nrow = ceiling(nrow(combos) / 2),
                    common.legend = TRUE, legend = "right")
  save_fig(fig4, "fig4_residual_associations", width = 9, height = 4.2 * ceiling(nrow(combos) / 2))
  
  # One figure per random level: overstory (top) and understory (bottom),
  # occurrence (left) and abundance (right); species ordered by mean elevation
  # of occurrence; only correlations with posterior support >= sig_level coloured
  assoc_panel_elev <- function(s, resp, lev) {
    d <- assoc %>% filter(stratum == s, response == resp, level == lev)
    if (nrow(d) == 0) return(ggplot() + theme_void())
    ord <- species_elev %>% filter(stratum == s) %>%
      mutate(species = sp_label(species)) %>%
      filter(species %in% d$sp1) %>% arrange(mean_elevation) %>% pull(species)
    ord <- c(ord, setdiff(unique(d$sp1), ord))
    n <- length(ord)
    d %>%
      mutate(sp1 = factor(sp1, levels = ord), sp2 = factor(sp2, levels = ord)) %>%
      ggplot(aes(x = sp1, y = sp2, fill = value)) +
      geom_tile(colour = "white", linewidth = if (n > 30) 0.05 else 0.2) +
      scale_fill_gradient2(low = "#4B3F8C", mid = "white", high = col_temp,
                           limits = c(-1, 1), na.value = "grey95", name = "Residual\ncorrelation") +
      coord_fixed() +
      labs(x = NULL, y = NULL,
           title = paste0(strata_labels[[s]], ": ", if (resp == "occurrence") "occurrence" else "abundance")) +
      theme_mora() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic",
                                       colour = "black", size = if (n > 30) 4.5 else 8),
            axis.text.y = element_text(face = "italic", colour = "black", size = if (n > 30) 4.5 else 8),
            axis.line = element_blank(), axis.ticks = element_blank())
  }
  for (lev in intersect(c("plot", "sample"), unique(assoc$level))) {
    panels <- list(assoc_panel_elev("overstory",  "occurrence", lev),
                   assoc_panel_elev("overstory",  "abundance",  lev),
                   assoc_panel_elev("understory", "occurrence", lev),
                   assoc_panel_elev("understory", "abundance",  lev))
    fig4_lev <- ggarrange(plotlist = panels, ncol = 2, nrow = 2, heights = c(0.6, 1.4),
                          labels = c("A", "B", "C", "D"), common.legend = TRUE, legend = "right")
    save_fig(fig4_lev, paste0("fig4_residual_associations_", lev), width = 9, height = 9)
  }
}

# -----------------------------------------------------------------------------
# 5. Variance partitioning and phylogenetic signal
# -----------------------------------------------------------------------------
if (nrow(vp) > 0) {
  p_vp <- vp %>%
    group_by(stratum, response, component) %>%
    summarise(proportion = mean(proportion), .groups = "drop") %>%
    mutate(stratum_lab = factor(strata_labels[stratum], levels = strata_labels),
           response = str_to_sentence(response)) %>%
    ggplot(aes(x = response, y = proportion, fill = component)) +
    geom_col(width = 0.7) +
    facet_wrap(~ stratum_lab) +
    scale_fill_viridis_d(option = "viridis", name = NULL) +
    scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    labs(x = NULL, y = "Mean proportion of explained variance") +
    theme_mora() + theme(legend.position = "right")
  p_rho <- rho %>%
    mutate(stratum_lab = factor(strata_labels[stratum], levels = strata_labels),
           response = str_to_sentence(response)) %>%
    ggplot(aes(x = response, y = mean, ymin = q05, ymax = q95)) +
    geom_pointrange(colour = "grey20") +
    facet_wrap(~ stratum_lab) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Phylogenetic signal (rho)") +
    theme_mora()
  fig5 <- ggarrange(p_vp, p_rho, ncol = 2, widths = c(1.6, 1), labels = c("A", "B"))
  save_fig(fig5, "fig5_variance_rho", width = 8, height = 3.2)
}

# -----------------------------------------------------------------------------
# 6. Occurrence change against abundance change
# -----------------------------------------------------------------------------
fig6 <- delta %>%
  filter(response %in% c("occurrence", "log10_abundance")) %>%
  select(stratum, stratum_lab, species, response, mean, MAT_mean) %>%
  pivot_wider(names_from = response, values_from = mean) %>%
  filter(!is.na(log10_abundance)) %>%
  ggplot(aes(x = occurrence, y = log10_abundance, colour = MAT_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.2) +
  facet_wrap(~ stratum_lab, scales = "free") +
  scale_colour_viridis_c(option = "plasma", name = "Thermal\naffinity (°C)") +
  scale_y_continuous(labels = fold_labels) +
  labs(x = "Change in occurrence probability",
       y = paste("Change in", abundance_label, "when present")) +
  theme_mora()
save_fig(fig6, "fig6_occurrence_vs_abundance", width = 6.5, height = 3.2)

# -----------------------------------------------------------------------------
# 7. Posterior predictive check (observed vs predicted, past and present)
# -----------------------------------------------------------------------------
# Per species and survey, observed prevalence and mean log10 abundance when
# present against the posterior mean prediction (with all random effects) and
# its 90% interval, from the tables written by 04. One panel per layer and
# response, each with its own axis range.
ppc <- read_res("posterior_predictive_check.csv")
if (nrow(ppc) > 0) {
  ppc_abu_lab <- list(overstory  = c(abu = "basal area", rel = "relative basal area",
                                     abucnt = "stem density", relcnt = "relative density"),
                      understory = c(abu = "cover", rel = "relative cover"))
  resp_code <- str_extract(models[["overstory"]], "[a-z]+$")
  col_time  <- c(past = "#F4A259", present = "#B5501C")
  fig7 <- ppc %>%
    mutate(panel = if_else(response == "occurrence",
                           paste0(strata_labels[stratum], ": occurrence (prevalence)"),
                           paste0(strata_labels[stratum], ": log10 ",
                                  map_chr(stratum, ~ ppc_abu_lab[[.x]][[resp_code]]), " when present")),
           panel = factor(panel, levels = unique(panel[order(match(stratum, names(models)),
                                                             response != "occurrence")])),
           time  = factor(time, levels = c("past", "present"))) %>%
    ggplot(aes(x = observed, y = pred_mean, colour = time)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = pred_q05, ymax = pred_q95), width = 0, alpha = 0.5) +
    geom_point(size = 1.8) +
    facet_wrap(~ panel, ncol = 2, scales = "free") +
    scale_colour_manual(values = col_time, labels = c(past = "Past (~1975)", present = "Present (~2015)"), name = NULL) +
    labs(x = "Observed", y = "Predicted (posterior mean, 90% interval)") +
    theme_mora() +
    theme(legend.position = "top", strip.text = element_text(colour = "black", face = "bold", size = 10, hjust = 0))
  save_fig(fig7, "figS_posterior_predictive_check", width = 6.5, height = 6.5)
}

cat("Figures written to", fig_dir, "\n")

# Open the trait-effects figure (macOS: Preview; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig3_file)
