# =============================================================================
# 05b_species_change_centered_figure.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the local computer, after
#           23_species_change_centered.R; copy outputs/centered/
#           species_change_centered.csv from the computing machine first)
# Author: Agostina Torres
#
# Species change between past and present, with understory abundance shown
# relative to the mean change across species (centered within each posterior
# sample, 23_species_change_centered.R),
# in the style of the species-change figure of 05_hmsc_figures_combined.R:
# change in occurrence probability (full model) and fold change in relative
# abundance when present, per species, overstory and understory.
#
#   - Species are ordered by the mean elevation of the plots where they occur
#     (lowest at the bottom), the same order in the occurrence and abundance
#     panels of a layer (as in 05).
#   - Species names are coloured by thermal affinity (MAT_mean; plasma, 0–12.5 °C),
#     from dark blue (cold-affiliated) to orange (warm-affiliated).
#   - A mountain profile behind the species names illustrates the elevational
#     order; the panel has no axes (row positions carry no values).
#   - The x range is shared by overstory and understory within each response,
#     so the zero lines are aligned between layers; the x axis is drawn in the
#     understory row only (centering of understory abundance stated in the caption).
#   - Supported changes (posterior support >= 0.90 in one direction) in full
#     colour, unsupported ones transparent (as in 05).
#   - Moisture affinity (MAP_mean, mm) uses the moisture palette from orange
#     (dry-affiliated) to blue (wet-affiliated). Two layouts are written:
#     names on both sides coloured by thermal (left) and moisture (right)
#     affinity, or black names with a right-hand panel of both affinities.
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/centered/species_change_centered.csv (23)
#          pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/
#            data_clean/species_climatic_affinity.csv
#          pipeline/01_data_prep/01_cleaning/outputs/plot_info.csv, <stratum>_long.csv
# Outputs: pipeline/02_analysis/02_hmsc/figures/m1_final/
#            ms1__fig1d_species_change_centered_m1_final_names.png
#            ms1__fig1d_species_change_centered_m1_final_panel.png
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)
library(patchwork)

# -----------------------------------------------------------------------------
# 0. Settings and style (as in 05)
# -----------------------------------------------------------------------------
tag        <- "m1_final"
fig_prefix <- "ms1__"
sig_level  <- 0.90
strata_labels <- c(overstory = "Overstory", understory = "Understory")
# Panel titles with the layer symbols used in the thermophilization figures
strata_titles <- c(overstory = "\U0001F332 Overstory", understory = "\U0001F33F Understory")
col_resp <- c(occurrence = "#E07B39", log10_abundance = "#9A3F12")
# Thermal affinity of species names: plasma palette over the range of the
# thermal-affinity figures, from dark blue (low MAT, cold-affiliated) to orange
# (high MAT, warm-affiliated); the yellow end is dropped so that all names are
# readable on white
col_affinity <- viridisLite::plasma(256, end = 0.80, direction = 1)
mat_limits   <- c(0, 12.5)
# Moisture affinity: palette of the moisture-affinity figures, dry (orange) to wet (blue)
col_moisture <- c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D")
row_heights  <- c(overstory = 0.75, understory = 1.6)
# How climatic affinities are shown (one figure is written per option):
#   "names": species names on both sides, coloured by thermal affinity (left)
#            and by moisture affinity (right)
#   "panel": species names in black and a right-hand panel with thermal and
#            moisture affinity as dots on their own axes
affinity_displays <- c("names", "panel")
# Colour keys (MAT, MAP) and support legend under the figure (TRUE), or the
# colour keys inside the overstory panels (FALSE)
keys_at_bottom <- TRUE

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))

theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = 0.6),
          axis.line.y = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          axis.text = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 11),
          axis.title.y = element_text(colour = "black", size = 11),
          legend.text = element_text(colour = "black"), legend.title = element_text(colour = "black"),
          plot.title = element_text(colour = "black", size = 11, face = "bold"),
          plot.margin = margin(12, 10, 10, 10))
}
fold_labels <- function(x) paste0("×", signif(10^x, 2))
sp_label    <- function(x) str_replace_all(x, "_", " ")

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
fig_dir    <- file.path(module_dir, "figures", tag)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Inputs
# -----------------------------------------------------------------------------
# Only understory abundance is centered: the change shared by all species there
# (common increase in relative cover) comes from the change in sampling design.
# Understory occurrence is already nearly centered by the reduction of the past
# data, and in the overstory, where the design did not change, the common change
# is interpretable; these panels show the uncentered change.
centered_panels <- tibble(stratum = "understory", response = "log10_abundance")
centered <- read_csv(file.path(module_dir, "outputs", "centered", "species_change_centered.csv"),
                     show_col_types = FALSE) %>%
  mutate(response = recode(response, abundance = "log10_abundance"),
         species  = str_replace_all(species, " ", "_")) %>%
  left_join(mutate(centered_panels, use_centered = TRUE), by = c("stratum", "response")) %>%
  filter(version == if_else(coalesce(use_centered, FALSE), "centered", "raw")) %>%
  select(-use_centered)
is_centered <- function(s, resp) any(centered_panels$stratum == s & centered_panels$response == resp)

affinity <- read_csv(here("pipeline", "01_data_prep", "05_species_distribution",
                          "02_integrate_gbif_chelsa", "data_clean", "species_climatic_affinity.csv"),
                     show_col_types = FALSE) %>%
  transmute(species = str_replace_all(species, " ", "_"), MAT_mean, MAP_mean)

# Mean elevation of the plots where each species occurs (either time), as in 05
cleaning_dir <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
plot_info <- read_csv(cleaning_dir("outputs", "plot_info.csv"), show_col_types = FALSE)
plot_key  <- intersect(c("plot_id", "plot", "Plot", "PlotID"), names(plot_info))[1]
elev_col  <- intersect(c("elevation", "elevation_m", "elev", "elev_m"), names(plot_info))[1]
plot_elev <- plot_info %>%
  transmute(plot_id = as.character(.data[[plot_key]]), elevation = as.numeric(.data[[elev_col]])) %>%
  filter(!is.na(elevation)) %>% distinct(plot_id, .keep_all = TRUE)
species_elev <- map_dfr(names(strata_labels), function(s) {
  read_csv(cleaning_dir("outputs", paste0(s, "_long.csv")), show_col_types = FALSE) %>%
    filter(occurrence > 0) %>%
    mutate(plot_id = as.character(plot_id)) %>%
    left_join(plot_elev, by = "plot_id") %>%
    group_by(species) %>%
    summarise(mean_elevation = mean(elevation, na.rm = TRUE), .groups = "drop") %>%
    mutate(species = str_replace_all(species, " ", "_"), stratum = s)
})

centered <- centered %>%
  left_join(species_elev, by = c("stratum", "species")) %>%
  left_join(affinity, by = "species") %>%
  mutate(supported = factor(direction != "not supported", levels = c(TRUE, FALSE),
                            labels = c(paste0("Supported (≥ ", sig_level, ")"), "Not supported")))
if (any(is.na(centered$mean_elevation)) || any(is.na(centered$MAT_mean))) {
  warning("Missing elevation or thermal affinity for: ",
          paste(unique(centered$species[is.na(centered$mean_elevation) | is.na(centered$MAT_mean)]), collapse = ", "))
}

# Shared x range per response, symmetric around zero (zero aligned between layers and centred in each panel)
x_lims <- centered %>% group_by(response) %>%
  summarise(lo = min(q05, na.rm = TRUE), hi = max(q95, na.rm = TRUE), .groups = "drop") %>%
  { setNames(purrr::map2(.$lo, .$hi, ~ c(-1, 1) * max(abs(c(.x, .y))) * 1.04), .$response) }   # symmetric, zero in the middle

# Shared range of moisture affinity (both strata), for the x axis and colours
map_lims <- range(centered$MAP_mean, na.rm = TRUE)
map_lims <- map_lims + c(-0.04, 0.04) * diff(map_lims)

# -----------------------------------------------------------------------------
# 2. Mountain profile (unit square), drawn behind the species names
# -----------------------------------------------------------------------------
# Ridge line from control points with random midpoint displacement (fixed seed),
# so the profile looks like a real silhouette; the left flank is cut by the
# panel edge. A few short strokes suggest rock faces.
displace <- function(x, y, iters, sd, decay = 0.55) {
  for (i in seq_len(iters)) {
    n  <- length(x)
    mx <- (x[-n] + x[-1]) / 2
    my <- (y[-n] + y[-1]) / 2 + rnorm(n - 1, 0, sd)
    x  <- c(rbind(x[-n], mx), x[n])
    y  <- c(rbind(y[-n], my), y[n])
    sd <- sd * decay
  }
  tibble(x = x, y = y)
}
set.seed(7)
ridge <- displace(c(0, 0.16, 0.30, 0.42, 0.52, 0.60, 0.70, 0.80, 0.90, 1.00),
                  c(0.58, 0.66, 0.80, 0.90, 1.00, 0.88, 0.70, 0.48, 0.24, 0.00), iters = 6, sd = 0.045)
strokes <- list(list(x = c(0.52, 0.50, 0.48), y = c(0.99, 0.93, 0.88)),
                list(x = c(0.56, 0.58, 0.57), y = c(0.94, 0.89, 0.84)))
rocks <- imap_dfr(strokes, ~ displace(.x$x, .x$y, iters = 3, sd = 0.015) %>% mutate(g = .y))

# -----------------------------------------------------------------------------
# 3. Panels
# -----------------------------------------------------------------------------
species_order <- function(s) {
  centered %>% filter(stratum == s, response == "occurrence") %>%
    arrange(mean_elevation) %>% pull(species)
}

# Species names (coloured by thermal affinity) over the mountain profile
label_panel <- function(s, coloured = TRUE) {
  ord <- species_order(s); n <- length(ord)
  labs_df <- centered %>% filter(stratum == s, response == "occurrence") %>%
    mutate(y = match(species, ord), label = sp_label(species))
  p <- ggplot() +
    geom_path(data = mutate(filter(ridge, x >= 0.30), y = 0.5 + y * n), aes(x, y), colour = "grey50", alpha = 0.6,
              linewidth = 0.45, lineend = "round", linejoin = "round") +
    geom_path(data = mutate(rocks, y = 0.5 + y * n), aes(x, y, group = g), colour = "grey50", alpha = 0.6,
              linewidth = 0.3) +
    (if (coloured) geom_text(data = labs_df, aes(x = 0.98, y = y, label = label, colour = MAT_mean),
                             hjust = 1, fontface = "italic", size = if (s == "overstory") 4.4 else 3.3)
     else geom_text(data = labs_df, aes(x = 0.98, y = y, label = label), colour = "black",
                    hjust = 1, fontface = "italic", size = if (s == "overstory") 4.4 else 3.3)) +
    scale_colour_gradientn(colours = col_affinity, limits = mat_limits, oob = scales::squish, guide = "none")
  # MAT colour key drawn inside the overstory panel, vertical, at the lower left
  # under the mountain flank (data coordinates, so it is not collected with the
  # support legend)
  if (!keys_at_bottom && s == "overstory" && coloured) {
    key_x <- c(0.33, 0.37); key_y <- c(1, 0.45 * n)
    # 0 °C at the top and 12.5 °C at the bottom
    key <- tibble(y = seq(key_y[1], key_y[2], length.out = 100), val = seq(mat_limits[2], mat_limits[1], length.out = 100))
    ticks <- tibble(val = c(0, 4, 8, 12)) %>%
      mutate(y = key_y[2] - (val - mat_limits[1]) / diff(mat_limits) * diff(key_y))
    p <- p +
      geom_tile(data = key, aes(x = mean(key_x), y = y, fill = val), width = diff(key_x),
                height = diff(key_y) / 99) +
      scale_fill_gradientn(colours = col_affinity, limits = mat_limits, guide = "none") +
      geom_text(data = ticks, aes(x = key_x[2] + 0.015, y = y, label = val), hjust = 0, size = 3.8) +
      annotate("text", x = key_x[1], y = key_y[2] + 0.6, label = "MAT (°C)", hjust = 0, size = 4.2)
  }
  p +
    # Only the right part of the profile is shown (x from 0.30); clipping is off
    # so that long species names are not cut at the panel edge
    coord_cartesian(xlim = c(0.30, 1), ylim = c(0.5, n + 0.5), expand = FALSE, clip = "off") +
    labs(x = NULL, y = NULL, title = strata_titles[[s]]) +
    theme_void() +
    theme(plot.title = element_text(colour = "black", size = 22, face = "bold", margin = margin(b = 14)),
          plot.margin = margin(12, 10, 10, 10))
}

# The x axis is drawn only in the bottom (understory) row: both layers share the
# same x range, so the overstory row needs no axis of its own
# (centering of understory abundance is explained in the caption)
change_panel <- function(s, resp, xlab, x_labels = waiver(), show_x = (s == "understory")) {
  ord <- species_order(s)
  d <- centered %>% filter(stratum == s, response == resp) %>%
    mutate(species = factor(sp_label(species), levels = sp_label(ord)))
  ggplot(d, aes(y = species, x = mean, alpha = supported)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.6, colour = col_resp[[resp]]) +
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 1.4, colour = col_resp[[resp]]) +
    geom_point(size = 2, colour = col_resp[[resp]]) +
    scale_alpha_manual(values = c(1, 0.3), name = NULL, drop = FALSE) +
    scale_x_continuous(limits = x_lims[[resp]], labels = x_labels) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0.5)) +   # rows span 0.5 to n + 0.5, as the label panel
    labs(x = if (show_x) xlab else NULL, y = NULL) +
    theme_mora() +
    theme(axis.text.y = element_blank(), legend.position = "top",
          axis.text.x = element_text(colour = "black", size = 15),
          axis.title.x = element_text(colour = "black", size = 20, hjust = 0.5)) +
    (if (!show_x) theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                        axis.text.x = element_blank()) else NULL) +
    # One support legend for both responses and layers, drawn in neutral grey
    # (occurrence and abundance differ only in the shade of orange)
    guides(alpha = guide_legend(override.aes = list(colour = "grey30")))
}

# Option "names": species names on the right, coloured by moisture affinity,
# with the MAP colour key at the lower right of the overstory panel
right_label_panel <- function(s) {
  ord <- species_order(s); n <- length(ord)
  labs_df <- centered %>% filter(stratum == s, response == "occurrence") %>%
    mutate(y = match(species, ord), label = sp_label(species))
  p <- ggplot() +
    geom_text(data = labs_df, aes(x = 0.02, y = y, label = label, colour = MAP_mean),
              hjust = 0, fontface = "italic", size = if (s == "overstory") 4.4 else 3.3) +
    scale_colour_gradientn(colours = col_moisture, limits = map_lims, guide = "none")
  if (!keys_at_bottom && s == "overstory") {
    key_x <- c(0.86, 0.90); key_y <- c(1, 0.45 * n)
    brks <- scales::breaks_pretty(3)(map_lims); brks <- brks[brks >= map_lims[1] & brks <= map_lims[2]]
    key <- tibble(y = seq(key_y[1], key_y[2], length.out = 100), val = seq(map_lims[1], map_lims[2], length.out = 100))
    ticks <- tibble(val = brks) %>% mutate(y = key_y[1] + (val - map_lims[1]) / diff(map_lims) * diff(key_y))
    p <- p +
      geom_tile(data = key, aes(x = mean(key_x), y = y, fill = val), width = diff(key_x), height = diff(key_y) / 99) +
      scale_fill_gradientn(colours = col_moisture, limits = map_lims, guide = "none") +
      geom_text(data = ticks, aes(x = key_x[2] + 0.015, y = y, label = val), hjust = 0, size = 3.8) +
      annotate("text", x = key_x[1], y = key_y[2] + 0.6, label = "MAP (mm)", hjust = 0, size = 4.2)
  }
  p + coord_cartesian(xlim = c(0, 1), ylim = c(0.5, n + 0.5), expand = FALSE, clip = "off") +
    labs(x = NULL, y = NULL) + theme_void() + theme(plot.margin = margin(12, 45, 10, 10))   # room for the MAP key
}

# Option "panel": thermal and moisture affinity of each species as dots, in two
# facets with their own axes (ranges shared by both strata), coloured with the
# thermal (plasma) and moisture palettes
affinity_panel <- function(s, show_x = (s == "understory")) {
  ord <- species_order(s)
  pal_T <- scales::col_numeric(col_affinity, domain = mat_limits)
  pal_P <- scales::col_numeric(col_moisture, domain = map_lims)
  d <- centered %>% filter(stratum == s, response == "occurrence") %>%
    transmute(species = factor(sp_label(species), levels = sp_label(ord)),
              `MAT (°C)` = MAT_mean, `MAP (mm)` = MAP_mean) %>%
    pivot_longer(-species, names_to = "affinity", values_to = "value") %>%
    mutate(affinity = factor(affinity, levels = c("MAT (°C)", "MAP (mm)"))) %>%
    group_by(affinity) %>%
    mutate(col = if (first(affinity) == "MAT (°C)") pal_T(pmin(pmax(value, mat_limits[1]), mat_limits[2]))
           else pal_P(value)) %>%
    ungroup()
  # invisible points fix each facet's range to the range shared by both layers
  lims <- tibble(affinity = factor(rep(levels(d$affinity), each = 2), levels = levels(d$affinity)),
                 value = c(mat_limits, map_lims), species = factor(sp_label(ord[1]), levels = sp_label(ord)))
  ggplot(d, aes(x = value, y = species)) +
    geom_blank(data = lims) +
    geom_point(aes(colour = col), size = if (s == "overstory") 3 else 2.2) +
    scale_colour_identity() +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0.5)) +
    facet_grid(. ~ affinity, scales = "free_x") +
    labs(x = if (show_x) "Climatic affinity" else NULL, y = NULL) +
    theme_mora() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "black", size = 11),
          axis.title.x = element_text(colour = "black", size = 16, hjust = 0.5),
          strip.text = element_text(colour = "black", size = 11, hjust = 0.5),
          panel.spacing.x = unit(10, "pt")) +
    (if (!show_x) theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                        axis.text.x = element_blank()) else NULL) +
    (if (s == "understory") theme(strip.text = element_blank()) else NULL)
}

# Moisture affinity of each species, in the same row order as the change panels
moisture_panel <- function(s, show_x = (s == "understory")) {
  ord <- species_order(s)
  d <- centered %>% filter(stratum == s, response == "occurrence") %>%
    mutate(species = factor(sp_label(species), levels = sp_label(ord)))
  ggplot(d, aes(y = species, x = MAP_mean, colour = MAP_mean)) +
    geom_segment(aes(x = map_lims[1], xend = MAP_mean, yend = species), colour = "grey85", linewidth = 0.4) +
    geom_point(size = if (s == "overstory") 3 else 2.2) +
    scale_colour_gradientn(colours = col_moisture, limits = map_lims, guide = "none") +
    scale_x_continuous(limits = map_lims, breaks = scales::breaks_pretty(3)) +
    scale_y_discrete(drop = FALSE, expand = c(0, 0.5)) +
    labs(x = if (show_x) "Moisture affinity\n(MAP, mm)" else NULL, y = NULL) +
    theme_mora() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "black", size = 13),
          axis.title.x = element_text(colour = "black", size = 16, hjust = 0.5)) +
    (if (!show_x) theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                        axis.text.x = element_blank()) else NULL)
}

# Horizontal colour keys for the bottom of the figure: thermal affinity under
# the species names on the left, moisture affinity under those on the right
colour_key <- function(cols, limits, breaks, title) {
  ggplot(tibble(x = seq(limits[1], limits[2], length.out = 200)), aes(x, 0, fill = x)) +
    geom_tile(height = 0.6, width = 1.3 * diff(limits) / 199) +
    scale_fill_gradientn(colours = cols, limits = limits, guide = "none") +
    scale_x_continuous(breaks = breaks, expand = c(0, 0)) +
    coord_cartesian(xlim = limits, ylim = c(-0.5, 0.5), clip = "off") +
    labs(title = title) +
    theme_void() +
    theme(axis.text.x = element_text(colour = "black", size = 16, margin = margin(t = 3)),
          plot.title = element_text(colour = "black", size = 17, hjust = 0.5, margin = margin(b = 5)),
          plot.margin = margin(4, 20, 4, 20))
}
# Support key drawn as a plot, so that it sits in the same row as the colour
# keys: thick and thin bars with a point, in full colour (supported) or
# transparent (not supported), as in the change panels
support_key <- function() {
  items <- tibble(x0 = c(0, 3.9), a = c(1, 0.3),
                  label = c(paste0("Supported (\u2265 ", sig_level, ")"), "Not supported"))
  ggplot(items) +
    geom_segment(aes(x = x0 - 0.35, xend = x0 + 0.35, y = 0, yend = 0, alpha = a), colour = "grey30", linewidth = 0.8) +
    geom_segment(aes(x = x0 - 0.17, xend = x0 + 0.17, y = 0, yend = 0, alpha = a), colour = "grey30", linewidth = 1.9) +
    geom_point(aes(x = x0, y = 0, alpha = a), colour = "grey30", size = 2.8) +
    geom_text(aes(x = x0 + 0.5, y = 0, label = label), hjust = 0, size = 17 / .pt) +
    scale_alpha_identity() +
    coord_cartesian(xlim = c(-0.5, 7), ylim = c(-0.5, 0.5), clip = "off") +
    labs(title = " ") +   # same title height as the colour keys, so the rows align
    theme_void() +
    theme(plot.title = element_text(size = 17, margin = margin(b = 5)), plot.margin = margin(4, 20, 4, 20))
}
keys_row <- function() {
  map_brks <- scales::breaks_pretty(3)(map_lims)
  map_brks <- map_brks[map_brks >= map_lims[1] & map_brks <= map_lims[2]]
  colour_key(col_affinity, mat_limits, c(0, 4, 8, 12), "Thermal affinity (MAT, \u00b0C)") +
    support_key() +
    colour_key(col_moisture, map_lims, map_brks, "Moisture affinity (MAP, mm)") +
    plot_layout(widths = c(0.9, 2.5, 0.9))
}

row_for <- function(s, display) {
  occ <- change_panel(s, "occurrence", "Change in occurrence\nprobability")
  abu <- change_panel(s, "log10_abundance", "Change in relative\nabundance when present", x_labels = fold_labels)
  if (display == "names") {
    label_panel(s, coloured = TRUE) + occ + abu + right_label_panel(s) + plot_layout(widths = c(0.9, 1.5, 1, 0.9))
  } else {
    label_panel(s, coloured = FALSE) + occ + abu + affinity_panel(s) + plot_layout(widths = c(0.9, 1.5, 1, 0.9))
  }
}

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
strata_present <- intersect(names(strata_labels), unique(centered$stratum))
for (display in affinity_displays) {
  rows    <- purrr::map(strata_present, row_for, display = display)
  heights <- row_heights[strata_present]
  if (keys_at_bottom && display == "names") {
    rows <- c(rows, list(keys_row())); heights <- c(heights, 0.22)
  }
  fig <- wrap_plots(rows, ncol = 1, heights = heights) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = list(c("A", "", "", "", "B", "", "", ""))) &
    theme(legend.position = if (keys_at_bottom && display == "names") "none" else "top", legend.text = element_text(colour = "black", size = 13),
          plot.tag = element_text(size = 18, face = "bold"))
  fig_file <- file.path(fig_dir, paste0(fig_prefix, "fig1d_species_change_centered_", tag, "_", display, ".png"))
  ggsave(fig_file, fig, device = "png", dpi = 300, width = 10.5,
         height = sum(row_heights[strata_present]) * 4.2 + if (keys_at_bottom) 1.8 else 0.8, scale = 1.3, bg = "white", limitsize = FALSE)
  cat("Figure written to", fig_file, "\n")
  open_fig(fig_file)
}
