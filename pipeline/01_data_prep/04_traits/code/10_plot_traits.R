# ------------------------------------------------------------------------------
# 10_plot_traits.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Phylogeny-based trait figures (the trait tables, coverage and validation
# are in 09_trait_imputation.R). Run after 08_integration.R, which completes
# overstory_traits.csv/understory_traits.csv (literature values and reversal
# of the phylogenetic imputation). Independent of 09_trait_imputation.R.
#
# Sections:
#   1. Phylogeny tree plot, tips colored by layer.
#   2. Phylogeny + trait rings per layer (circular ggtree layout with one
#      gheatmap() ring and one legend per trait; see below).
#   3. Species-by-elevation trait plots: species ordered by mean occurrence
#      elevation, trait value +/- SD (from raw TRY records) on x, faceted
#      by trait, per layer.
#
# Figures -> figures/ (created if missing).
#
# Figure conventions: ggsave(), PNG, bg = "white", dpi = 300, theme_ipsum()
# without grid, viridis scales wherever colour/fill maps to data.
# mora_theme() (section 3) adds finite black x/y axis lines;
# mora_tree_theme() (ggtree plots, sections 1-2) omits axis lines, text and
# ticks, since ggtree manages its own coordinates.
#
# Traits shown (longevity, wood density and growth form are not shown; see
# PLOT_NUMERIC_TRAITS/PLOT_CATEGORICAL_TRAITS below):
#   numeric:     seed_mass, SLA
#   categorical: woodiness
#
# Notes:
#   - plot_info.csv (01_cleaning/outputs/) has a plot-ID column and an
#     elevation column; their names are detected from a candidate list and
#     printed at run time (PLOT_KEY_COL/ELEV_COL below).
#   - The offset and width of the heatmap rings (section 2) are set
#     relative to the tree's x range, with RING_WIDTH / RING_START_OFFSET /
#     RING_GAP below.
#
# Section 2 layout: a circular tree layout uses the full circumference for
# the tip labels (120 understory species), which leaves more room per label
# than a rectangular layout. ggtree sets the panel limits from the tree
# geometry alone, so the limits are expanded (hexpand()/vexpand()) to keep
# the labels inside the panel.
#
# Section 2 rings and legends:
#   - Each trait has its own heatmap ring and legend (gheatmap() +
#     new_scale_fill() once per trait). Rings carry no column names
#     (gheatmap(..., colnames = FALSE)), so that they do not overlap the tip
#     labels; legend titles identify the traits and legends are placed at
#     the top of the panel.
#   - Climatic affinity (MAT_mean, MAP_mean from
#     05_species_distribution/02_integrate_gbif_chelsa/data_clean/
#     species_climatic_affinity.csv) and mean occurrence elevation are shown
#     in a separate tree.
#   - seed_mass and SLA are shown as z-scores of the log value
#     (z = scale(log(x))): both are positive, right-skewed traits, and
#     z-scoring the raw values would let one extreme species (e.g. the seed
#     mass of Corylus cornuta) dominate the colour range. This matches the
#     log scale of the genus-imputation validation (r_squared_log,
#     09_trait_imputation.R).
#   - No panel titles; the legends identify the content of each panel.
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggtree",     quietly = TRUE)) BiocManager::install("ggtree", update = FALSE, ask = FALSE)
if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale")
if (!requireNamespace("hrbrthemes", quietly = TRUE)) install.packages("hrbrthemes")
if (!requireNamespace("ragg",       quietly = TRUE)) install.packages("ragg")

library(tidyverse)
library(ape)
library(ggtree)
library(ggnewscale)
library(hrbrthemes)
library(ragg)
library(here)

this_module_dir     <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)
cleaning_dir        <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
outputs_dir         <- function(...) this_module_dir("outputs", ...)
figures_dir         <- function(...) this_module_dir("figures", ...)
# Climatic-affinity module, for the MAT/MAP rings of section 2 (produced
# by 01_combine_gbif_chelsa.R).
climaffinity_dir    <- function(...) here("pipeline", "01_data_prep", "05_species_distribution",
                                          "02_integrate_gbif_chelsa", "data_clean", ...)

dir.create(figures_dir(), recursive = TRUE, showWarnings = FALSE)

VIRIDIS_OPTION <- "D"
QUALITY_DPI    <- 300
GGSAVE_SCALE   <- 1.8

# Base font size, as in 09_trait_imputation.R.
BASE_FONT_SIZE <- 20

mora_theme <- function(base_size = BASE_FONT_SIZE) {
  theme_ipsum(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title    = element_text(size = rel(1.3), face = "bold"),
      plot.subtitle = element_text(size = rel(1.05)),
      axis.title    = element_text(size = rel(1.05)),
      axis.text     = element_text(size = rel(1)),
      strip.text    = element_text(size = rel(1.05), face = "bold"),
      legend.text   = element_text(size = rel(1)),
      legend.title  = element_text(size = rel(1.05))
    )
}

# Variant for the ggtree plots (sections 1-2): same grid removal, white
# background and larger text, without axis lines/text/ticks, since ggtree
# manages its own coordinates and a numeric axis would cross the tree.
mora_tree_theme <- function(base_size = BASE_FONT_SIZE) {
  theme_ipsum(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line  = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title   = element_text(size = rel(1.4), face = "bold"),
      # Larger legend text, title and keys, in proportion to the tip labels
      # and branch widths of the tree figures.
      legend.text   = element_text(size = rel(1.7)),
      legend.title  = element_text(size = rel(1.8), face = "bold"),
      legend.key.size = unit(1.5, "lines")
    )
}

with_retry <- function(expr_fn, label, tries = 5, wait_sec = 2) {
  for (i in seq_len(tries)) {
    result <- try(expr_fn(), silent = TRUE)
    if (!inherits(result, "try-error")) return(invisible(result))
    if (i < tries) {
      cat("  (", label, ": write failed, attempt", i, "of", tries, "-- retrying in", wait_sec, "s)\n", sep = "")
      Sys.sleep(wait_sec)
    }
  }
  stop(label, ": could not write after ", tries, " attempts -- check the file isn't open ",
       "elsewhere and that Polybox isn't stuck syncing that folder.\n",
       "Underlying error: ", attr(result, "condition")$message)
}

save_png <- function(plot, path, width, height, bg = "white", scale = GGSAVE_SCALE) {
  with_retry(function() {
    ggsave(filename = path, plot = plot, device = ragg::agg_png, bg = bg,
           width = width, height = height, dpi = QUALITY_DPI, scale = scale,
           units = "in", limitsize = FALSE)
  }, label = basename(path))
}

# Display trait set for section 3 (elevation plots) and as the categorical
# default for section 2.
PLOT_NUMERIC_TRAITS     <- c("seed_mass", "SLA")
PLOT_CATEGORICAL_TRAITS <- c("woodiness")
# Section 2 functional-trait rings: seed_mass and SLA (longevity is not
# used in the analyses). Woodiness (categorical) is added for the
# understory only; see plot_functional_tree() below.
SECTION2_NUMERIC_TRAITS <- PLOT_NUMERIC_TRAITS

load(outputs_dir("MORA_combined_phylo_tree.rda"))  # loads C, phy
cat("Tree tips:", length(phy$tip.label), "\n")

overstory_traits_raw  <- read_csv(outputs_dir("overstory_traits.csv"),  show_col_types = FALSE)
understory_traits_raw <- read_csv(outputs_dir("understory_traits.csv"), show_col_types = FALSE)

drop_zero_occurrence_species <- function(long_df) {
  long_df %>%
    group_by(species) %>%
    summarise(total_occ = sum(occurrence, na.rm = TRUE), .groups = "drop") %>%
    filter(total_occ > 0) %>%
    pull(species)
}

overstory_long  <- read_csv(cleaning_dir("outputs", "overstory_long.csv"),  show_col_types = FALSE)
understory_long <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE)
overstory_sp  <- drop_zero_occurrence_species(overstory_long)
understory_sp <- drop_zero_occurrence_species(understory_long)

# ==============================================================================
# 1. Phylogeny tree plot, tips colored by layer
# ==============================================================================
cat("\n--- Section 1: phylogeny tree plot ---\n")

stratum_lookup <- tibble(species = phy$tip.label) %>%
  mutate(
    in_over  = species %in% overstory_sp,
    in_under = species %in% understory_sp,
    stratum  = case_when(
      in_over & in_under  ~ "both",
      in_over             ~ "overstory",
      in_under            ~ "understory",
      TRUE                ~ NA_character_
    )
  ) %>%
  select(label = species, stratum)

# Branch width (size = 1 on ggtree(); the default ~0.5 is a hairline at
# this number of tips), tip label size 3.2 and a canvas large enough for
# labels not to overlap.
p_tree <- ggtree(phy, layout = "circular", size = 1) %<+% stratum_lookup +
  geom_tiplab(aes(color = stratum), size = 3.2, fontface = "italic", align = TRUE, offset = 0.5) +
  scale_color_viridis_d(option = VIRIDIS_OPTION, name = "Stratum", na.translate = FALSE) +
  ggtitle(paste0("MORA combined phylogeny (", length(phy$tip.label), " species)")) +
  mora_tree_theme() +
  theme(legend.position = "right")

save_png(p_tree, figures_dir("MORA_combined_phylo_tree.png"), width = 26, height = 26)
cat("Saved: figures/MORA_combined_phylo_tree.png\n")

# ==============================================================================
# 2. Phylogeny + trait rings, per layer (functional-trait tree and
#    climatic-affinity tree)
# ==============================================================================
cat("\n--- Section 2: combined phylogeny + trait heatmap (overstory | understory) ---\n")

# Climatic affinity (MAT_mean deg C, MAP_mean mm), produced by
# 01_combine_gbif_chelsa.R.
species_climatic_affinity <- read_csv(climaffinity_dir("species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)
cat("Climatic affinity loaded for", nrow(species_climatic_affinity), "species\n")

# Ring widths: a circular ring of a given width fraction covers more
# area than a linear column, since it is swept around the circumference.
RING_WIDTH        <- 0.07  # fraction of tree diameter, per single ring
RING_START_OFFSET <- 0.08  # fraction of tree diameter, gap from tips to the first ring
RING_GAP          <- 0.015 # fraction of tree diameter, gap between consecutive rings

# Tip label and heatmap font sizes (label clipping is handled by the panel
# expansion in the tree builder below).
TIPLAB_SIZE    <- 9
BRANCH_LINEWIDTH <- 1.4   # ggtree default is 0.5
# HEATMAP_SCALE = 1, as in 09_trait_imputation.R: ggsave()'s `scale`
# shrinks every plot element as it increases, so text and lines are sized
# directly and the canvas width/height are set to fit them.
HEATMAP_SCALE <- 1

# Mean occurrence elevation per species (climatic-affinity tree):
# elevation is a fixed plot attribute (the past-survey plot_info.csv value
# is used for both periods), averaged across
# all plots x periods where the species was recorded (occurrence == 1),
# pooling both surveys.
plot_info <- read_csv(cleaning_dir("outputs", "plot_info.csv"), show_col_types = FALSE)
plot_elevation <- plot_info %>%
  filter(period == "past", !is.na(elevation)) %>%
  distinct(plot_id, elevation)

mean_occurrence_elevation <- function(long_df) {
  long_df %>%
    filter(occurrence == 1) %>%
    inner_join(plot_elevation, by = "plot_id") %>%
    group_by(species) %>%
    summarise(mean_elevation = mean(elevation, na.rm = TRUE), .groups = "drop")
}

overstory_mean_elevation  <- mean_occurrence_elevation(overstory_long)
understory_mean_elevation <- mean_occurrence_elevation(understory_long)

# Ring widths: a circular ring of a given width fraction covers more
# area than a linear column, since it is swept around the circumference.
RING_WIDTH        <- 0.07  # fraction of tree diameter, per single ring
RING_START_OFFSET <- 0.08  # fraction of tree diameter, gap from tips to the first ring
RING_GAP          <- 0.015 # fraction of tree diameter, gap between consecutive rings

# Tip label and heatmap font sizes (label clipping is handled by the panel
# expansion below).
TIPLAB_SIZE    <- 9
BRANCH_LINEWIDTH <- 1.4   # ggtree default is 0.5
# HEATMAP_SCALE = 1, as in 09_trait_imputation.R: ggsave()'s `scale`
# shrinks every plot element as it increases, so text and lines are sized
# directly and the canvas width/height are set to fit them.
HEATMAP_SCALE <- 1

# Horizontal colour bar used for every continuous ring below; with
# legend.box = "horizontal" in theme(), the legends are laid out side by
# side. Bars of 10 x 0.7 cm make the gradient readable, and legend text and
# title sizes are set here (12/13 pt) rather than inherited from
# mora_tree_theme() (legend.text = rel(1.7)), which is too large for the
# tick labels of a horizontal bar.
horizontal_colorbar <- guide_colorbar(direction = "horizontal", title.position = "top",
                                      barwidth = unit(10, "cm"), barheight = unit(0.7, "cm"))

# Ring palettes, in three families: grey (white-black) for elevation, a
# neutral physical quantity; viridis "inferno" (black -> red -> yellow)
# for thermal affinity and SLA; viridis "mako" (dark blue/green -> pale)
# for moisture affinity and seed mass. Two breaks per legend
# (scales::breaks_pretty), the low and high ends, are enough to read a
# gradient legend.
continuous_fill_scale <- function(palette, name, na_value = "grey85") {
  breaks_fn <- scales::breaks_pretty(n = 2)
  switch(palette,
         grey      = scale_fill_gradient(low = "white", high = "grey10", name = name,
                                         na.value = na_value, guide = horizontal_colorbar, breaks = breaks_fn),
         redyellow = scale_fill_viridis_c(option = "inferno", name = name,
                                          na.value = na_value, guide = horizontal_colorbar, breaks = breaks_fn),
         greenblue = scale_fill_viridis_c(option = "mako", name = name,
                                          na.value = na_value, guide = horizontal_colorbar, breaks = breaks_fn),
         stop("Unknown palette: ", palette)
  )
}

# Ring-tree builder, used for the two trees per layer (functional traits
# and climatic affinity), each with its own set of rings.
#
# `ring_specs` is a list of per-ring specs, each one of:
#   list(values = <1-col data frame, rownames = species>, title = <legend
#        title>, discrete = FALSE, option = <viridis letter>,
#        direction = 1/-1, na_value = "grey85")
#   list(values = ..., title = ..., discrete = TRUE, option = <viridis
#        letter>)
# applied in order, outward from the tree tips.
build_ring_tree <- function(stratum_tree, ring_specs) {
  tip_order <- stratum_tree$tip.label
  
  p_current <- ggtree(stratum_tree, layout = "circular", size = BRANCH_LINEWIDTH, open.angle = 15)
  tree_width <- max(p_current$data$x, na.rm = TRUE)
  # Units: gheatmap()'s `offset` is in data units (the scale of the tree's
  # x axis), whereas `width` is a fraction of the tree width, applied
  # internally, so it is passed as the plain fraction (RING_WIDTH, e.g.
  # 0.07), not multiplied by tree_width. ring_offset stays in data units
  # and is advanced by each ring's rendered width (tree_width * RING_WIDTH).
  ring_offset    <- tree_width * RING_START_OFFSET
  ring_width_abs <- tree_width * RING_WIDTH  # for ring_offset bookkeeping only
  ring_gap       <- tree_width * RING_GAP
  
  for (spec in ring_specs) {
    ring_df <- spec$values[tip_order, , drop = FALSE]
    colnames(ring_df) <- spec$title
    p_current <- (p_current + new_scale_fill()) %>%
      gheatmap(ring_df, offset = ring_offset, width = RING_WIDTH, colnames = FALSE, color = NULL)
    p_current <- if (isTRUE(spec$discrete)) {
      p_current + scale_fill_viridis_d(option = spec$option, name = spec$title, na.translate = FALSE,
                                       guide = guide_legend(direction = "horizontal", title.position = "top"))
    } else {
      p_current + continuous_fill_scale(spec$palette, spec$title, spec$na_value)
    }
    ring_offset <- ring_offset + ring_width_abs + ring_gap
  }
  
  outer_edge <- ring_offset - ring_gap  # last ring's own outer edge
  
  # Tip labels added LAST, after every ring, with geom_tiplab2() -- the
  # ggtree label layer meant for adding labels after ring layers like
  # gheatmap(), since gheatmap() returns a coordinate-transformed plot that
  # the ordinary geom_tiplab() does not handle correctly.
  TIPLAB_OFFSET <- outer_edge + tree_width * 0.02
  p_final <- p_current +
    geom_tiplab2(size = TIPLAB_SIZE, fontface = "italic", align = TRUE, offset = TIPLAB_OFFSET)
  
  # Panel-limit expansion: ggtree/gheatmap set the x/y limits from the tree
  # and ring geometry alone, before the labels of geom_tiplab2() are added
  # beyond them, so labels would be clipped at the panel edge.
  # hexpand()/vexpand() extend the limits in data units on all four sides;
  # applied last, after the labels, so that the expansion is sized against
  # the final extent of the plot.
  PANEL_EXPAND_FRACTION <- 0.6
  p_final <- p_final +
    hexpand(PANEL_EXPAND_FRACTION, direction = 1)  + hexpand(PANEL_EXPAND_FRACTION, direction = -1) +
    vexpand(PANEL_EXPAND_FRACTION, direction = 1)  + vexpand(PANEL_EXPAND_FRACTION, direction = -1)
  
  p_final +
    mora_tree_theme() +
    theme(legend.position = "top", legend.box = "horizontal",
          legend.box.just = "top", legend.spacing.x = unit(1, "cm"),
          legend.key.width = unit(1.4, "cm"),
          # Explicit sizes, overriding mora_tree_theme()'s legend.text/
          # legend.title (rel(1.7)/rel(1.8)), which are too large for a bar
          # this wide.
          legend.text  = element_text(size = 15),
          legend.title = element_text(size = 17, face = "bold"),
          # Every axis.* element is blanked explicitly, because ggtree draws a
          # radial branch-length scale in the circular layout that the
          # shorthand axis.text/axis.ticks do not remove.
          axis.text.x  = element_blank(), axis.text.y  = element_blank(),
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.ticks.length = unit(0, "pt"))
}

# Shared setup (tree pruning, tip order, log-z trait rings) used by both
# the functional and the climatic tree for a given layer.
prep_stratum_tree <- function(traits, stratum_label) {
  sp_on_tree <- intersect(traits$species, phy$tip.label)
  n_dropped <- nrow(traits) - length(sp_on_tree)
  if (n_dropped > 0) cat(stratum_label, "-- dropping", n_dropped, "species absent from the phylogeny\n")
  
  stratum_tree <- keep.tip(phy, sp_on_tree)
  list(tree = stratum_tree, sp_on_tree = sp_on_tree, n_tips = length(stratum_tree$tip.label))
}

# Functional-trait tree (right-hand panel): seed mass and SLA (z-score of
# the log value; see header), plus woodiness for the understory only
# (overstory species are all woody trees, so the ring would carry no
# information there).
TRAIT_DISPLAY_SPECS <- list(
  seed_mass = list(title = "Seed mass", palette = "greenblue"),
  SLA       = list(title = "SLA",       palette = "redyellow")
)

plot_functional_tree <- function(traits, stratum_label) {
  numeric_present <- intersect(SECTION2_NUMERIC_TRAITS, names(traits))
  categorical_present <- if (stratum_label == "overstory") character(0) else intersect(PLOT_CATEGORICAL_TRAITS, names(traits))
  if (length(numeric_present) == 0) { cat("No numeric traits for", stratum_label, "-- skipping.\n"); return(list(plot = NULL, n_tips = 0)) }
  
  prepped <- prep_stratum_tree(traits, stratum_label)
  sp_on_tree <- prepped$sp_on_tree
  
  numeric_rings <- traits %>%
    filter(species %in% sp_on_tree) %>%
    select(species, all_of(numeric_present)) %>%
    column_to_rownames("species") %>%
    mutate(across(everything(), ~ as.numeric(scale(log(.x)))))
  
  ring_specs <- map(numeric_present, function(trait_name) {
    spec <- TRAIT_DISPLAY_SPECS[[trait_name]]
    list(values = numeric_rings[, trait_name, drop = FALSE], title = spec$title,
         discrete = FALSE, palette = spec$palette, na_value = "grey85")
  })
  
  if (length(categorical_present) > 0) {
    categorical_values <- traits %>%
      filter(species %in% sp_on_tree) %>%
      select(species, all_of(categorical_present)) %>%
      column_to_rownames("species") %>%
      mutate(across(everything(), as.factor))
    ring_specs <- c(ring_specs, list(list(values = categorical_values, title = str_to_sentence(categorical_present),
                                          discrete = TRUE, option = VIRIDIS_OPTION)))
  }
  
  list(plot = build_ring_tree(prepped$tree, ring_specs), n_tips = prepped$n_tips)
}

# Climatic-affinity tree (left-hand panel): MAT and MAP affinity and mean
# occurrence elevation, shown as plain z-scores (scale(x)). Unlike
# seed_mass/SLA they are not log-transformed, since MAT_mean can be at or
# below zero for cold-affinity species. Species without an affinity or
# elevation are shown as grey (NA) cells.
CLIMAFFINITY_SPECS <- list(
  MAT_mean       = list(title = "Thermal aff. (z)",  palette = "redyellow"),
  MAP_mean       = list(title = "Moisture aff. (z)", palette = "greenblue"),
  mean_elevation = list(title = "Mean elev. (z)",     palette = "grey")
)

plot_climatic_tree <- function(traits, stratum_label, mean_elevation_df) {
  prepped <- prep_stratum_tree(traits, stratum_label)
  sp_on_tree <- prepped$sp_on_tree
  
  climaffinity_values <- species_climatic_affinity %>%
    filter(species %in% sp_on_tree) %>%
    full_join(mean_elevation_df %>% filter(species %in% sp_on_tree), by = "species") %>%
    column_to_rownames("species")
  n_no_affinity  <- sum(is.na(climaffinity_values$MAT_mean) | is.na(climaffinity_values$MAP_mean))
  n_no_elevation <- sum(is.na(climaffinity_values$mean_elevation))
  if (n_no_affinity  > 0) cat(" ", stratum_label, "--", n_no_affinity,  "species with no computed climatic affinity (blank in the MAT/MAP rings)\n")
  if (n_no_elevation > 0) cat(" ", stratum_label, "--", n_no_elevation, "species with no computed mean occurrence elevation (blank in that ring)\n")
  
  climaffinity_z <- climaffinity_values %>%
    mutate(across(everything(), ~ as.numeric(scale(.x))))
  
  ring_specs <- map(names(CLIMAFFINITY_SPECS), function(var) {
    spec <- CLIMAFFINITY_SPECS[[var]]
    list(values = climaffinity_z[, var, drop = FALSE], title = spec$title,
         discrete = FALSE, palette = spec$palette, na_value = "grey85")
  })
  
  list(plot = build_ring_tree(prepped$tree, ring_specs), n_tips = prepped$n_tips)
}

if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
library(ggpubr)

# One figure per layer, with the climatic-affinity tree on the left and the
# functional-trait tree on the right.
save_stratum_figure <- function(traits, stratum_label, mean_elevation_df) {
  climatic_result   <- plot_climatic_tree(traits, stratum_label, mean_elevation_df)
  functional_result <- plot_functional_tree(traits, stratum_label)
  
  p_combined <- ggarrange(climatic_result$plot, functional_result$plot, ncol = 2, widths = c(1, 1))
  
  out_path <- figures_dir(str_glue("{stratum_label}_trait_phylogeny_heatmap.png"))
  # A circular tree needs an equal-aspect (square) panel, or the circle
  # renders as an ellipse -- each of the two panels is saved as a
  # PANEL_SIZE x PANEL_SIZE square, so the combined canvas is 2*PANEL_SIZE
  # wide and PANEL_SIZE tall.
  n_tips <- max(climatic_result$n_tips, functional_result$n_tips)
  PANEL_SIZE <- max(24, n_tips * 0.26)
  save_png(p_combined, out_path, width = PANEL_SIZE * 2, height = PANEL_SIZE, scale = HEATMAP_SCALE)
  cat("Saved:", out_path, "\n")
}

save_stratum_figure(overstory_traits_raw,  "overstory",  overstory_mean_elevation)
save_stratum_figure(understory_traits_raw, "understory", understory_mean_elevation)

# ==============================================================================
# 3. Traits by elevation, per layer
# ==============================================================================
cat("\n--- Section 3: traits by elevation ---\n")

plot_info <- read_csv(cleaning_dir("outputs", "plot_info.csv"), show_col_types = FALSE)
cat("plot_info.csv columns:", paste(names(plot_info), collapse = ", "), "\n")
PLOT_KEY_COL <- intersect(c("plot", "plot_id", "Plot", "PlotID"), names(plot_info))[1]
ELEV_COL     <- intersect(c("elevation", "elevation_m", "elev", "elev_m"), names(plot_info))[1]
if (is.na(PLOT_KEY_COL) || is.na(ELEV_COL)) {
  stop("Couldn't find a plot-ID / elevation column in plot_info.csv -- set PLOT_KEY_COL/ELEV_COL explicitly (see columns printed above).")
}
cat("Using plot key column:", PLOT_KEY_COL, "| elevation column:", ELEV_COL, "\n")

plot_elev <- plot_info %>%
  transmute(plot = .data[[PLOT_KEY_COL]], elevation = as.numeric(.data[[ELEV_COL]])) %>%
  filter(!is.na(elevation)) %>%
  distinct(plot, .keep_all = TRUE)

compute_species_elevation <- function(long_df) {
  plot_key <- intersect(c("plot", "plot_id", "Plot", "PlotID"), names(long_df))[1]
  if (is.na(plot_key)) stop("Couldn't find a plot-ID column in this long table.")
  long_df %>%
    filter(occurrence > 0) %>%
    rename(plot = all_of(plot_key)) %>%
    left_join(plot_elev, by = "plot") %>%
    group_by(species) %>%
    summarise(mean_elevation = mean(elevation, na.rm = TRUE), n_plots = n(), .groups = "drop")
}

overstory_elev  <- compute_species_elevation(overstory_long)
understory_elev <- compute_species_elevation(understory_long)

raw_obj <- readRDS(this_module_dir("tmp_files", "trait_matching_intermediate.rds"))
list2env(raw_obj, envir = environment())

numeric_raw_tables <- list(seed_mass = seed_mass_raw, SLA = sla_raw)

summarise_raw_trait <- function(raw_df, trait_label) {
  raw_df %>%
    group_by(species) %>%
    summarise(trait_mean = mean(value, na.rm = TRUE),
              trait_sd   = if (n() >= 2) sd(value, na.rm = TRUE) else NA_real_,
              n_records  = n(), .groups = "drop") %>%
    mutate(trait = trait_label)
}
raw_trait_long <- map2_dfr(numeric_raw_tables, names(numeric_raw_tables), summarise_raw_trait)

load_final_traits <- function(traits) {
  present_traits <- intersect(PLOT_NUMERIC_TRAITS, names(traits))
  traits %>% select(species, all_of(present_traits)) %>%
    pivot_longer(cols = -species, names_to = "trait", values_to = "final_value")
}

build_plot_data <- function(final_traits, species_elev, stratum_label) {
  final_traits %>%
    left_join(raw_trait_long, by = c("species", "trait")) %>%
    mutate(trait_mean = coalesce(trait_mean, final_value), n_records = coalesce(n_records, 0L)) %>%
    filter(!is.na(trait_mean)) %>%
    inner_join(species_elev %>% select(species, mean_elevation), by = "species") %>%
    filter(!is.na(mean_elevation)) %>%
    mutate(stratum = stratum_label)
}

overstory_plot_data  <- build_plot_data(load_final_traits(overstory_traits_raw),  overstory_elev,  "overstory")
understory_plot_data <- build_plot_data(load_final_traits(understory_traits_raw), understory_elev, "understory")

# seed_mass facet, display-only cap at 90 mg: species above it (Corylus
# cornuta, 826 mg, a single TRY record) are left out of this panel so that
# the x axis is not dominated by one extreme value. The value is kept in
# understory_traits.csv and in all tables and models.
SEED_MASS_DISPLAY_CAP <- 90
understory_plot_data <- understory_plot_data %>%
  filter(!(trait == "seed_mass" & trait_mean > SEED_MASS_DISPLAY_CAP))

plot_traits_by_elevation <- function(plot_data, stratum_label, out_path) {
  species_order <- plot_data %>% distinct(species, mean_elevation) %>% arrange(mean_elevation) %>% pull(species)
  plot_data <- plot_data %>% mutate(species = factor(species, levels = species_order))
  n_species <- length(species_order)
  
  # Plain black, no per-trait color -- trait identity is already given by
  # the facet label.
  p <- ggplot(plot_data, aes(x = trait_mean, y = species)) +
    geom_pointrange(aes(xmin = trait_mean - trait_sd, xmax = trait_mean + trait_sd),
                    size = 0.25, linewidth = 0.4, na.rm = TRUE, color = "black") +
    facet_wrap(~ trait, scales = "free_x", nrow = 1) +
    labs(title = paste0(str_to_title(stratum_label), ": trait values by species, ordered by mean occurrence elevation"),
         subtitle = "Species ordered bottom (low elevation) to top (high elevation); error bars = +/- 1 SD across raw TRY records (point only where < 2 records)",
         x = "Trait value", y = NULL) +
    mora_theme() +
    theme(axis.text.y = element_text(face = "italic", size = rel(0.75)),
          strip.text = element_text(face = "bold"))
  
  save_png(p, out_path, width = 11, height = max(6, n_species * 0.14))
}

plot_traits_by_elevation(overstory_plot_data,  "overstory",  figures_dir("overstory_traits_by_elevation.png"))
plot_traits_by_elevation(understory_plot_data, "understory", figures_dir("understory_traits_by_elevation.png"))
cat("Saved: figures/overstory_traits_by_elevation.png, figures/understory_traits_by_elevation.png\n")

cat("\n========================================\n")
cat("10_plot_traits.R COMPLETE\n")
cat("Figures -> figures/MORA_combined_phylo_tree.png, figures/combined_trait_phylogeny_heatmap.png,\n")
cat("           figures/<stratum>_traits_by_elevation.png\n")
cat("========================================\n")
