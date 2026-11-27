# ------------------------------------------------------------------------------
# 04_species_lifeform_contribution_integration.R
# Module: pipeline/02_analysis/01_thermophilization/03_integration_over_under
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: combined figure of contributions to the change in the community
# indices:
#   - Row 1: contribution to CTI change -- overstory by species (left),
#            understory by life form (right)
#   - Row 2: contribution to CMI change -- same split
# Overstory is shown by species (05_species_contribution.R), since all its
# species are trees. Understory is shown by life form
# (07_species_contribution_by_lifeform.R); the species-level understory
# ranking is in 02_understory/figures/.
#
# STYLE (as in 01_thermo_xero_lag_integration.R):
#   - theme_ipsum() with the grid blanked, white panel background, finite
#     black axis lines; "Overstory"/"Understory" shown once as a header row.
#   - Colour encodes the variable: coral for CTI, ocean blue for CMI.
#     Understory bars: full colour = positive contribution, lighter tint =
#     negative. Overstory species: coloured by MAT or MAP affinity with the
#     same palettes and limits as Figure 3 (plasma for MAT, dry-to-wet
#     moisture palette for MAP), with a continuous colour bar. Point shape shows the direction of the mean
#     change in relative basal area (up triangle = increase, down triangle =
#     decrease, circle = little change), because the sign of a contribution
#     depends on both this direction and on whether the species is warmer
#     (wetter) or colder (drier) than the community where it occurs.
#
# A second version groups the understory by woodiness
# (understory_woodiness_contribution.csv) with the same overstory column,
# because life form has a large "Unknown" class while woodiness is almost
# fully resolved.
#
# Output:
#   - figures/species_lifeform_contribution_overstory_understory.png
#   - figures/species_woodiness_contribution_overstory_understory.png
#   - figures/species_contribution_overstory_understory.png (species level in
#     both strata)
#   - figures/species_lifeform_contribution_six_panels.png (overstory species,
#     understory species and understory life forms, for CTI and CMI; reads
#     understory_lifeform/total_contribution_by_process.csv from 07)
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(hrbrthemes)
library(viridisLite)

# 0a. Module directories --------------------------------------------------------
overstory_dir   <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)
understory_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "02_understory", ...)
this_module_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "03_integration_over_under", ...)

dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

N_SPECIES_SHOWN <- 10  # overstory species panel: top and bottom 10 by community contribution

# Mean change in relative basal area below which a species is shown as
# "Little change" (proportion of plot basal area)
CHANGE_THRESHOLD <- 0.01

HEADER_OVERSTORY  <- "\U0001F332 Overstory"
HEADER_UNDERSTORY <- "\U0001F33F Understory"

# Color encodes the VARIABLE, as in 01-03: coral for CTI/MAT (temperature),
# deep ocean blue for CMI/MAP (precipitation). Within a panel, the sign of the
# contribution is shown by shade: the full color for positive, a lighter tint
# of the same color for negative.
THERMO_COLOR <- "#F26B5B"  # coral: temperature / CTI / MAT
XERO_COLOR   <- "#1A759F"  # ocean blue: precipitation / CMI / MAP

# Species affinity palettes, identical to Figure 3
# (02_hmsc/05b_species_change_centered_figure.R) so that a species has the same
# colour in both figures: plasma for thermal affinity over 0-12.5 degC (yellow
# end dropped for readability on white), and the moisture palette from dry
# (orange) to wet (blue) over the MAP range used in Figure 3
COL_MAT    <- viridisLite::plasma(256, end = 0.80, direction = 1)
MAT_LIMITS <- c(0, 12.5)
COL_MAP    <- c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D")
MAP_LIMITS <- c(800, 2700)  # must equal map_lims in 05b_species_change_centered_figure.R

# Mix a color toward white (lighter tint: negative bars, low end of the
# affinity gradient) or toward black (darker shade: high end of the gradient)
lighten_color <- function(col, amount = 0.55) {
  m <- col2rgb(col) * (1 - amount) + 255 * amount
  rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

darken_color <- function(col, amount = 0.25) {
  m <- col2rgb(col) * (1 - amount)
  rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

base_theme <- theme_ipsum(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line.x      = element_line(color = "black", linewidth = 0.3),
    axis.line.y      = element_line(color = "black", linewidth = 0.3),
    # theme_ipsum() sets axis.title.x / axis.title.y as their own elements, which
    # take precedence over the parent axis.title, so both are set directly here
    axis.title.x     = element_text(size = 14, hjust = 0.5),
    axis.title.y     = element_text(size = 14),
    axis.text        = element_text(size = 12),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.text       = element_text(size = 13),
    plot.title        = element_blank()  # no per-panel title -- see the column header below
  )

# Colorbar length is 5 x legend.key.width in ggplot2, so 1.1 cm gives a ~5.5 cm
# bar: long enough that the tick labels (e.g. 1000, 1500, 2000, 2500 mm) do not
# overlap, short enough to stay within the panel width
species_theme  <- base_theme + theme(legend.key.width  = unit(1.1, "cm"),
                                     legend.key.height = unit(0.4, "cm"),
                                     legend.text       = element_text(size = 11),
                                     legend.title      = element_text(size = 12))
lifeform_theme <- base_theme

# 1. Load -------------------------------------------------------------------------
require_file <- function(path) {
  if (!file.exists(path)) stop("Expected file not found: ", path,
                               "\n(check that the relevant stratum's own 05_species_contribution.R or ",
                               "07_species_contribution_by_lifeform.R has been run).")
  path
}

overstory_species_CTI <- read_csv(require_file(overstory_dir("outputs", "overstory_species_contribution_CTI.csv")), show_col_types = FALSE)
overstory_species_CMI <- read_csv(require_file(overstory_dir("outputs", "overstory_species_contribution_CMI.csv")), show_col_types = FALSE)
understory_lifeform   <- read_csv(require_file(understory_dir("outputs", "understory_lifeform_contribution.csv")), show_col_types = FALSE)
understory_species_CTI <- read_csv(require_file(understory_dir("outputs", "understory_species_contribution_CTI.csv")), show_col_types = FALSE)
understory_species_CMI <- read_csv(require_file(understory_dir("outputs", "understory_species_contribution_CMI.csv")), show_col_types = FALSE)
understory_lifeform_process <- read_csv(require_file(understory_dir("outputs", "understory_lifeform_contribution_by_process.csv")), show_col_types = FALSE)
understory_total_process    <- read_csv(require_file(understory_dir("outputs", "understory_total_contribution_by_process.csv")), show_col_types = FALSE) %>%
  rename(life_form = group)

# Life forms represented by very few species and plots are not shown (sedges:
# one species in five plots); their contribution is negligible and is still
# included in the total change in the index
LIFEFORMS_EXCLUDED <- c("Sedge")
understory_lifeform           <- filter(understory_lifeform, !life_form %in% LIFEFORMS_EXCLUDED)
understory_lifeform_process   <- filter(understory_lifeform_process, !life_form %in% LIFEFORMS_EXCLUDED)
understory_woodiness  <- read_csv(require_file(understory_dir("outputs", "understory_woodiness_contribution.csv")), show_col_types = FALSE) %>%
  mutate(woodiness = recode(woodiness, woody = "Woody", `non-woody` = "Non-woody"))

# 2. Panel-building functions -----------------------------------------------------
build_species_panel <- function(df, affinity_col, palette, limits, legend_name, x_lab, n_show = N_SPECIES_SHOWN) {
  df_show <- bind_rows(
    slice_max(df, order_by = community_contribution, n = n_show),
    slice_min(df, order_by = community_contribution, n = n_show)
  ) %>%
    distinct(species, .keep_all = TRUE)
  
  # Direction of the mean change in relative basal area across plots
  df_show %>%
    mutate(
      species = fct_reorder(species, community_contribution),
      change  = case_when(abs(mean_rel_weight_change) < CHANGE_THRESHOLD ~ "Little change",
                          mean_rel_weight_change > 0 ~ "Increased",
                          TRUE ~ "Decreased"),
      change  = factor(change, levels = c("Increased", "Decreased", "Little change"))
    ) %>%
    ggplot(aes(x = community_contribution, y = species,
               color = .data[[affinity_col]], fill = .data[[affinity_col]])) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = community_contribution, y = species, yend = species), linewidth = 0.6) +
    geom_point(aes(shape = change), size = 2.8) +
    scale_shape_manual(values = c(Increased = 24, Decreased = 25, `Little change` = 21),
                       drop = FALSE, name = "Relative basal area",
                       guide = guide_legend(title.position = "top",
                                            override.aes = list(fill = "grey40", colour = "grey40"))) +
    # Affinity palette shared with Figure 3, applied to segment, point outline and fill
    scale_colour_gradientn(colours = palette, limits = limits, oob = scales::squish,
                           name = legend_name, aesthetics = c("colour", "fill"),
                           guide = guide_colorbar(direction = "horizontal", title.position = "top")) +
    labs(x = x_lab, y = NULL) +
    species_theme +
    theme(legend.box = "vertical")
}

build_lifeform_panel <- function(df, group_col, contribution_col, x_lab, base_color) {
  df %>%
    mutate(
      .group    = fct_reorder(.data[[group_col]], .data[[contribution_col]]),
      direction = if_else(.data[[contribution_col]] > 0, "Positive", "Negative")
    ) %>%
    ggplot(aes(x = .data[[contribution_col]], y = .group, fill = direction)) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
    # show.legend = TRUE with limits + drop = FALSE draws both keys, filled,
    # even when a panel has only one sign (e.g. all-positive CTI contributions)
    geom_col(width = 0.6, show.legend = TRUE) +
    scale_fill_manual(values = c("Positive" = base_color, "Negative" = lighten_color(base_color)),
                      limits = c("Negative", "Positive"), drop = FALSE, name = NULL) +
    labs(x = x_lab, y = NULL) +
    lifeform_theme
}

# 3. Overstory panels (species level; identical in both figure versions) ----------
p_CTI_over  <- build_species_panel(overstory_species_CTI, "MAT_mean", COL_MAT, MAT_LIMITS, "MAT affinity (°C)", "Contribution to mean CTI change (°C)")
p_CMI_over  <- build_species_panel(overstory_species_CMI, "MAP_mean", COL_MAP, MAP_LIMITS, "MAP affinity (mm)", "Contribution to mean CMI change (mm)")

# 4. Column header (Overstory / Understory, shown once, bigger) -----------------
column_headers <- wrap_elements(grid::textGrob(HEADER_OVERSTORY,  gp = grid::gpar(fontsize = 20, fontface = "bold"))) |
  wrap_elements(grid::textGrob(HEADER_UNDERSTORY, gp = grid::gpar(fontsize = 20, fontface = "bold")))

# 5. Two figure versions: understory grouped by life form, and by woodiness ----
#    (woodiness is the better-covered trait, see 07_species_contribution_by_lifeform.R)
make_figure <- function(understory_df, group_col, file_name) {
  p_CTI_under <- build_lifeform_panel(understory_df, group_col, "community_contribution_CTI", "Contribution to mean CTI change (°C)", THERMO_COLOR)
  p_CMI_under <- build_lifeform_panel(understory_df, group_col, "community_contribution_CMI", "Contribution to mean CMI change (mm)", XERO_COLOR)
  
  p_grid <- (p_CTI_over | p_CTI_under) / (p_CMI_over | p_CMI_under)
  
  p_final <- column_headers / p_grid +
    plot_layout(heights = c(0.07, 1))
  print(p_final)
  
  ggsave(this_module_dir("figures", file_name), p_final,
         width = 13, height = 10.5, dpi = 300, bg = "white")
  invisible(p_final)
}

p_by_lifeform  <- make_figure(understory_lifeform,  "life_form", "species_lifeform_contribution_overstory_understory.png")
p_by_woodiness <- make_figure(understory_woodiness, "woodiness", "species_woodiness_contribution_overstory_understory.png")

# 5b. Species-level version for both layers ------------------------------------
#     Understory species ranked as in the overstory (top and bottom
#     N_SPECIES_SHOWN by community contribution), with the same affinity
#     palettes and the same shapes for the direction of the mean change in
#     relative cover. The understory mean change in relative cover is shifted
#     upward for species that keep cover in the present, because fewer species
#     had cover per plot in the present survey (Table S3.1), so the shapes are
#     read as a relative indication of direction.
p_CTI_under_sp <- build_species_panel(understory_species_CTI, "MAT_mean", COL_MAT, MAT_LIMITS, "MAT affinity (°C)", "Contribution to mean CTI change (°C)") +
  scale_shape_manual(values = c(Increased = 24, Decreased = 25, `Little change` = 21),
                     drop = FALSE, name = "Relative cover",
                     guide = guide_legend(title.position = "top",
                                          override.aes = list(fill = "grey40", colour = "grey40"))) +
  # Same affinity scale as the overstory panel above, so its colour bar is not repeated
  guides(colour = "none", fill = "none")
p_CMI_under_sp <- build_species_panel(understory_species_CMI, "MAP_mean", COL_MAP, MAP_LIMITS, "MAP affinity (mm)", "Contribution to mean CMI change (mm)") +
  scale_shape_manual(values = c(Increased = 24, Decreased = 25, `Little change` = 21),
                     drop = FALSE, name = "Relative cover",
                     guide = guide_legend(title.position = "top",
                                          override.aes = list(fill = "grey40", colour = "grey40"))) +
  guides(colour = "none", fill = "none")

p_species_both <- column_headers /
  ((p_CTI_over | p_CTI_under_sp) / (p_CMI_over | p_CMI_under_sp)) +
  plot_layout(heights = c(0.06, 1))
print(p_species_both)
ggsave(this_module_dir("figures", "species_contribution_overstory_understory.png"), p_species_both,
       width = 14, height = 13, dpi = 300, bg = "white")

# 5c. Six-panel version: thermal (CTI) and moisture (CMI) contributions in
#     columns; overstory species, understory species and understory life forms
#     in rows. Life-form bars are split into the four processes that shift an
#     index (gains and losses of species warmer/colder or wetter/drier than the
#     community; 07_species_contribution_by_lifeform.R): the gross
#     contribution of each process is shown as wide semi-transparent bars and
#     the net contribution as a narrow dark bar over them. Processes are told
#     apart by hue, not by shade, because shade distinguishes overstory from
#     understory in the other figures: warm hues
#     (CTI) or blues (CMI) for the processes that raise the index, cool hues
#     (CTI) or browns (CMI) for those that lower it -----------------------------
PROCESS_COLORS <- list(
  CTI = c("Gain of warm-affiliated" = "#F89540", "Loss of cold-affiliated" = "#CC4778",
          "Gain of cold-affiliated" = "#7E03A8", "Loss of warm-affiliated" = "#0D0887"),
  CMI = c("Gain of wet-affiliated"  = "#1A759F", "Loss of dry-affiliated"  = "#5AB4AC",
          "Gain of dry-affiliated"  = "#C98B2B", "Loss of wet-affiliated"  = "#8C510A")
)

build_process_panel <- function(df_process, idx, x_lab) {
  df <- df_process %>% filter(index == idx)
  net <- df %>% group_by(life_form) %>% summarise(net = sum(community_contribution), .groups = "drop")
  # Whole understory on top, life forms below ranked by their net contribution
  # to the CTI change; the same order is used in the CTI and CMI panels so
  # that the two panels can be read row by row
  ord <- df_process %>% filter(index == "CTI") %>% group_by(life_form) %>%
    summarise(net = sum(community_contribution), .groups = "drop")
  group_order <- c(ord %>% filter(life_form != "All understory") %>% arrange(net) %>% pull(life_form),
                   intersect("All understory", ord$life_form))
  group_order <- c(setdiff(net$life_form, group_order), group_order)
  df  <- mutate(df,  life_form = factor(life_form, levels = group_order),
                process   = factor(process, levels = names(PROCESS_COLORS[[idx]])))
  net <- mutate(net, life_form = factor(life_form, levels = group_order))
  ggplot(df, aes(x = community_contribution, y = life_form, fill = process)) +
    # Gross contribution of each process: wide, semi-transparent bars; positive
    # and negative parts stack on their own side of zero
    geom_col(width = 0.75, alpha = 0.55) +
    # Net contribution: narrow solid bar drawn over the gross bars of the same group
    geom_col(data = net, aes(x = net, y = life_form, alpha = "Net contribution"), inherit.aes = FALSE,
             width = 0.3, fill = "grey15") +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = PROCESS_COLORS[[idx]], drop = FALSE, name = NULL,
                      guide = guide_legend(nrow = 2, byrow = TRUE, order = 2,
                                           override.aes = list(alpha = 0.55))) +
    scale_alpha_manual(values = c("Net contribution" = 1), name = NULL,
                       guide = guide_legend(order = 1, override.aes = list(fill = "grey15"))) +
    labs(x = x_lab, y = NULL) +
    lifeform_theme +
    theme(legend.text = element_text(size = 11), legend.box = "vertical",
          legend.spacing.y = unit(0.1, "cm"))
}

# The whole-understory bar sets the axis range (its gross contributions are
# several times larger than any life form's), so it is shown only when
# SHOW_UNDERSTORY_TOTAL is TRUE; its values are reported in the text
SHOW_UNDERSTORY_TOTAL <- FALSE
understory_process_all <- if (SHOW_UNDERSTORY_TOTAL) {
  bind_rows(understory_total_process, understory_lifeform_process)
} else {
  understory_lifeform_process
}
p_CTI_under_lf <- build_process_panel(understory_process_all, "CTI", "Contribution to mean CTI change (°C)")
p_CMI_under_lf <- build_process_panel(understory_process_all, "CMI", "Contribution to mean CMI change (mm)")

header_text <- function(label, size = 20, rot = 0) {
  wrap_elements(grid::textGrob(label, rot = rot, gp = grid::gpar(fontsize = size, fontface = "bold")))
}

# Each row is a separate patchwork, so the column widths are set in every row:
# a narrow column for the row labels, next to the panels
row_widths <- plot_layout(widths = c(0.13, 1, 1))
col_header <- (plot_spacer() | header_text("Thermal affinity (CTI)") | header_text("Moisture affinity (CMI)")) + row_widths
# Row labels on two lines so that the shortest row (life forms) is not clipped
row_over   <- (header_text(paste0(HEADER_OVERSTORY, "\nspecies"), 17, 90)     | p_CTI_over     | p_CMI_over)     + row_widths
row_sp     <- (header_text(paste0(HEADER_UNDERSTORY, "\nspecies"), 17, 90)    | p_CTI_under_sp | p_CMI_under_sp) + row_widths
row_lf     <- (header_text(paste0(HEADER_UNDERSTORY, "\nlife forms"), 17, 90) | p_CTI_under_lf | p_CMI_under_lf) + row_widths

# Row heights roughly proportional to the number of bars (15 overstory
# species, up to 20 understory species, 5 life forms) plus the legends
p_six <- (col_header / row_over / row_sp / row_lf) +
  plot_layout(heights = c(0.07, 1, 1.25, 0.75), widths = c(0.05, 1, 1)) +
  plot_annotation(tag_levels = list(c("", "", "", "A", "B", "", "C", "D", "", "E", "F")))
print(p_six)
ggsave(this_module_dir("figures", "species_lifeform_contribution_six_panels.png"), p_six,
       width = 15, height = 20, dpi = 300, bg = "white")

# 6. Summary --------------------------------------------------------------------
cat("\n========================================\n")
cat("SPECIES (OVERSTORY) / LIFE FORM + WOODINESS (UNDERSTORY) CONTRIBUTION INTEGRATION -- COMPLETE\n")
cat("========================================\n")
cat("Saved:\n")
cat(" - figures/species_lifeform_contribution_overstory_understory.png\n")
cat(" - figures/species_woodiness_contribution_overstory_understory.png\n")
cat(" - figures/species_contribution_overstory_understory.png\n")
cat(" - figures/species_lifeform_contribution_six_panels.png\n")

# 7. Open the saved figure --------------------------------------------------------
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(this_module_dir("figures", "species_lifeform_contribution_six_panels.png"))