# ------------------------------------------------------------------------------
# 03_cti_cmi_past_present_elevation_integration.R
# Module: pipeline/02_analysis/01_thermophilization/03_integration_over_under
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: combine the overstory and understory "index vs elevation, past vs
# present" figures (each produced separately by
# 01_overstory/code/03_CTI_CMI_past_present_elevation.R and
# 02_understory/code/03_CTI_CMI_past_present_elevation.R) into one 2x2 panel:
#   - Row 1: CTI -- overstory (left), understory (right)
#   - Row 2: CMI -- overstory (left), understory (right)
# Reads each layer's own CTI_CMI_quantiles.csv (from 01_compute_cti_cmi*.R);
# does not recompute anything.
#
# STYLE (as in 01_thermo_xero_lag_integration.R):
#   - theme_ipsum() with panel.grid.major/minor blanked and a white panel
#     background; finite black x/y axis lines.
#   - Color encodes the VARIABLE, as in 01/02/04: coral for the CTI row (MAT,
#     temperature), deep ocean blue for the CMI row (MAP, precipitation), the
#     same in the overstory and understory column. Period is shown by shade
#     within the row's color: lighter tint = past (~1975), full color =
#     present (~2015). A neutral grey light/dark legend explains the period
#     encoding for all four panels.
#   - "Overstory"/"Understory" shown once, bigger/bolder, as a single
#     header row above the whole 2x2 grid (not repeated as a title on every
#     panel), with the same optional tree/understory pictogram used in
#     01_thermo_xero_lag_integration.R.
#
# PRIMARY WEIGHT PER LAYER: same as 01_thermo_xero_lag_integration.R --
# basal_area_spp_perc for overstory, cover for understory.
#
# Output:
#   - figures/CTI_CMI_past_present_elevation_overstory_understory.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(hrbrthemes)
library(cowplot)
library(viridisLite)

# 0a. Module directories --------------------------------------------------------
overstory_dir   <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)
understory_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "02_understory", ...)
this_module_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "03_integration_over_under", ...)

dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

OVERSTORY_WEIGHT  <- "basal_area_spp_perc"
UNDERSTORY_WEIGHT <- "cover"

# Coral for temperature (CTI) and deep ocean blue for precipitation (CMI), same
# two colors as in every script of this module. Past = lighter tint of the
# row's color, present = the full color.
THERMO_COLOR <- "#F26B5B"  # coral: temperature / CTI
XERO_COLOR   <- "#1A759F"  # ocean blue: precipitation / CMI

lighten_color <- function(col, amount = 0.55) {
  m <- col2rgb(col) * (1 - amount) + 255 * amount
  rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

PERIOD_PAST    <- "Past (~1975)"
PERIOD_PRESENT <- "Present (~2015)"

HEADER_OVERSTORY  <- "\U0001F332 Overstory"
HEADER_UNDERSTORY <- "\U0001F33F Understory"

finite_axis_theme <- theme_ipsum(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line.x      = element_line(color = "black", linewidth = 0.3),
    axis.line.y      = element_line(color = "black", linewidth = 0.3),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 15),
    axis.text        = element_text(size = 13),
    legend.position  = "none",
    plot.title       = element_blank()
  )

# 1. Load each layer's CTI/CMI quantiles, restrict to primary weight --------
require_file <- function(path) {
  if (!file.exists(path)) stop("Expected file not found: ", path,
                               "\n(check that this stratum's own 01_compute_CTI_CMI.R has been run).")
  path
}

overstory_cq <- read_csv(require_file(overstory_dir("data_clean", "overstory_CTI_CMI_quantiles.csv")), show_col_types = FALSE) %>%
  filter(weight_var == OVERSTORY_WEIGHT) %>%
  mutate(period = factor(period, levels = c("past", "present"), labels = c("Past (~1975)", "Present (~2015)")))

understory_cq <- read_csv(require_file(understory_dir("data_clean", "understory_CTI_CMI_quantiles.csv")), show_col_types = FALSE) %>%
  filter(weight_var == UNDERSTORY_WEIGHT) %>%
  mutate(period = factor(period, levels = c("past", "present"), labels = c("Past (~1975)", "Present (~2015)")))

# 2. Panel-building function ------------------------------------------------------
build_index_panel <- function(df, median_col, y_lab, base_color) {
  period_colors <- c(lighten_color(base_color), base_color)
  names(period_colors) <- c(PERIOD_PAST, PERIOD_PRESENT)
  ggplot(df, aes(x = elevation, y = .data[[median_col]], color = period, fill = period)) +
    geom_point(alpha = 0.35, size = 1.4) +
    geom_smooth(method = "loess", span = 0.75, se = TRUE, linewidth = 1, alpha = 0.15) +
    scale_color_manual(values = period_colors, name = NULL) +
    scale_fill_manual(values = period_colors, name = NULL) +
    labs(x = NULL, y = y_lab) +
    finite_axis_theme
}

# 3. Build the four panels ---------------------------------------------------------
p_CTI_over  <- build_index_panel(overstory_cq,  "CTI_median", "Community temperature index (°C)", THERMO_COLOR)
p_CTI_under <- build_index_panel(understory_cq, "CTI_median", "Community temperature index (°C)", THERMO_COLOR)
p_CMI_over  <- build_index_panel(overstory_cq,  "CMI_median", "Community moisture index (mm)", XERO_COLOR)
p_CMI_under <- build_index_panel(understory_cq, "CMI_median", "Community moisture index (mm)", XERO_COLOR)

p_CTI_over
p_CTI_under
p_CMI_over
p_CMI_under

# 4. Shared legend (Past vs Present), extracted once -----------------------------
# Neutral grey light/dark pair, since the period encoding applies to both the
# coral and the blue rows
legend_df <- tibble(
  x      = rep(c(1, 2), 2),
  y      = rep(c(1, 2), 2),
  period = factor(rep(c(PERIOD_PAST, PERIOD_PRESENT), each = 2), levels = c(PERIOD_PAST, PERIOD_PRESENT))
)
p_legend_source <- ggplot(legend_df, aes(x = x, y = y, color = period)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = setNames(c(lighten_color("grey35"), "grey35"), c(PERIOD_PAST, PERIOD_PRESENT)), name = NULL) +
  theme_ipsum(base_size = 15) +
  theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 15))
shared_legend <- get_legend(p_legend_source)

# 5. Column header (Overstory / Understory, shown once, bigger) + grid + shared
#    x-axis label, all stacked together ------------------------------------------
column_headers <- wrap_elements(grid::textGrob(HEADER_OVERSTORY,  gp = grid::gpar(fontsize = 20, fontface = "bold"))) |
  wrap_elements(grid::textGrob(HEADER_UNDERSTORY, gp = grid::gpar(fontsize = 20, fontface = "bold")))

p_grid <- (p_CTI_over | p_CTI_under) / (p_CMI_over | p_CMI_under)

x_axis_label <- wrap_elements(grid::textGrob("Elevation (m)", gp = grid::gpar(fontsize = 17)))

p_index_integration <- wrap_elements(shared_legend) / column_headers / p_grid / x_axis_label +
  plot_layout(heights = c(0.06, 0.07, 1, 0.04))
p_index_integration

ggsave(this_module_dir("figures", "CTI_CMI_past_present_elevation_overstory_understory.png"), p_index_integration,
       width = 12, height = 11.5, dpi = 300, bg = "white")

# 6. Summary --------------------------------------------------------------------
cat("\n========================================\n")
cat("CTI/CMI PAST-PRESENT-ELEVATION INTEGRATION (OVERSTORY x UNDERSTORY) -- COMPLETE\n")
cat("========================================\n")
cat("Saved:\n")
cat(" - figures/CTI_CMI_past_present_elevation_overstory_understory.png\n")
