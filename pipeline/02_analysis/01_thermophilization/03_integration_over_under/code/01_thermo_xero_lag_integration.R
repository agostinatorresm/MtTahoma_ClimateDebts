# ------------------------------------------------------------------------------
# 01_thermo_xero_lag_integration.R
# Module: pipeline/02_analysis/01_thermophilization/03_integration_over_under
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: combine the overstory and understory "observed vs expected" rate-
# vs-elevation figures (each produced separately by
# 01_overstory/code/02_thermo_xero_elevation_lag.R and
# 02_understory/code/02_thermo_xero_elevation_lag.R) into a single 2x2 panel:
#   - Row 1: thermophilization (CTI) -- overstory (left), understory (right)
#   - Row 2: xerophilization (CMI)   -- overstory (left), understory (right)
# A first column summarises, for each index and layer, the distribution of
# plot-level observed and expected rates (violins) with the mean expected rate,
# and, in a zoom panel with its own axis, the posterior mean observed rate with
# its 90% and 50% credible intervals (Bayesian model A of 07_rate_models_bayes.R),
# because the elevational curves alone do not show the overall rates reported
# in the Results.
# Reads the elevation_curves and thermo_xero_rates files each layer's own
# module 01/02 scripts already produce, and the coefficients table of 07; does
# not recompute anything.
#
# STYLE:
#   - theme_ipsum(), with panel.grid.major/minor blanked and a white panel
#     background (blanking only the parent "panel.grid" element is not fully
#     inherited by theme_ipsum()); finite black x/y axis lines.
#   - Colour encodes the variable (CTI vs CMI), not the layer: coral for
#     thermophilization (temperature), deep ocean blue for xerophilization
#     (precipitation/moisture), the same colour in the overstory and
#     understory columns of each row.
#   - "Overstory"/"Understory" shown once as a header row above the 2x2
#     grid, with a small pictogram in each header (tree for overstory,
#     herb for understory).
#
# The shared legend (built once from a dummy plot with
# cowplot::get_legend() and placed below the grid) has 3 keys:
# thermophilization observed (coral), xerophilization observed (ocean blue),
# expected/space-for-time (black dashed lines in B, grey violins in A).
#
# PRIMARY WEIGHT PER LAYER: overstory uses basal_area_spp_perc (dominance-
# weighted), understory uses cover -- same as each layer's own 02 script.
#
# Inputs (besides the layer files):
#   - <stratum>/data_clean/<stratum>_thermo_xero_expected.csv (plot-level lag)
#   - outputs/rate_models_bayes_coefficients.csv (07_rate_models_bayes.R)
# Output:
#   - figures/thermo_xero_lag_overstory_understory.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(hrbrthemes)
library(cowplot)      # get_legend()

# 0a. Module directories --------------------------------------------------------
overstory_dir   <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)
understory_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "02_understory", ...)
this_module_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "03_integration_over_under", ...)

dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

OVERSTORY_WEIGHT  <- "basal_area_spp_perc"
UNDERSTORY_WEIGHT <- "cover"

# Coral for temperature (thermophilization/CTI) and deep ocean blue for
# precipitation/moisture (CMI/xerophilization) -- same color used for both
# layers within a row; the row (variable), not the column (stratum), is what
# color encodes. The blue sits between teal and blue (blue-dominant), a
# deep-water tone. Same two colors in every script of this module.
THERMO_COLOR <- "#F26B5B"  # coral: temperature / CTI
XERO_COLOR   <- "#1A759F"  # ocean blue: precipitation / CMI

# 0b. Shared theme: theme_ipsum(), grid blanked (major + minor), white
#     panel background, finite black axis lines, no legend per panel (one
#     shared legend is built in section 4) ---------------------------------------
finite_axis_theme <- theme_ipsum(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line.x      = element_line(color = "black", linewidth = 0.3),
    axis.line.y      = element_line(color = "black", linewidth = 0.3),
    axis.title.x     = element_blank(),  # shown once, centered below the grid -- see section 5
    axis.title.y     = element_text(size = 18),
    axis.text        = element_text(size = 17),
    legend.position  = "none",
    plot.title       = element_blank()   # no per-panel title -- see the column header in section 5
  )

# 1. Load each layer's curve + rate files, checking existence explicitly ------
require_file <- function(path) {
  if (!file.exists(path)) stop("Expected file not found: ", path,
                               "\n(check that this stratum's own 01_compute_CTI_CMI.R and ",
                               "02_thermo_xero_elevation_lag.R have both been run).")
  path
}

overstory_curves  <- read_csv(require_file(overstory_dir("outputs", "overstory_elevation_curves.csv")), show_col_types = FALSE)
understory_curves <- read_csv(require_file(understory_dir("outputs", "understory_elevation_curves.csv")), show_col_types = FALSE)

overstory_rates <- read_csv(require_file(overstory_dir("outputs", "overstory_thermo_xero_rates.csv")), show_col_types = FALSE) %>%
  filter(weight_var == OVERSTORY_WEIGHT)
understory_rates <- read_csv(require_file(understory_dir("outputs", "understory_thermo_xero_rates.csv")), show_col_types = FALSE) %>%
  filter(weight_var == UNDERSTORY_WEIGHT)

# Mean observed rate (model A) and mean lag (model A2) per layer and index,
# with credible intervals
mean_rates <- read_csv(require_file(this_module_dir("outputs", "rate_models_bayes_coefficients.csv")), show_col_types = FALSE) %>%
  filter(term == "Intercept", str_detect(model, " A mean rate$| A2 mean lag$")) %>%
  mutate(stratum  = factor(str_to_sentence(word(model, 1)), levels = c("Overstory", "Understory")),
         variable = word(model, 2),
         quantity = factor(if_else(str_detect(model, "A2 mean lag"), "Lag", "Observed rate"),
                           levels = c("Observed rate", "Lag")))
stopifnot("Expected 8 rows (2 strata x 2 indices x rate and lag)" = nrow(mean_rates) == 8,
          "90% and 50% intervals missing: rerun 07_rate_models_bayes.R" = all(c("lower_90", "lower_50") %in% names(mean_rates)))
cat("\nMean rates and lags (90% and 50% CI):\n")
print(mean_rates %>% select(stratum, variable, quantity, mean, lower_90, upper_90, lower_50, upper_50, pd))

# 2. Plot-level rates and panel-building functions -------------------------------
# Plot-level observed and expected rates of each layer, long format. The
# plot-level lag is read from <stratum>_thermo_xero_expected.csv (the same
# values used for model A2 of 07_rate_models_bayes.R; lag = expected -
# observed), and the expected rate is the observed rate plus the lag.
read_lag <- function(dir_fun, st, weight) {
  e <- read_csv(require_file(dir_fun("data_clean", paste0(st, "_thermo_xero_expected.csv"))), show_col_types = FALSE)
  if ("weight_var" %in% names(e)) e <- filter(e, weight_var == weight)
  select(e, plot_id, lag_CTI, lag_CMI)
}
overstory_rates  <- left_join(overstory_rates,  read_lag(overstory_dir,  "overstory",  OVERSTORY_WEIGHT),  by = "plot_id")
understory_rates <- left_join(understory_rates, read_lag(understory_dir, "understory", UNDERSTORY_WEIGHT), by = "plot_id")

plot_rates <- map_dfr(list(Overstory = overstory_rates, Understory = understory_rates), function(df) {
  map_dfr(c("CTI", "CMI"), function(v) {
    tibble(variable = v,
           Observed = df[[paste0("rate_", v)]],
           Expected = df[[paste0("rate_", v)]] + df[[paste0("lag_", v)]])
  })
}, .id = "stratum") %>%
  pivot_longer(c(Observed, Expected), names_to = "type", values_to = "rate") %>%
  filter(!is.na(rate)) %>%
  mutate(stratum = factor(stratum, levels = c("Overstory", "Understory")),
         type    = factor(type, levels = c("Observed", "Expected")))

# Check: mean expected rate across plots vs the mean of the expected curve
curves_all <- bind_rows(Overstory = overstory_curves, Understory = understory_curves, .id = "stratum")
exp_check <- plot_rates %>% filter(type == "Expected") %>%
  group_by(stratum, variable) %>% summarise(mean_expected_plots = mean(rate), .groups = "drop") %>%
  left_join(curves_all %>% group_by(stratum, variable) %>%
              summarise(mean_expected_curve = mean(rate_expected, na.rm = TRUE), .groups = "drop"),
            by = c("stratum", "variable"))
cat("\nMean expected rate, plots vs curve:\n"); print(exp_check)

# Summary panels (first column), in two parts:
#   - violins: distribution of plot-level observed (coloured) and expected
#     (grey) rates per layer; the vertical distance between them is the lag;
#     dark point on the grey violin: mean expected rate across plots.
#     Black point on the observed violin: posterior mean observed rate. A
#     dashed box marks the range enlarged in the zoom panel (arrow).
#   - zoom panel with its own y axis: posterior mean observed rate from the
#     intercept-only model rate ~ 1 (model A of 07_rate_models_bayes.R; the
#     intercept is the mean rate across plots and f its probability of
#     direction, i.e. support for a mean rate different from zero) with 90% (thin) and 50% (thick) credible
#     intervals and its probability of direction f. These intervals are much
#     narrower than the spread among plots, so they are not visible at the
#     scale of the violins.
FUNNEL_X <- 0.8

# Common y limits per index for the violins and the elevation panels (the
# y-axis title is shown once per row, on the violins)
y_limits <- function(variable_name) {
  curves <- bind_rows(overstory_curves, understory_curves) %>% filter(variable == variable_name)
  r <- range(c(filter(plot_rates, variable == variable_name)$rate, curves$rate_lo, curves$rate_hi), na.rm = TRUE)
  r + c(-0.04, 0.04) * diff(r)
}

# y range of the zoom panel (90% intervals and zero, with margins); also drawn
# as a dashed box on the violins to show which part is enlarged
# The zero of the zoom panel is placed at the same height as the zero of the
# violins: the zoom limits keep the same proportion below and above zero as
# the violin limits, and are just wide enough for the 90% intervals and labels.
zoom_range <- function(variable_name) {
  mr <- filter(mean_rates, variable == variable_name, quantity == "Observed rate")
  lo <- min(mr$lower_90, 0); hi <- max(mr$upper_90, 0); r <- hi - lo
  lo <- lo - 0.15 * r; hi <- hi + 0.30 * r
  yl <- y_limits(variable_name)
  p  <- -yl[1] / diff(yl)                      # position of zero in the violin panel
  k  <- max(hi / (1 - p), -lo / p)             # total height of the zoom range
  c(-p * k, (1 - p) * k)
}

build_violin_panel <- function(variable_name, point_color, y_lab, title = NULL) {
  d  <- filter(plot_rates, variable == variable_name)
  mr <- filter(mean_rates, variable == variable_name, quantity == "Observed rate")
  me <- d %>% filter(type == "Expected") %>% group_by(stratum) %>% summarise(mean = mean(rate), .groups = "drop")
  zr <- zoom_range(variable_name)
  yl <- y_limits(variable_name)
  dodge <- position_dodge(width = 0.8)
  ggplot(d, aes(x = stratum, y = rate, fill = type)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_violin(position = dodge, width = 0.75, colour = NA, alpha = 0.5, scale = "width") +
    # posterior mean observed rate (model A) on the observed violin
    geom_point(data = mr, aes(x = as.numeric(stratum) - 0.2, y = mean), inherit.aes = FALSE,
               size = 2.2, colour = "black") +
    # mean expected rate across plots on the expected violin
    geom_point(data = me, aes(x = as.numeric(stratum) + 0.2, y = mean), inherit.aes = FALSE,
               size = 2.2, colour = "grey20") +
    # small box (thin solid grey line) around the two mean observed rates: range enlarged in
    # the zoom panel (on the right), joined to the top and bottom of the zoom
    # panel by two lines (funnel); FUNNEL_X sets how far right the lines reach
    # beyond the panel edge, in x-axis units (adjust if they do not meet the
    # zoom panel)
    annotate("rect", xmin = 0.68, xmax = 1.92, ymin = zr[1], ymax = zr[2],
             fill = NA, colour = "grey50", linewidth = 0.25) +
    annotate("segment", x = 1.92, xend = 2.6 + FUNNEL_X, y = zr[2], yend = yl[2],
             colour = "grey50", linewidth = 0.25) +
    annotate("segment", x = 1.92, xend = 2.6 + FUNNEL_X, y = zr[1], yend = yl[1],
             colour = "grey50", linewidth = 0.25) +
    coord_cartesian(xlim = c(0.4, 2.6), ylim = yl, expand = FALSE, clip = "off") +
    scale_fill_manual(values = c(Observed = point_color, Expected = "grey55"), name = NULL) +
    scale_x_discrete(labels = c(Overstory = "\U0001F332", Understory = "\U0001F33F")) +
    labs(x = NULL, y = y_lab, title = title) +
    finite_axis_theme +
    theme(axis.text.x = element_text(size = 22, colour = "black"), legend.position = "none",
          axis.title.y = element_text(size = 20)) +
    theme(plot.title = element_text(size = 23, face = "bold", hjust = 0, margin = margin(b = 10)),
          plot.title.position = "plot")
}

build_zoom_panel <- function(variable_name, point_color, y_lab, title = NULL) {
  mr <- filter(mean_rates, variable == variable_name, quantity == "Observed rate")
  ggplot(mr, aes(x = stratum, y = mean)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_linerange(aes(ymin = lower_90, ymax = upper_90), colour = point_color,
                   linewidth = 0.9, lineend = "round") +
    geom_linerange(aes(ymin = lower_50, ymax = upper_50), colour = point_color,
                   linewidth = 2.4, lineend = "round") +
    geom_point(size = 2.6, colour = point_color) +
    geom_text(aes(y = upper_90, label = sprintf("f = %.2f", pd)), vjust = -0.9, size = 5.3) +
    scale_x_discrete(labels = c(Overstory = "\U0001F332", Understory = "\U0001F33F")) +
    scale_y_continuous(position = "right", limits = zoom_range(variable_name), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = title) +
    finite_axis_theme +
    theme(axis.text.x = element_text(size = 22, colour = "black"),
          plot.margin = margin(10, 5, 5, 4)) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0, margin = margin(b = 2)),
          plot.title.position = "panel")
}

build_mean_panel <- function(variable_name, point_color, y_lab, zoom_lab, titles = list(NULL, NULL)) {
  (build_violin_panel(variable_name, point_color, y_lab, titles[[1]]) |
     build_zoom_panel(variable_name, point_color, zoom_lab, titles[[2]])) +
    plot_layout(widths = c(1, 0.6))
}

build_lag_panel <- function(curve_df, variable_name, points_df, rate_col, point_color, y_lab, title = NULL, subtitle = NULL,
                            show_y_title = TRUE) {
  curve_sub <- curve_df %>% filter(variable == variable_name)
  
  ggplot() +
    geom_ribbon(data = curve_sub, aes(x = elevation, ymin = rate_lo, ymax = rate_hi),
                fill = point_color, alpha = 0.15) +
    geom_line(data = curve_sub, aes(x = elevation, y = rate),
              color = point_color, linewidth = 1) +
    geom_line(data = curve_sub, aes(x = elevation, y = rate_expected),
              color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_point(data = points_df, aes(x = elevation, y = .data[[rate_col]]),
               color = point_color, alpha = 0.4, size = 1.4) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    coord_cartesian(ylim = y_limits(variable_name), expand = FALSE) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    finite_axis_theme +
    theme(plot.title = element_text(size = 23, face = "bold", hjust = 0, margin = margin(b = 10)),
          plot.subtitle = element_text(size = 21, face = "bold", hjust = 0, margin = margin(b = 4)),
          plot.title.position = "plot")
}

# 3. Build the panels (same color within a row, across columns) -----------------
p_CTI_mean <- build_mean_panel("CTI", THERMO_COLOR, "Thermal affinity shift\nCTI (°C/decade)",
                               "Thermal affinity shift\nCTI (°C/decade)", titles = list("A   Across plots", NULL))
p_CMI_mean <- build_mean_panel("CMI", XERO_COLOR,   "Moisture affinity shift\nCMI (mm/decade)",
                               "Moisture affinity shift\nCMI (mm/decade)")
p_CTI_over   <- build_lag_panel(overstory_curves,  "CTI", overstory_rates,  "rate_CTI", THERMO_COLOR, "Thermal affinity shift\nCTI (°C/decade)",
                                title = "B   Along elevation", subtitle = "\U0001F332 Overstory")
p_CTI_under  <- build_lag_panel(understory_curves, "CTI", understory_rates, "rate_CTI", THERMO_COLOR, "Thermal affinity shift\nCTI (°C/decade)",
                                title = " ", subtitle = "\U0001F33F Understory", show_y_title = FALSE)
p_CMI_over   <- build_lag_panel(overstory_curves,  "CMI", overstory_rates,  "rate_CMI", XERO_COLOR,   "Moisture affinity shift\nCMI (mm/decade)")
p_CMI_under  <- build_lag_panel(understory_curves, "CMI", understory_rates, "rate_CMI", XERO_COLOR,   "Moisture affinity shift\nCMI (mm/decade)",
                                show_y_title = FALSE)


# 4. Shared legend: 3 keys (thermo observed, xero observed, expected), color +
#    linetype mapped to the same factor so they merge into one legend ------------
legend_levels <- c("Observed CTI rate", "Observed CMI rate", "Expected rate (space-for-time; dashed in B, grey in A)")

legend_dummy_df <- tibble(
  x    = rep(c(1, 2), 3),
  y    = rep(c(1, 2), 3),
  type = factor(rep(legend_levels, each = 2), levels = legend_levels)
)

p_legend_dummy <- ggplot(legend_dummy_df, aes(x = x, y = y, color = type, linetype = type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Observed CTI rate" = THERMO_COLOR,
    "Observed CMI rate"   = XERO_COLOR,
    "Expected rate (space-for-time; dashed in B, grey in A)"         = "black"
  ), name = NULL) +
  scale_linetype_manual(values = c(
    "Observed CTI rate" = "solid",
    "Observed CMI rate"   = "solid",
    "Expected rate (space-for-time; dashed in B, grey in A)"         = "dashed"
  ), name = NULL) +
  theme_ipsum(base_size = 18) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 19, margin = margin(r = 28)),
        legend.key.width = unit(1.6, "cm"))

shared_legend <- get_legend(p_legend_dummy)

# 5. Grid, shared x-axis label and shared legend, stacked together; the section
#    titles (A across plots, B along elevation) are the titles of the first
#    panel of each section -----------------------------------------------------

# wrap_plots() keeps the two-part summary panels nested as single cells
p_grid <- wrap_plots(list(p_CTI_mean, p_CTI_over, p_CTI_under, p_CMI_mean, p_CMI_over, p_CMI_under),
                     ncol = 3, widths = c(0.75, 1, 1))

x_axis_label <- (plot_spacer() | wrap_elements(grid::textGrob("Elevation (m)", gp = grid::gpar(fontsize = 21)))) +
  plot_layout(widths = c(0.75, 2))

# Shared legend under the figure (TRUE) or above it (FALSE)
legend_at_bottom <- TRUE
p_lag_integration <- if (legend_at_bottom) {
  p_grid / x_axis_label / wrap_elements(shared_legend) +
    plot_layout(heights = c(1, 0.04, 0.06))
} else {
  wrap_elements(shared_legend) / p_grid / x_axis_label +
    plot_layout(heights = c(0.06, 1, 0.04))
}
p_lag_integration

fig_file <- this_module_dir("figures", "thermo_xero_lag_overstory_understory.png")
ggsave(fig_file, p_lag_integration,
       width = 16.5, height = 11.5, dpi = 300, bg = "white")

# 6. Summary --------------------------------------------------------------------
cat("\n========================================\n")
cat("THERMO/XERO LAG INTEGRATION (OVERSTORY x UNDERSTORY) -- COMPLETE\n")
cat("========================================\n")
cat("Saved:\n")
cat(" - figures/thermo_xero_lag_overstory_understory.png\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
