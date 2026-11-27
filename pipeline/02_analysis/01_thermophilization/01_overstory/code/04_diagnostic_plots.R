# ------------------------------------------------------------------------------
# 04_diagnostic_plots.R
# Module: pipeline/02_analysis/01_thermophilization/01_overstory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: five complementary diagnostic figures for the overstory CTI/CMI
# results, beyond the elevation-lag figures in 02_thermo_xero_elevation_lag.R
# and the past/present-by-elevation figures in
# 03_CTI_CMI_past_present_elevation.R:
#
#   1. Trajectories: each plot's CTI/CMI change shown as a vertical arrow at
#      its own elevation (past -> present), to see plot-to-plot heterogeneity
#      that an aggregated LOESS curve averages away.
#   2. Climatic niche biplot: CTI vs CMI together (not each against elevation
#      separately), past and present, connected per plot -- shows how the
#      community moves in combined temperature x moisture space.
#   3. Climatic debt scatter: each plot's community response rate
#      (rate_CTI/rate_CMI) against its own measured climate-change rate
#      (climate_rate_MAT/climate_rate_MAP), with a 1:1 line -- points below
#      the line are responding more slowly than their local climate changed
#      (a visual version of debt_CTI/debt_CMI from 01_compute_cti_cmi.R).
#   4. Rate distributions: histograms of rate_CTI/rate_CMI across plots, to
#      see whether the (small) mean rates reported in 01_compute_cti_cmi.R
#      reflect a broadly shifted distribution or a mix of large opposing
#      shifts that average out.
#   5. Weighting-choice comparison: rate_CTI/rate_CMI boxplots across the
#      three community-weighting variables (basal_area_spp_perc, count,
#      density) computed in 01_compute_cti_cmi.R, to show whether the result
#      is robust to that choice.
#
# All figures use the primary weight variable (basal_area_spp_perc) except
# figure 5, which compares all three weight variables by design.
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)

# 0a. Module directories --------------------------------------------------------
this_module_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)

dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("tmp_files"), recursive = TRUE, showWarnings = FALSE)

PRIMARY_WEIGHT <- "basal_area_spp_perc"

PAST_COLOR    <- "#0A2463"  # navy/marine blue
PRESENT_COLOR <- "#FF6B5B"  # coral red
POINT_COLOR   <- "#D2691E"  # terracotta, matching 02_thermo_xero_elevation_lag.R

# 1. Load ------------------------------------------------------------------------
cti_cmi_quantiles <- read_csv(this_module_dir("data_clean", "overstory_CTI_CMI_quantiles.csv"), show_col_types = FALSE)
thermo_xero_rates <- read_csv(this_module_dir("outputs", "overstory_thermo_xero_rates.csv"), show_col_types = FALSE)

cq <- cti_cmi_quantiles %>% filter(weight_var == PRIMARY_WEIGHT)
rates_primary <- thermo_xero_rates %>% filter(weight_var == PRIMARY_WEIGHT)

# ================================================================================
# 1. TRAJECTORIES: per-plot CTI/CMI change as a vertical arrow at its elevation --
# ================================================================================
traj_wide <- cq %>%
  select(plot_id, period, elevation, CTI_median, CMI_median) %>%
  pivot_wider(names_from = period, values_from = c(CTI_median, CMI_median), names_sep = ".") %>%
  drop_na(CTI_median.past, CTI_median.present, CMI_median.past, CMI_median.present) %>%
  mutate(
    CTI_warmed = CTI_median.present > CTI_median.past,
    CMI_wetted = CMI_median.present > CMI_median.past
  )

# Diagnostic: some plots have valid CTI/CMI in both periods but no elevation
# -- this happens when a plot_id has vegetation records for period == "past"
# in overstory_long.csv but no matching period == "past" row in
# plot_info.csv, so 01_compute_cti_cmi.R's plot_fixed table (which sources
# elevation only from past-period plot_info.csv rows) has nothing to join.
# These plots are written out for review and dropped from anything plotted
# against elevation (they would otherwise show up as an unexplained "NA"
# facet/panel).
plots_missing_elevation <- traj_wide %>% filter(is.na(elevation)) %>% pull(plot_id)

if (length(plots_missing_elevation) > 0) {
  cat(length(plots_missing_elevation),
      "plot(s) have CTI/CMI in both periods but no elevation (no period == \"past\" row in plot_info.csv):\n")
  cat(paste(" -", plots_missing_elevation), sep = "\n")
  write_csv(tibble(plot_id = plots_missing_elevation),
            this_module_dir("tmp_files", "overstory_plots_missing_elevation.csv"))
  cat("Written to tmp_files/overstory_plots_missing_elevation.csv for review.\n")
  cat("Dropped from the trajectory/biplot figures below (elevation is required for both).\n")
  traj_wide <- traj_wide %>% filter(!is.na(elevation))
}

plot_trajectories <- function(df, y_past, y_pres, direction_col, y_lab, title_lab) {
  ggplot(df) +
    geom_segment(
      aes(x = elevation, xend = elevation, y = .data[[y_past]], yend = .data[[y_pres]],
          color = .data[[direction_col]]),
      arrow = arrow(length = unit(0.06, "inches")), alpha = 0.6, linewidth = 0.5
    ) +
    scale_color_manual(values = c("TRUE" = PRESENT_COLOR, "FALSE" = PAST_COLOR), guide = "none") +
    labs(x = "Elevation (m)", y = y_lab, title = title_lab) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "plain"))
}

p_traj_CTI <- plot_trajectories(
  traj_wide, "CTI_median.past", "CTI_median.present", "CTI_warmed",
  "CTI (°C)", "CTI: past → present per plot\n(coral = warmed, navy = cooled)"
)

p_traj_CMI <- plot_trajectories(
  traj_wide, "CMI_median.past", "CMI_median.present", "CMI_wetted",
  "CMI (mm)", "CMI: past → present per plot\n(coral = wetted, navy = dried)"
)

p_trajectories <- p_traj_CTI + p_traj_CMI + plot_layout(ncol = 2)
p_trajectories

ggsave(this_module_dir("figures", "overstory_CTI_CMI_trajectories_elevation.png"), p_trajectories,
       width = 12, height = 5.5, dpi = 300, bg = "white")

# ================================================================================
# 2. CLIMATIC NICHE BIPLOT: CTI vs CMI, past and present, connected per plot -----
# Faceted into 3 elevation bands (low/mid/high, by tercile -- roughly equal
# numbers of plots per panel) -- with all 324 plots on one panel the
# past-present arrows overlapped too much to read (a "spaghetti" plot);
# splitting by elevation keeps each panel legible and lets the direction of
# the community's shift in CTI x CMI space be compared across the gradient.
# ================================================================================
N_ELEV_BANDS <- 3

elev_band_breaks <- quantile(traj_wide$elevation, probs = seq(0, 1, length.out = N_ELEV_BANDS + 1), na.rm = TRUE)
elev_band_breaks_rounded <- round(elev_band_breaks)

elev_band_labels <- paste0(
  c("Low", "Mid", "High")[seq_len(N_ELEV_BANDS)], " elevation (",
  head(elev_band_breaks_rounded, -1), "–", elev_band_breaks_rounded[-1], " m)"
)

make_elev_band <- function(elevation) {
  cut(elevation, breaks = elev_band_breaks, include.lowest = TRUE, labels = elev_band_labels)
}

traj_wide <- traj_wide %>% mutate(elevation_band = make_elev_band(elevation))

# Restrict the plotted points to the same balanced (both-periods) plots as
# traj_wide -- cq on its own also includes single-period plots (a plot with
# data in only one period), which can sit at an elevation outside the
# breaks derived from traj_wide and would otherwise show up as a stray
# "NA" facet.
cq_banded <- cq %>%
  filter(plot_id %in% traj_wide$plot_id) %>%
  mutate(elevation_band = make_elev_band(elevation))

n_dropped_unbalanced <- n_distinct(cq$plot_id) - n_distinct(cq_banded$plot_id)
if (n_dropped_unbalanced > 0) {
  cat(n_dropped_unbalanced, "single-period (unbalanced) plots excluded from the biplot -- it only shows plots with both periods.\n")
}

# Sanity check: which plot sits in the highest elevation band, so its
# CTI/CMI values can be reviewed directly if that panel looks unusual.
highest_band <- levels(traj_wide$elevation_band)[max(as.integer(traj_wide$elevation_band), na.rm = TRUE)]
cat("\nPlot(s) in the highest elevation band (", highest_band, "):\n", sep = "")
print(traj_wide %>% filter(elevation_band == highest_band) %>%
        select(plot_id, elevation, CTI_median.past, CTI_median.present, CMI_median.past, CMI_median.present))

p_biplot <- ggplot() +
  geom_segment(
    data = traj_wide,
    aes(x = CTI_median.past, y = CMI_median.past, xend = CTI_median.present, yend = CMI_median.present),
    arrow = arrow(length = unit(0.05, "inches")), color = "grey50", alpha = 0.5, linewidth = 0.35
  ) +
  geom_point(
    data = cq_banded, aes(x = CTI_median, y = CMI_median, color = period), size = 1.8, alpha = 0.7
  ) +
  scale_color_manual(
    values = c(past = PAST_COLOR, present = PRESENT_COLOR),
    labels = c("Past (~1975)", "Present (~2015)"), name = NULL
  ) +
  facet_wrap(~elevation_band, scales = "free", nrow = 1) +
  labs(
    x = "Community temperature index, CTI (°C)",
    y = "Community moisture index, CMI (mm)",
    title = "Overstory community position in temperature x moisture space, past vs present",
    subtitle = paste0("Faceted into ", N_ELEV_BANDS, " elevation bands (low/mid/high, by tercile)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 14, face = "plain"),
    strip.text = element_text(size = 10)
  )
p_biplot

ggsave(this_module_dir("figures", "overstory_climate_biplot_CTI_CMI.png"), p_biplot,
       width = 13, height = 5.5, dpi = 300, bg = "white")

# ================================================================================
# 3. CLIMATIC DEBT SCATTER: community rate vs the plot's own climate-change rate -
# ================================================================================
plot_debt_scatter <- function(df, x_col, y_col, x_lab, y_lab, title_lab) {
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(color = POINT_COLOR, alpha = 0.5, size = 1.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    labs(x = x_lab, y = y_lab, title = title_lab) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 12, face = "plain"))
}

p_debt_CTI <- plot_debt_scatter(
  rates_primary, "climate_rate_MAT", "rate_CTI",
  "Climate warming rate (MAT, °C/decade)", "Community thermophilization rate (°C/decade)",
  "CTI vs measured climate-warming rate\n(1:1 dashed line = perfect tracking)"
)

p_debt_CMI <- plot_debt_scatter(
  rates_primary, "climate_rate_MAP", "rate_CMI",
  "Climate precipitation-change rate (MAP, mm/decade)", "Community response rate (CMI, mm/decade)",
  "CMI vs measured climate precipitation-change rate\n(1:1 dashed line = perfect tracking)"
)

p_debt <- p_debt_CTI + p_debt_CMI + plot_layout(ncol = 2)
p_debt

ggsave(this_module_dir("figures", "overstory_climatic_debt_scatter.png"), p_debt,
       width = 12, height = 5.5, dpi = 300, bg = "white")

# ================================================================================
# 4. RATE DISTRIBUTIONS: histograms of rate_CTI / rate_CMI across plots ---------
# ================================================================================
plot_rate_distribution <- function(df, rate_col, x_lab, title_lab) {
  ggplot(df, aes(x = .data[[rate_col]])) +
    geom_histogram(fill = POINT_COLOR, color = "white", bins = 30) +
    geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
    labs(x = x_lab, y = "Number of plots", title = title_lab) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 13, face = "plain"))
}

p_dist_CTI <- plot_rate_distribution(rates_primary, "rate_CTI", "Thermophilization rate (°C/decade)", "Distribution of CTI rate across plots")
p_dist_CMI <- plot_rate_distribution(rates_primary, "rate_CMI", "Moisture-affinity shift, CMI (mm/decade)", "Distribution of CMI rate across plots")

p_distributions <- p_dist_CTI + p_dist_CMI + plot_layout(ncol = 2)
p_distributions

ggsave(this_module_dir("figures", "overstory_rate_distributions.png"), p_distributions,
       width = 12, height = 5, dpi = 300, bg = "white")

# ================================================================================
# 5. WEIGHTING-CHOICE COMPARISON: rate_CTI / rate_CMI across the 3 weight_vars ---
# ================================================================================
plot_weight_comparison <- function(df, rate_col, y_lab, title_lab) {
  ggplot(df, aes(x = weight_var, y = .data[[rate_col]], fill = weight_var)) +
    geom_boxplot(alpha = 0.75, outlier.alpha = 0.3) +
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
    scale_fill_viridis_d(guide = "none") +
    labs(x = NULL, y = y_lab, title = title_lab) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(size = 13, face = "plain"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
}

p_weight_CTI <- plot_weight_comparison(thermo_xero_rates, "rate_CTI", "Thermophilization rate (°C/decade)", "CTI rate by community-weighting choice")
p_weight_CMI <- plot_weight_comparison(thermo_xero_rates, "rate_CMI", "Moisture-affinity shift, CMI (mm/decade)", "CMI rate by community-weighting choice")

p_weight_comparison <- p_weight_CTI + p_weight_CMI + plot_layout(ncol = 2)
p_weight_comparison

ggsave(this_module_dir("figures", "overstory_weight_var_comparison.png"), p_weight_comparison,
       width = 12, height = 5.5, dpi = 300, bg = "white")

# 6. Summary --------------------------------------------------------------------
cat("Saved:\n")
cat(" - figures/overstory_CTI_CMI_trajectories_elevation.png\n")
cat(" - figures/overstory_climate_biplot_CTI_CMI.png\n")
cat(" - figures/overstory_climatic_debt_scatter.png\n")
cat(" - figures/overstory_rate_distributions.png\n")
cat(" - figures/overstory_weight_var_comparison.png\n")
