# ------------------------------------------------------------------------------
# 02_thermo_xero_elevation_lag.R
# Module: pipeline/02_analysis/01_thermophilization/02_understory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: compare the OBSERVED thermophilization/xerophilization rate (from
# 01_compute_cti_cmi_understory.R) against an EXPECTED rate derived from a space-for-time
# substitution (Fadrique et al. 2018-style climatic-debt approach), and plot
# both against elevation. Understory counterpart of
# 01_overstory/code/02_thermo_xero_elevation_lag.R -- identical method, only
# the module directory and the primary weight variable change (cover instead
# of basal_area_spp_perc; see 01_compute_cti_cmi_understory.R for why the understory
# module uses cover/occurrence in place of the overstory's basal_area/count/
# density).
#
# EXPECTED RATE (space-for-time substitution): a linear relationship between
# each plot's baseline (past) community index (CTI or CMI) and its baseline
# climate (MAT_past for CTI, MAP_past for CMI) is fit across plots (a LOESS
# relationship is kept as a check; it cannot predict outside the past
# climate range, where many plots fall for MAP) -- this captures "what community composition looks like at a given
# climate", using spatial variation among plots as a proxy for temporal
# change. That same fitted relationship is then used to predict, for each
# plot, the CTI/CMI value "expected" at its actual PRESENT-day climate
# (MAT_present / MAP_present), still using the plot's OWN elevation-driven
# position in the relationship. The expected rate is the resulting expected
# delta, divided by that plot's actual year_gap_decades.
#
# OBSERVED RATE: taken directly from understory_thermo_xero_rates.csv
# (01_compute_cti_cmi_understory.R): the actual measured change in community CTI/CMI
# between the past and present surveys, per plot, per actual year gap.
#
# CLIMATIC LAG: expected_rate - observed_rate, per plot. A positive lag means
# the community changed less than expected from a full climate-tracking
# community (a "debt" -- composition lagging climate); a negative lag means
# the community changed more than the space-for-time relationship predicts.
#
# WEIGHT VARIABLE: figures use cover (PRIMARY_WEIGHT below), the abundance-
# weighted CWM, analogous to the overstory module's dominance-weighted
# choice. occurrence was also computed in 01_compute_cti_cmi_understory.R and is
# available in the same input files for a robustness check, but is not
# plotted here by default.
#
# BOOTSTRAP COMPONENT: the observed elevation-smoothed curve (with its 95% CI
# ribbon) uses only the bootstrap MEDIAN (CTI_median/CMI_median). The 95% CI
# shown here comes from LOESS standard errors (delta-method propagation
# across the past/present fits), not from the q05/q95 bootstrap components.
#
# Output:
#   - data_clean/understory_thermo_xero_expected.csv: per plot x weight_var --
#     expected_CTI_past/present, rate_CTI_expected, lag_CTI (and the CMI
#     equivalents)
#   - outputs/understory_elevation_curves.csv: elevation-grid LOESS curves
#     (observed rate + 95% CI, expected rate) for CTI and CMI, primary weight
#     variable only -- the data behind the figures below
#   - figures/understory_thermo_CTI_elevation_lag.png
#   - figures/understory_thermo_CMI_elevation_lag.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)

# 0a. Module directories --------------------------------------------------------
this_module_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "02_understory", ...)
climate_dir     <- function(...) here("pipeline", "01_data_prep", "03_climate", "02_explore_climate", ...)

for (d in c("data_clean", "outputs", "figures")) {
  dir.create(this_module_dir(d), recursive = TRUE, showWarnings = FALSE)
}

PRIMARY_WEIGHT <- "cover"
LOESS_SPAN <- 0.75

# 1. Load inputs -----------------------------------------------------------------
cti_cmi_quantiles <- read_csv(this_module_dir("data_clean", "understory_CTI_CMI_quantiles.csv"), show_col_types = FALSE)
thermo_xero_rates <- read_csv(this_module_dir("outputs", "understory_thermo_xero_rates.csv"), show_col_types = FALSE)

climate_wide <- read_csv(climate_dir("data_clean", "climate_past_present_wide.csv"), show_col_types = FALSE) %>%
  select(plot_id, MAT_past, MAT_present, MAP_past, MAP_present)

cat("Using weight_var =", PRIMARY_WEIGHT, "for all figures/curves in this script.\n")

# 2. Restrict to the primary weight variable, join climate -----------------------
cq <- cti_cmi_quantiles %>%
  filter(weight_var == PRIMARY_WEIGHT) %>%
  left_join(climate_wide, by = "plot_id")

past_df <- cq %>% filter(period == "past")
pres_df <- cq %>% filter(period == "present")

rates <- thermo_xero_rates %>%
  filter(weight_var == PRIMARY_WEIGHT) %>%
  left_join(climate_wide, by = "plot_id")

# Curves are fitted on the plots with an index in both surveys (the plots
# with a rate), so past and present curves describe the same set of plots
past_df <- past_df %>% filter(plot_id %in% rates$plot_id)
pres_df <- pres_df %>% filter(plot_id %in% rates$plot_id)
cat("Plots used for the elevation curves:", n_distinct(past_df$plot_id), "\n")

mean_decades <- mean(rates$year_gap_decades, na.rm = TRUE)
cat("Mean year_gap_decades across balanced plots:", round(mean_decades, 3),
    "-- used only for the elevation-smoothed OBSERVED curve below.",
    "Per-plot rates elsewhere (rates$rate_CTI, rates$rate_CMI, and the expected/lag",
    "values computed in section 4) all use each plot's own actual year_gap_decades.\n")

# 3. Observed rate vs elevation: LOESS(CTI_median ~ elevation) and
#    LOESS(CMI_median ~ elevation), fit separately for past and present,
#    predicted on a common elevation grid, with 95% CI via SE propagation ------
elev_range <- range(c(past_df$elevation, pres_df$elevation), na.rm = TRUE)
elev_grid <- seq(elev_range[1], elev_range[2], length.out = 100)

compute_observed_curve <- function(response_var) {
  loess_past <- loess(as.formula(paste(response_var, "~ elevation")), data = past_df, span = LOESS_SPAN)
  loess_pres <- loess(as.formula(paste(response_var, "~ elevation")), data = pres_df, span = LOESS_SPAN)
  
  pred_past <- predict(loess_past, newdata = data.frame(elevation = elev_grid), se = TRUE)
  pred_pres <- predict(loess_pres, newdata = data.frame(elevation = elev_grid), se = TRUE)
  
  delta   <- pred_pres$fit - pred_past$fit
  rate    <- delta / mean_decades
  se_rate <- sqrt(pred_past$se.fit^2 + pred_pres$se.fit^2) / mean_decades
  
  tibble(
    elevation = elev_grid,
    rate_observed = rate,
    rate_observed_lo = rate - 1.96 * se_rate,
    rate_observed_hi = rate + 1.96 * se_rate
  )
}

observed_CTI <- compute_observed_curve("CTI_median")
observed_CMI <- compute_observed_curve("CMI_median")

# 4. Expected rate: space-for-time substitution (CTI~MAT_past, CMI~MAP_past) -----
# Spatial relationship fit ONLY on past-period data (the elevation/climate
# gradient as it existed in ~1975), then used to predict what CTI/CMI "should"
# be at each plot's actual present-day climate.
# Main: linear relationship (index ~ past climate), which also predicts for
# plots whose present climate lies outside the past range. Check: LOESS
# relationship (span LOESS_SPAN), which does not extrapolate (NA outside the
# past range).
lm_CTI_MAT    <- lm(CTI_median ~ MAT_past, data = past_df)
lm_CMI_MAP    <- lm(CMI_median ~ MAP_past, data = past_df)
loess_CTI_MAT <- loess(CTI_median ~ MAT_past, data = past_df, span = LOESS_SPAN)
loess_CMI_MAP <- loess(CMI_median ~ MAP_past, data = past_df, span = LOESS_SPAN)
cat("Space-for-time slopes (past plots): CTI per degC MAT =", round(coef(lm_CTI_MAT)[2], 3),
    "| CMI per mm MAP =", round(coef(lm_CMI_MAP)[2], 3), "\n")

expected_delta <- function(fit, x_past, x_present, var) {
  nd_past <- setNames(data.frame(x_past), var)
  nd_pres <- setNames(data.frame(x_present), var)
  as.numeric(predict(fit, newdata = nd_pres)) - as.numeric(predict(fit, newdata = nd_past))
}

expected_by_plot <- rates %>%
  select(plot_id, elevation, year_gap_decades, MAT_past, MAT_present, MAP_past, MAP_present,
         rate_CTI, rate_CMI) %>%
  mutate(
    delta_CTI_expected       = expected_delta(lm_CTI_MAT, MAT_past, MAT_present, "MAT_past"),
    rate_CTI_expected        = delta_CTI_expected / year_gap_decades,
    lag_CTI                  = rate_CTI_expected - rate_CTI,
    rate_CTI_expected_loess  = expected_delta(loess_CTI_MAT, MAT_past, MAT_present, "MAT_past") / year_gap_decades,
    lag_CTI_loess            = rate_CTI_expected_loess - rate_CTI,
    
    delta_CMI_expected       = expected_delta(lm_CMI_MAP, MAP_past, MAP_present, "MAP_past"),
    rate_CMI_expected        = delta_CMI_expected / year_gap_decades,
    lag_CMI                  = rate_CMI_expected - rate_CMI,
    rate_CMI_expected_loess  = expected_delta(loess_CMI_MAP, MAP_past, MAP_present, "MAP_past") / year_gap_decades,
    lag_CMI_loess            = rate_CMI_expected_loess - rate_CMI
  )

cat("\nExpected-vs-observed rates computed for", nrow(expected_by_plot), "plots.\n")
# LOESS check: plots whose present climate lies outside the past range have
# no LOESS expectation
cat("Plots without a LOESS expectation (present climate outside the past range): CTI",
    sum(is.na(expected_by_plot$rate_CTI_expected_loess)), "| CMI",
    sum(is.na(expected_by_plot$rate_CMI_expected_loess)), "\n")
cat("Mean lag, LOESS check on plots within the past range: CTI",
    round(mean(expected_by_plot$lag_CTI_loess, na.rm = TRUE), 4), "| CMI",
    round(mean(expected_by_plot$lag_CMI_loess, na.rm = TRUE), 4), "\n")
cat("Mean lag_CTI (expected - observed, degC/decade):", round(mean(expected_by_plot$lag_CTI, na.rm = TRUE), 4), "\n")
cat("Mean lag_CMI (expected - observed, mm/decade):  ", round(mean(expected_by_plot$lag_CMI, na.rm = TRUE), 4), "\n")

write_csv(expected_by_plot, this_module_dir("data_clean", "understory_thermo_xero_expected.csv"))

# 5. Expected rate smoothed vs elevation, for the same grid used in section 3 ----
expected_CTI_smooth <- loess(rate_CTI_expected ~ elevation, data = expected_by_plot, span = LOESS_SPAN)
expected_CMI_smooth <- loess(rate_CMI_expected ~ elevation, data = expected_by_plot, span = LOESS_SPAN)

observed_CTI <- observed_CTI %>%
  rename(rate = rate_observed, rate_lo = rate_observed_lo, rate_hi = rate_observed_hi) %>%
  mutate(rate_expected = predict(expected_CTI_smooth, newdata = data.frame(elevation = elev_grid)))
observed_CMI <- observed_CMI %>%
  rename(rate = rate_observed, rate_lo = rate_observed_lo, rate_hi = rate_observed_hi) %>%
  mutate(rate_expected = predict(expected_CMI_smooth, newdata = data.frame(elevation = elev_grid)))

elevation_curves <- bind_rows(
  observed_CTI %>% mutate(variable = "CTI"),
  observed_CMI %>% mutate(variable = "CMI")
)

write_csv(elevation_curves, this_module_dir("outputs", "understory_elevation_curves.csv"))

# 6. Figures: observed LOESS curve + CI ribbon + expected line + per-plot points -
POINT_COLOR <- "#2E8B57"  # sea green, distinct from overstory's terracotta, keeps
# overstory/understory figures visually distinguishable
# when placed side by side in the manuscript

make_lag_figure <- function(curve_df, points_df, rate_col, y_lab, title_lab) {
  ggplot() +
    geom_ribbon(data = curve_df, aes(x = elevation, ymin = rate_lo, ymax = rate_hi),
                fill = POINT_COLOR, alpha = 0.15) +
    geom_line(data = curve_df, aes(x = elevation, y = rate),
              color = POINT_COLOR, linewidth = 1) +
    geom_line(data = curve_df, aes(x = elevation, y = rate_expected),
              color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_point(data = points_df, aes(x = elevation, y = .data[[rate_col]]),
               color = POINT_COLOR, alpha = 0.4, size = 1.6) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    labs(x = "Elevation (m)", y = y_lab, title = title_lab) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(size = 15, face = "plain"))
}

p_CTI <- make_lag_figure(
  observed_CTI, rates, "rate_CTI",
  "Thermophilization rate (°C / decade)",
  "Understory CTI: observed (green) vs expected (dashed black)"
)
p_CTI

p_CMI <- make_lag_figure(
  observed_CMI, rates, "rate_CMI",
  "Moisture-affinity shift, CMI (mm / decade)",
  "Understory CMI: observed (green) vs expected (dashed black)"
)
p_CMI

ggsave(this_module_dir("figures", "understory_thermo_CTI_elevation_lag.png"), p_CTI,
       width = 8, height = 6, dpi = 300, bg = "white")

ggsave(this_module_dir("figures", "understory_thermo_CMI_elevation_lag.png"), p_CMI,
       width = 8, height = 6, dpi = 300, bg = "white")

# 7. Summary --------------------------------------------------------------------
cat("\n========================================\n")
cat("OBSERVED VS EXPECTED THERMO/XERO RATES vs ELEVATION -- COMPLETE\n")
cat("========================================\n")
cat("Saved:\n")
cat(" - data_clean/understory_thermo_xero_expected.csv (per-plot expected rates + lag)\n")
cat(" - outputs/understory_elevation_curves.csv (elevation-grid curves behind the figures)\n")
cat(" - figures/understory_thermo_CTI_elevation_lag.png\n")
cat(" - figures/understory_thermo_CMI_elevation_lag.png\n")
