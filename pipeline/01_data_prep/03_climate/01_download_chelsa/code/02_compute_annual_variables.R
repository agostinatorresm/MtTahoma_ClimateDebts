# ------------------------------------------------------------------------------
# 02_compute_annual_variables.R
# Module: pipeline/01_data_prep/03_climate/02_explore_climate
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: process the monthly CHELSA CSVs produced by
# 01_download_chelsa/code/01_chelsa_download_monthly.R (tas, tasmin, pr; full 1979-
# 2021 monthly record; one row per plot_id x month x year) into:
#   1. A full annual time series per plot (one row per plot_id x year), for
#      describing how climate has changed at the site continuously over the
#      whole record -- not just at the two survey windows.
#   2. A past (~1975 survey) vs. present (~2015 survey) window summary table
#      with delta columns, for the climate-change description figures and as
#      candidate HMSC covariates.
#
# For each plot and year, this computes:
#   - MAT   : annual mean temperature (mean of the 12 monthly tas values)
#   - MAP   : annual precipitation (sum of the 12 monthly pr values)
#   - Seasonal mean temperature and total precipitation (DJF/MAM/JJA/SON)
#   - PAS   : precipitation as snow, approximated as the sum of monthly
#             precipitation in months where mean monthly temperature <= 0 C,
#             within each calendar year. This is a threshold approximation;
#             ClimateNA/ClimateWNA (Wang et al. 2012) instead apply
#             month-specific logistic functions of monthly temperature
#   - MAT_min : annual mean of monthly minimum temperature (from tasmin)
#
# Season definition note: DJF/MAM/JJA/SON are assigned by calendar month, with
# December grouped into the DJF of the SAME calendar year (not the following
# winter, which would pull December into the next year's window). This is a
# simplification of the meteorological convention (which spans Dec of year Y-1
# with Jan-Feb of year Y).
#
# Coverage note: CHELSA's monthly archive has no tas/tasmin for 1979 (pr for
# 1979 is available), so 1979's MAT/MAT_min/PAS/seasonal temperature are NA
# for every plot while MAP/seasonal precipitation are not -- carried through
# as real NAs in the annual series rather than dropped, so the gap is visible
# in the trend plot instead of silently interpolated over.
#
# Past/present windows (must match 01_download_chelsa/code/01_chelsa_download_monthly.R's
# coverage): past = 1979-1983 (nearest available window to the ~1975 survey;
# 4 years of temperature data within it due to the 1979 gap), present =
# 2013-2017 (centered on the ~2015 survey). Years outside both windows are
# kept in the full annual series (for the trend plot) but are not part of
# the past/present summary table.
#
# Output (data_clean/):
#   - climate_annual_long.csv: one row per plot_id x year, full 1979-2021
#     series (all computed variables + a window column: "past"/"present"/NA)
#   - climate_past_present_wide.csv: one row per plot_id, past/present means
#     and delta columns, joined with plot elevation (from 01_download_chelsa's
#     coords_lookup.csv) for the climate-change description and PCA
# ------------------------------------------------------------------------------

# 0. Libraries -------------------------------------------------------------------
library(tidyverse)
library(here)

# 0a. Module directories ----------------------------------------------------------
download_dir <- function(...) here("pipeline", "01_data_prep", "03_climate", "01_download_chelsa", ...)
this_module_dir <- function(...) here("pipeline", "01_data_prep", "03_climate", "02_explore_climate", ...)

dir.create(this_module_dir("data_clean"), recursive = TRUE, showWarnings = FALSE)

# 1. Windows for the past/present summary (subset of the full annual series) -----
past_window    <- 1979:1983
present_window <- 2013:2017

window_lookup <- bind_rows(
  tibble(year = past_window,    window = "past"),
  tibble(year = present_window, window = "present")
)

# 2. Read all available monthly tas, tasmin, pr -----------------------------------
read_variable_long <- function(var, data_dir) {
  files <- list.files(data_dir, pattern = paste0("^", var, "_\\d{4}\\.csv$"), full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No files found for ", var, " in ", data_dir)
    return(tibble())
  }
  
  map_dfr(files, read_csv, show_col_types = FALSE)
}

assign_season <- function(month) {
  case_when(
    month %in% c(12, 1, 2) ~ "DJF",
    month %in% c(3, 4, 5)  ~ "MAM",
    month %in% c(6, 7, 8)  ~ "JJA",
    month %in% c(9, 10, 11) ~ "SON",
    TRUE ~ NA_character_
  )
}

# sum(x, na.rm = TRUE) of an all-NA vector silently returns 0, not NA -- unlike
# mean(), which returns NaN (treated as NA downstream) in that case. This
# matters for the plots with no valid coordinates (all 12 months NA for every
# year): their precipitation would otherwise come out as a false 0 mm instead
# of a missing value, contaminating the past/present summaries and plots.
sum_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

tas_long <- read_variable_long("tas", download_dir("data_raw", "chelsa_monthly", "tas")) %>%
  rename(tas = value)

tasmin_long <- read_variable_long("tasmin", download_dir("data_raw", "chelsa_monthly", "tasmin")) %>%
  rename(tasmin = value)

pr_long <- read_variable_long("pr", download_dir("data_raw", "chelsa_monthly", "pr")) %>%
  rename(pr = value)

cat("Years with tas data:", n_distinct(tas_long$year), "\n")
cat("Years with pr data:", n_distinct(pr_long$year), "\n")

# 3. Annual variables (per plot, per year) -----------------------------------------
annual_tas <- tas_long %>%
  group_by(plot_id, year) %>%
  summarise(MAT = mean(tas, na.rm = TRUE), .groups = "drop")

annual_tasmin <- tasmin_long %>%
  group_by(plot_id, year) %>%
  summarise(MAT_min = mean(tasmin, na.rm = TRUE), .groups = "drop")

annual_pr <- pr_long %>%
  group_by(plot_id, year) %>%
  summarise(MAP = sum_or_na(pr), .groups = "drop")

# 4. Seasonal variables (per plot, per year, per season) ---------------------------
seasonal_tas <- tas_long %>%
  mutate(season = assign_season(month)) %>%
  filter(!is.na(season)) %>%
  group_by(plot_id, year, season) %>%
  summarise(temp = mean(tas, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = season, values_from = temp, names_prefix = "temp_")

seasonal_pr <- pr_long %>%
  mutate(season = assign_season(month)) %>%
  filter(!is.na(season)) %>%
  group_by(plot_id, year, season) %>%
  summarise(precip = sum_or_na(pr), .groups = "drop") %>%
  pivot_wider(names_from = season, values_from = precip, names_prefix = "precip_")

# 5. Precipitation as snow (PAS): sum of monthly pr where monthly tas <= 0 --------
pas_annual <- tas_long %>%
  inner_join(pr_long, by = c("plot_id", "year", "month")) %>%
  filter(tas <= 0) %>%
  group_by(plot_id, year) %>%
  summarise(PAS = sum(pr, na.rm = TRUE), .groups = "drop")

# Plots/years with a tas value but no month <= 0 C get PAS = 0, not NA; a
# plot/year with no tas at all (the 1979 gap) keeps PAS = NA via the join
# below instead of being forced to 0.
# Note: the MAT-based filter has to happen BEFORE distinct(plot_id, year),
# since distinct() without .keep_all = TRUE drops every other column
# (including MAT) -- filtering first keeps MAT available to filter on.
pas_annual <- annual_tas %>%
  filter(!is.na(MAT)) %>%
  distinct(plot_id, year) %>%
  left_join(pas_annual, by = c("plot_id", "year")) %>%
  mutate(PAS = replace_na(PAS, 0))

# 6. Combine into one full annual series (all years, all plots) -------------------
climate_annual_long <- annual_pr %>%
  full_join(annual_tas, by = c("plot_id", "year")) %>%
  full_join(annual_tasmin, by = c("plot_id", "year")) %>%
  full_join(seasonal_tas, by = c("plot_id", "year")) %>%
  full_join(seasonal_pr, by = c("plot_id", "year")) %>%
  full_join(pas_annual, by = c("plot_id", "year")) %>%
  left_join(window_lookup, by = "year") %>%
  arrange(plot_id, year)

write_csv(climate_annual_long, this_module_dir("data_clean", "climate_annual_long.csv"))

cat("\nFull annual series:", nrow(climate_annual_long), "rows (plot x year)\n")
cat("Year range:", min(climate_annual_long$year), "-", max(climate_annual_long$year), "\n")

# Note: a plot can be missing MAT in every year (e.g. the plots with no
# lat/lon coordinates, which never got a CHELSA extraction), which shows up
# as "missing MAT" in every year even though the underlying issue is
# plot-specific, not year-specific -- so this reports both counts rather
# than just which years have any NA, to avoid the misleading impression of
# a systematic per-year data gap.
n_plots_total <- n_distinct(climate_annual_long$plot_id)
plots_missing_mat_always <- climate_annual_long %>%
  group_by(plot_id) %>%
  summarise(all_na = all(is.na(MAT)), .groups = "drop") %>%
  filter(all_na) %>%
  pull(plot_id)

cat(length(plots_missing_mat_always), "of", n_plots_total,
    "plots have no MAT value in any year (no valid coordinates):\n")
if (length(plots_missing_mat_always) > 0) {
  cat("  ", paste(plots_missing_mat_always, collapse = ", "), "\n")
}
cat("Years with at least one plot without MAT (1979 always, because of the CHELSA gap;",
    "otherwise the plots without coordinates listed above):\n")
cat("  ", paste(sort(unique(climate_annual_long$year[is.na(climate_annual_long$MAT)])), collapse = ", "), "\n")

# 7. Past/present window summary (mean across each window's years, per plot) -------
climate_by_window <- climate_annual_long %>%
  filter(!is.na(window)) %>%
  group_by(plot_id, window) %>%
  summarise(across(c(MAT, MAT_min, MAP, PAS,
                     starts_with("temp_"), starts_with("precip_")),
                   ~ mean(.x, na.rm = TRUE)),
            n_years_temp   = sum(!is.na(MAT)),
            n_years_precip = sum(!is.na(MAP)),
            n_years_PAS    = sum(!is.na(PAS)),
            .groups = "drop")

climate_past_present_wide <- climate_by_window %>%
  pivot_wider(
    id_cols = plot_id,
    names_from = window,
    values_from = c(MAT, MAT_min, MAP, PAS,
                    starts_with("temp_"), starts_with("precip_"),
                    n_years_temp, n_years_precip, n_years_PAS)
  )

delta_vars <- c("MAT", "MAT_min", "MAP", "PAS",
                "temp_DJF", "temp_MAM", "temp_JJA", "temp_SON",
                "precip_DJF", "precip_MAM", "precip_JJA", "precip_SON")

for (v in delta_vars) {
  past_col <- paste0(v, "_past")
  present_col <- paste0(v, "_present")
  if (all(c(past_col, present_col) %in% names(climate_past_present_wide))) {
    climate_past_present_wide[[paste0("delta_", v)]] <-
      climate_past_present_wide[[present_col]] - climate_past_present_wide[[past_col]]
  }
}

# 8. Join plot elevation (from coords_lookup.csv, if it has one) -------------------
coords_lookup <- read_csv(download_dir("outputs", "coords_lookup.csv"), show_col_types = FALSE)

if ("elevation" %in% names(coords_lookup)) {
  climate_past_present_wide <- climate_past_present_wide %>%
    left_join(coords_lookup %>% select(plot_id, elevation), by = "plot_id")
  cat("\nElevation added to climate_past_present_wide.csv.\n")
} else {
  cat("\nWARNING: coords_lookup.csv has no elevation column -- climate_past_present_wide.csv is saved without it.\n")
  cat("See find_col() in 01_download_chelsa/code/01_chelsa_download_monthly.R for the accepted column names.\n")
}

write_csv(climate_past_present_wide, this_module_dir("data_clean", "climate_past_present_wide.csv"))

cat("\n========================================\n")
cat("PROCESSING COMPLETE\n")
cat("========================================\n")
cat("Plots processed:", n_distinct(climate_annual_long$plot_id), "\n")
cat("Saved in pipeline/01_data_prep/03_climate/02_explore_climate/data_clean/\n")
cat("  - climate_annual_long.csv (full series 1979-2021, for the trend figure)\n")
cat("  - climate_past_present_wide.csv (past/present summary + deltas + elevation, for the zoom panel and the PCA)\n")
