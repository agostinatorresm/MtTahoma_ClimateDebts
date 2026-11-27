# ------------------------------------------------------------------------------
# 05_check_results_numbers.R
# Module: pipeline (diagnostic)
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: print the numbers needed to write the Results section on
# climate-affinity shifts (thermophilization / xerophilization), from the
# outputs of 01_overstory and 02_understory.
# Diagnostic only; nothing is written to disk.
#
# Sections:
#   1. Index values and observed rates per layer (mean, 95% CI, share of
#      plots with a positive rate, test against zero)
#   2. Rates by elevation band
#   3. Elevation curves: where the 95% interval excludes zero; observed vs
#      expected
#   4. Climate change at the plots and climatic lag (expected - observed)
#   5. Community affinity tails (cold/dry vs warm/wet tail rates)
#   6. Species and group contributions (community contribution)
# The shared vs frequent pool comparison is reported by the sensitivity
# module, not here.
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

thermo_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", ...)

section <- function(title) cat("\n==============================\n", title, "\n==============================\n", sep = "")

read_if_exists <- function(path) {
  if (!file.exists(path)) { cat("  [missing]", path, "\n"); return(NULL) }
  read_csv(path, show_col_types = FALSE)
}

STRATA <- tibble(
  stratum = c("overstory", "understory"),
  dir     = c("01_overstory", "02_understory"),
  weight  = c("basal_area_spp_perc", "cover")
)

mean_ci <- function(x) {
  x <- x[!is.na(x)]
  tt <- t.test(x)
  tibble(n = length(x), mean = mean(x), lo95 = tt$conf.int[1], hi95 = tt$conf.int[2],
         median = median(x), pct_positive = 100 * mean(x > 0),
         p_ttest = tt$p.value, p_wilcoxon = suppressWarnings(wilcox.test(x)$p.value))
}

rates <- map(set_names(STRATA$stratum), function(st) {
  d <- STRATA %>% filter(stratum == st)
  read_if_exists(thermo_dir(d$dir, "outputs", paste0(st, "_thermo_xero_rates.csv"))) %>%
    { if (is.null(.)) NULL else filter(., weight_var == d$weight) }
})

# 1. Index values and observed rates --------------------------------------------------
section("1. INDEX VALUES AND OBSERVED RATES (primary weight)")

for (st in names(rates)) {
  d <- rates[[st]]; if (is.null(d)) next
  cat("\n---", st, "---\n")
  cat("Mean CTI_median: past", round(mean(d$CTI_median.past, na.rm = TRUE), 2),
      "| present", round(mean(d$CTI_median.present, na.rm = TRUE), 2), "degC\n")
  cat("Mean CMI_median: past", round(mean(d$CMI_median.past, na.rm = TRUE), 1),
      "| present", round(mean(d$CMI_median.present, na.rm = TRUE), 1), "mm\n")
  cat("\nThermophilization rate (degC/decade):\n"); print(mean_ci(d$rate_CTI), width = Inf)
  cat("\nMoisture-affinity shift rate (mm/decade):\n"); print(mean_ci(d$rate_CMI), width = Inf)
  cat("\nCorrelation between rate_CTI and rate_CMI across plots (Spearman):",
      round(cor(d$rate_CTI, d$rate_CMI, method = "spearman", use = "complete.obs"), 3), "\n")
}

# 2. Rates by elevation band ------------------------------------------------------------
section("2. RATES BY ELEVATION BAND (250 m)")

for (st in names(rates)) {
  d <- rates[[st]]; if (is.null(d)) next
  cat("\n---", st, "---\n")
  d %>%
    mutate(band = cut(elevation, breaks = seq(0, 2500, by = 250), include.lowest = TRUE)) %>%
    group_by(band) %>%
    summarise(n = n(), mean_rate_CTI = mean(rate_CTI, na.rm = TRUE),
              mean_rate_CMI = mean(rate_CMI, na.rm = TRUE), .groups = "drop") %>%
    filter(n > 0) %>% print(n = 20)
  cat("Spearman rate_CTI ~ elevation:",
      round(cor(d$rate_CTI, d$elevation, method = "spearman", use = "complete.obs"), 3),
      "| rate_CMI ~ elevation:",
      round(cor(d$rate_CMI, d$elevation, method = "spearman", use = "complete.obs"), 3), "\n")
}

# 3. Elevation curves ---------------------------------------------------------------------
section("3. ELEVATION CURVES (95% interval vs zero; observed vs expected)")

for (i in seq_len(nrow(STRATA))) {
  st <- STRATA$stratum[i]
  cur <- read_if_exists(thermo_dir(STRATA$dir[i], "outputs", paste0(st, "_elevation_curves.csv")))
  if (is.null(cur)) next
  cat("\n---", st, "---\n")
  cur %>%
    group_by(variable) %>%
    summarise(
      elev_min = min(elevation), elev_max = max(elevation),
      pct_grid_ci_above_zero = 100 * mean(rate_lo > 0),
      pct_grid_ci_below_zero = 100 * mean(rate_hi < 0),
      mean_observed = mean(rate, na.rm = TRUE),
      mean_expected = mean(rate_expected, na.rm = TRUE),
      n_grid_expected_NA = sum(is.na(rate_expected)),
      pct_grid_expected_above_observed = 100 * mean(rate_expected > rate, na.rm = TRUE),
      .groups = "drop") %>% print(width = Inf)
  # Elevation ranges where the observed 95% interval excludes zero
  for (v in unique(cur$variable)) {
    sub <- cur %>% filter(variable == v) %>% arrange(elevation)
    sig <- sub %>% filter(rate_lo > 0 | rate_hi < 0)
    cat(v, ": interval excludes zero between",
        if (nrow(sig) > 0) paste(round(min(sig$elevation)), "and", round(max(sig$elevation)), "m") else "-- never",
        "\n")
  }
}

# 4. Climate change and climatic lag --------------------------------------------------------
section("4. CLIMATE CHANGE AT THE PLOTS AND CLIMATIC LAG")

for (i in seq_len(nrow(STRATA))) {
  st <- STRATA$stratum[i]
  d <- rates[[st]]
  exp_df <- read_if_exists(thermo_dir(STRATA$dir[i], "data_clean", paste0(st, "_thermo_xero_expected.csv")))
  cat("\n---", st, "---\n")
  if (!is.null(d) && "climate_rate_MAT" %in% names(d)) {
    cat("Measured warming rate (degC/decade): mean", round(mean(d$climate_rate_MAT, na.rm = TRUE), 3),
        "| range", paste(round(range(d$climate_rate_MAT, na.rm = TRUE), 3), collapse = " to "), "\n")
    cat("Measured precipitation trend (mm/decade): mean", round(mean(d$climate_rate_MAP, na.rm = TRUE), 1),
        "| range", paste(round(range(d$climate_rate_MAP, na.rm = TRUE), 1), collapse = " to "), "\n")
    cat("Spearman observed rate_CTI ~ measured warming rate:",
        round(cor(d$rate_CTI, d$climate_rate_MAT, method = "spearman", use = "complete.obs"), 3), "\n")
    cat("Spearman observed rate_CMI ~ measured precipitation trend:",
        round(cor(d$rate_CMI, d$climate_rate_MAP, method = "spearman", use = "complete.obs"), 3), "\n")
  }
  if (!is.null(exp_df)) {
    if ("weight_var" %in% names(exp_df)) exp_df <- filter(exp_df, weight_var == STRATA$weight[i])
    cat("\nExpected thermophilization rate (degC/decade):\n"); print(mean_ci(exp_df$rate_CTI_expected), width = Inf)
    cat("\nLag CTI (expected - observed, degC/decade):\n");    print(mean_ci(exp_df$lag_CTI), width = Inf)
    cat("\nExpected moisture-affinity shift (mm/decade):\n");  print(mean_ci(exp_df$rate_CMI_expected), width = Inf)
    cat("\nLag CMI (expected - observed, mm/decade):\n");      print(mean_ci(exp_df$lag_CMI), width = Inf)
  }
}

# 5. Community affinity tails ----------------------------------------------------------------
section("5. COMMUNITY AFFINITY EDGES (weighted q10 / q50 / q90 of species affinities)")

for (i in seq_len(nrow(STRATA))) {
  st <- STRATA$stratum[i]
  e <- read_if_exists(thermo_dir(STRATA$dir[i], "outputs", paste0(st, "_affinity_edge_summary.csv")))
  if (!is.null(e)) { cat("\n---", st, "---\n"); print(e, n = 20, width = Inf) }
}

# 6. Species and group contributions ------------------------------------------------------------
# Ranked by community contribution (sum over plots / all plots with the index),
# i.e. the share of the mean index change attributable to each species.
section("6. SPECIES AND GROUP CONTRIBUTIONS (community contribution)")

show_top <- function(df, label, n = 5) {
  cols <- c("species", "n_plots", "community_contribution", "mean_contribution", "MAT_mean", "MAP_mean")
  cat("\n", label, "\n", sep = "")
  cat("  largest positive:\n"); print(df %>% slice_max(community_contribution, n = n) %>% select(any_of(cols)), width = Inf)
  cat("  largest negative:\n"); print(df %>% slice_min(community_contribution, n = n) %>% select(any_of(cols)), width = Inf)
  cat("  sum of community contributions (= mean index change):", signif(sum(df$community_contribution, na.rm = TRUE), 3), "\n")
}

for (i in seq_len(nrow(STRATA))) {
  st <- STRATA$stratum[i]
  for (v in c("CTI", "CMI")) {
    df <- read_if_exists(thermo_dir(STRATA$dir[i], "outputs", paste0(st, "_species_contribution_", v, ".csv")))
    if (!is.null(df)) show_top(df, paste(st, v))
  }
}

for (f in c("understory_lifeform_contribution.csv", "understory_woodiness_contribution.csv")) {
  df <- read_if_exists(thermo_dir("02_understory", "outputs", f))
  if (!is.null(df)) { cat("\n", f, "\n", sep = ""); print(df, n = 30, width = Inf) }
}

cat("\n========================================\n")
cat("RESULTS NUMBERS CHECK -- COMPLETE\n")
cat("========================================\n")
