# ------------------------------------------------------------------------------
# 01_compute_cti_cmi.R
# Module: pipeline/02_analysis/01_thermophilization/01_overstory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Community temperature index (CTI) and community moisture index (CMI) per
# overstory plot and survey, and their rates of change between surveys.
#
# CTI and CMI are abundance-weighted means of species' thermal (MAT) and
# moisture (MAP) affinities. Uncertainty in the affinities is propagated with
# the 100 resampled affinity values per species (02_integrate_gbif_chelsa):
# the index is computed once per resample, and the median and the 5th and
# 95th percentiles across resamples are kept. The percentiles describe the
# uncertainty of the community index due to uncertainty in the affinities;
# they are not tails of the community.
#
# Species: shared pool with climatic affinity
# (02_species_pool/outputs/overstory_species_pool_broad_full_clim_cov.csv),
# i.e. species recorded in both surveys. Vegetation, pool and affinity tables
# all use the vegetation names, so no name translation is needed.
#
# Weights: basal area as % of plot area (main) and stem count, computed in
# parallel. Stem density (count / plot area) gives the same weighted mean as
# count, since all stems of a plot share the same area, and is not computed.
#
# Survey interval: present-day year from 01_cleaning/outputs/
# present_census_year.csv (census year, or the mean of the two census years
# where both censuses are used); past year from plot_info.csv. Rates are
# (present - past) / interval in decades. A positive CMI rate indicates a
# shift towards species with wetter affinities.
#
# Climatic debt: climate change rate (Δ MAT, Δ MAP from
# 03_climate/02_explore_climate) minus the community rate; positive values
# indicate that the community lags behind climate change.
#
# Outputs:
#   data_clean/overstory_CTI_CMI_boot.csv       index per resample
#   data_clean/overstory_CTI_CMI_quantiles.csv  median, q05, q95 per plot x
#                                                 survey x weight
#   outputs/overstory_thermo_xero_rates.csv     rates and debt per plot x weight
#   tmp_files/overstory_species_match_check.csv species excluded, with reason
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(Hmisc)   # wtd.mean

this_module_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)
cleaning_dir     <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
affinity_dir     <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
climate_dir      <- function(...) here("pipeline", "01_data_prep", "03_climate", "02_explore_climate", ...)
species_pool_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)

for (d in c("data_clean", "outputs", "tmp_files")) {
  dir.create(this_module_dir(d), recursive = TRUE, showWarnings = FALSE)
}

stratum     <- "overstory"
weight_vars <- c("basal_area_spp_perc", "count")
N_RESAMPLES_EXPECTED <- 100

# 1. Inputs ---------------------------------------------------------------------
veg_long  <- read_csv(cleaning_dir("outputs", "overstory_long.csv"), show_col_types = FALSE)
plot_info <- read_csv(cleaning_dir("outputs", "plot_info.csv"), show_col_types = FALSE)
present_year <- read_csv(cleaning_dir("outputs", "present_census_year.csv"), show_col_types = FALSE) %>%
  filter(stratum == !!stratum) %>%
  select(plot_id, year_present = year, census_years)
affinity_resamples <- read_csv(affinity_dir("data_clean", "species_climatic_affinity_resamples.csv"),
                               show_col_types = FALSE)
species_pool <- read_csv(species_pool_dir("outputs", paste0(stratum, "_species_pool_broad_full_clim_cov.csv")),
                         show_col_types = FALSE) %>%
  pull(species)

cat("Records:", nrow(veg_long), "| species in shared pool with affinity:", length(species_pool), "\n")

# 2. Plot attributes and survey interval -----------------------------------------
# Elevation and coordinates are fixed plot attributes recorded with the past survey
plot_fixed <- plot_info %>%
  filter(period == "past") %>%
  distinct(plot_id, .keep_all = TRUE) %>%
  select(plot_id, elevation, lat, long, year_past = year)

year_gap <- plot_fixed %>%
  select(plot_id, year_past) %>%
  inner_join(present_year, by = "plot_id") %>%
  mutate(year_gap_decades = (year_present - year_past) / 10)

missing_year <- setdiff(unique(veg_long$plot_id), year_gap$plot_id)
if (length(missing_year) > 0) {
  cat("Plots without a past or present survey year (no rate):", paste(missing_year, collapse = ", "), "\n")
}

# 3. Species ---------------------------------------------------------------------
veg_species      <- unique(veg_long$species)
affinity_species <- unique(affinity_resamples$species)
species_matched  <- intersect(intersect(veg_species, affinity_species), species_pool)

tibble(species = setdiff(veg_species, species_matched)) %>%
  mutate(reason = case_when(!species %in% species_pool     ~ "not in the shared pool",
                            !species %in% affinity_species ~ "no climatic affinity",
                            TRUE                           ~ "other")) %>%
  write_csv(this_module_dir("tmp_files", paste0(stratum, "_species_match_check.csv")))

missing_pool <- setdiff(species_pool, species_matched)
if (length(missing_pool) > 0) {
  stop("Pool species without vegetation records or affinity: ", paste(missing_pool, collapse = ", "))
}
cat("Species used:", length(species_matched), "\n")

rep_counts <- unique(count(affinity_resamples, species)$n)
if (!all(rep_counts == N_RESAMPLES_EXPECTED)) {
  cat("Note: affinity resamples per species:", paste(rep_counts, collapse = ", "), "\n")
}

veg_matched <- veg_long %>% filter(species %in% species_matched, occurrence == 1)

# 4. CTI and CMI per resample ------------------------------------------------------
compute_cwm_boot <- function(weight_var) {
  veg_matched %>%
    select(plot_id, period, species, weight = all_of(weight_var)) %>%
    filter(!is.na(weight), weight > 0) %>%
    inner_join(affinity_resamples, by = "species", relationship = "many-to-many") %>%
    group_by(plot_id, period, rep) %>%
    summarise(CTI_boot  = Hmisc::wtd.mean(MAT_resample, weights = weight, na.rm = TRUE),
              CMI_boot  = Hmisc::wtd.mean(MAP_resample, weights = weight, na.rm = TRUE),
              n_species = n_distinct(species),
              .groups   = "drop") %>%
    mutate(weight_var = weight_var)
}
cti_cmi_boot <- map_dfr(weight_vars, compute_cwm_boot)
write_csv(cti_cmi_boot, this_module_dir("data_clean", paste0(stratum, "_CTI_CMI_boot.csv")))

# 5. Median and 5th / 95th percentiles across resamples ---------------------------
cti_cmi_quantiles <- cti_cmi_boot %>%
  group_by(plot_id, period, weight_var) %>%
  summarise(CTI_q05    = quantile(CTI_boot, 0.05, na.rm = TRUE),
            CTI_median = quantile(CTI_boot, 0.50, na.rm = TRUE),
            CTI_q95    = quantile(CTI_boot, 0.95, na.rm = TRUE),
            CMI_q05    = quantile(CMI_boot, 0.05, na.rm = TRUE),
            CMI_median = quantile(CMI_boot, 0.50, na.rm = TRUE),
            CMI_q95    = quantile(CMI_boot, 0.95, na.rm = TRUE),
            n_species  = first(n_species),
            .groups    = "drop") %>%
  left_join(plot_fixed %>% select(-year_past), by = "plot_id")
write_csv(cti_cmi_quantiles, this_module_dir("data_clean", paste0(stratum, "_CTI_CMI_quantiles.csv")))

# 6. Rates for plots with both surveys -------------------------------------------
thermo_xero_rates <- cti_cmi_quantiles %>%
  group_by(plot_id, weight_var) %>%
  filter(n_distinct(period) == 2) %>%
  ungroup() %>%
  pivot_wider(id_cols     = c(plot_id, weight_var, elevation, lat, long),
              names_from  = period,
              values_from = c(CTI_q05, CTI_median, CTI_q95, CMI_q05, CMI_median, CMI_q95, n_species),
              names_sep   = ".") %>%
  inner_join(year_gap, by = "plot_id") %>%
  mutate(delta_CTI = CTI_median.present - CTI_median.past,
         delta_CMI = CMI_median.present - CMI_median.past,
         rate_CTI  = delta_CTI / year_gap_decades,   # degC per decade
         rate_CMI  = delta_CMI / year_gap_decades)   # mm per decade

# 7. Climatic debt ------------------------------------------------------------------
climate_deltas <- read_csv(climate_dir("data_clean", "climate_past_present_wide.csv"),
                           show_col_types = FALSE) %>%
  select(plot_id, delta_MAT, delta_MAP)
thermo_xero_rates <- thermo_xero_rates %>%
  left_join(climate_deltas, by = "plot_id") %>%
  mutate(climate_rate_MAT = delta_MAT / year_gap_decades,
         climate_rate_MAP = delta_MAP / year_gap_decades,
         debt_CTI = climate_rate_MAT - rate_CTI,
         debt_CMI = climate_rate_MAP - rate_CMI)

write_csv(thermo_xero_rates, this_module_dir("outputs", paste0(stratum, "_thermo_xero_rates.csv")))

# 8. Summary ---------------------------------------------------------------------
cat("\nPlots with both surveys and a survey interval, per weight:\n")
print(count(thermo_xero_rates, weight_var))
cat("Survey interval (decades): mean", round(mean(year_gap$year_gap_decades), 2),
    "| range", paste(round(range(year_gap$year_gap_decades), 2), collapse = "-"), "\n")
print(thermo_xero_rates %>%
        group_by(weight_var) %>%
        summarise(n_plots = n(),
                  mean_rate_CTI = mean(rate_CTI, na.rm = TRUE),
                  mean_rate_CMI = mean(rate_CMI, na.rm = TRUE),
                  .groups = "drop"))
