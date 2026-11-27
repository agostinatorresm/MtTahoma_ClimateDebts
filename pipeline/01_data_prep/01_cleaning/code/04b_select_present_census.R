# =============================================================================
# 04b_select_present_census.R
# Author: Agostina Torres
#
# Defines which present-day census records are used per plot and layer.
# Most plots were resurveyed once (2014-2019); some plots were surveyed in
# two years.
#
#   - Plots with a documented choice of present-day year
#     (present_year_documented.csv, written by 01_clean_plot_ids.R from the
#     data provider's per-plot notes): only the census of that year is used.
#   - Plots surveyed in two years without a documented choice: both censuses
#     are used and averaged by the build scripts (05_build_overstory.R,
#     06_build_understory.R). These plots are listed in the output so that
#     they can be reported and examined in sensitivity analyses.
#   - Plots surveyed in two years without a documented choice, where one
#     understory census is incomplete (fewer than four cover subplots) and
#     the other complete: only the complete census is used, since averaging
#     would combine cover from different sampling efforts.
#   - Plots surveyed once: that census.
#
# The present-day year of a plot is the census year, or the mean of the two
# census years for averaged plots (used for the survey interval).
#
# Inputs : pipeline/01_data_prep/01_cleaning/outputs/tree_veg_plotid_clean.rda
#          pipeline/01_data_prep/01_cleaning/outputs/present_year_documented.csv
# Outputs: pipeline/01_data_prep/01_cleaning/outputs/
#            present_census_selected.csv  layer x plot x census year records
#                                         used by the build scripts
#            present_census_year.csv      one row per layer x plot: present-
#                                         day year and how it was defined
#            check_present_census_completeness.csv
#                                         every plot x year census with its
#                                         completeness measures
# =============================================================================

library(here)
library(tidyverse)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
subplot_codes <- c("QN", "QE", "QS", "QW")

clean_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
load(file.path(clean_dir, "tree_veg_plotid_clean.rda"))   # tree_R_clean, veg_R_clean, ...

documented <- read_csv(file.path(clean_dir, "present_year_documented.csv"),
                       show_col_types = FALSE) %>%
  select(plot_id, documented_year = year)

# -----------------------------------------------------------------------------
# 1. Censuses per plot, with completeness measures
# -----------------------------------------------------------------------------
census_overstory <- tree_R_clean %>%
  filter(!is.na(plot_id)) %>%
  group_by(plot_id, Year) %>%
  summarise(n_stems = n(), n_species = n_distinct(Species), .groups = "drop") %>%
  mutate(stratum = "overstory")

census_understory <- veg_R_clean %>%
  filter(!is.na(plot_id)) %>%
  mutate(Replicate = str_squish(Replicate)) %>%
  group_by(plot_id, Year) %>%
  summarise(n_subplots = n_distinct(Replicate[Replicate %in% subplot_codes]),
            n_records  = sum(Replicate %in% subplot_codes),
            n_species  = n_distinct(Species), .groups = "drop") %>%
  mutate(stratum = "understory")

# -----------------------------------------------------------------------------
# 2. Selection
# -----------------------------------------------------------------------------
census_all <- bind_rows(census_overstory, census_understory) %>%
  left_join(documented, by = "plot_id") %>%
  group_by(stratum, plot_id) %>%
  mutate(n_years       = n_distinct(Year),
         has_doc_year  = any(Year == documented_year, na.rm = TRUE),
         # Understory census incomplete (fewer than four cover subplots)
         # while the other census of the plot is complete
         incomplete    = stratum == "understory" & coalesce(n_subplots, 0L) < length(subplot_codes),
         has_complete  = any(!incomplete),
         selection = case_when(n_years == 1                 ~ "single census",
                               has_doc_year                 ~ "documented year",
                               any(incomplete) & has_complete ~ "complete census",
                               TRUE                         ~ "two censuses averaged"),
         used = case_when(selection == "documented year" ~ Year == documented_year,
                          selection == "complete census" ~ !incomplete,
                          TRUE                           ~ TRUE)) %>%
  ungroup() %>%
  select(stratum, plot_id, Year, documented_year, n_years, selection, used,
         n_stems, n_subplots, n_records, n_species)

# Plots with two censuses
census_all %>%
  filter(n_years > 1) %>%
  arrange(stratum, selection, plot_id, Year) %>%
  print(n = Inf, width = Inf)

census_all %>% distinct(stratum, plot_id, selection) %>% count(stratum, selection)

# Documented year without records in a layer (both censuses are then used)
census_all %>%
  filter(n_years > 1, !is.na(documented_year), selection == "two censuses averaged") %>%
  distinct(stratum, plot_id, documented_year)

present_census_selected <- census_all %>%
  filter(used) %>%
  select(stratum, plot_id, year = Year, selection)

present_census_year <- present_census_selected %>%
  group_by(stratum, plot_id, selection) %>%
  summarise(census_years = paste(sort(unique(year)), collapse = ", "),
            year = mean(year), .groups = "drop")

write_csv(present_census_selected, file.path(clean_dir, "present_census_selected.csv"))
write_csv(present_census_year,     file.path(clean_dir, "present_census_year.csv"))
write_csv(census_all, file.path(clean_dir, "check_present_census_completeness.csv"))
