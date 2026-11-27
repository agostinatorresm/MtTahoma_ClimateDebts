# -------------------------------------------------------------------------
# Author: Agostina Torres
# Project: MORAclimshift
# Script: 05_build_overstory.R
# Description:
# Builds long- and wide-format overstory (Tree) tables, past (1975) and
# present combined, following the basal-area method of the earlier
# cleaning script (CLIMreassembly/MORAcleaning/scripts/clean_overstory.R),
# applied to the resolved plots (plot_id_list_tree.csv) and species
# (spp_lib_current/spp_lib_legacy, canopy_broad == "overstory").
#
# Basal area (as in the original script):
#   basal_area_class_cm2 = pi * (diam_class / 2)^2
#   basal_area_class_m2  = basal_area_class_cm2 / 10000
#   basal_area_class_weighted = basal_area_class_m2 * num_ind
# diam_class is the midpoint of each 10-cm DBH class (e.g. DBH0 -> 5 cm),
# not the exact DBH of the individual. This applies to both periods:
# present-day DBH values are binned into the same DBH0-DBH120 classes as in
# 1975 before computing basal area, so that both surveys follow the same
# method (Step 6).
#
# PRESENT-DAY CENSUS: 13 plots were surveyed in two years. Where the data
# provider documented which census to keep, only that census is used; for the
# remaining plots both censuses are used and basal area and count are
# averaged across the census years in which the species was recorded
# (present_census_selected.csv, 04b_select_present_census.R). The present-day
# year is the census year, or the mean of the two years for averaged plots
# (present_census_year.csv).
#
# Outputs (pipeline/01_data_prep/01_cleaning/outputs/):
#   overstory_long.csv       - one row per plot x species x period, with
#                               occurrence/count/basal_area_m2/
#                               basal_area_perc/density
#   overstory_wide_occ.csv   - presence/absence matrix (plot x period rows,
#                               species columns)
#   overstory_wide_count.csv - abundance (individual count) matrix
#   overstory_wide_ba_m2.csv - basal area (m2) matrix
#   overstory_wide_ba_perc.csv - basal area (% of plot area) matrix
#   overstory_ageclass.csv   - individual count by plot x species x DBH
#                               class (diagnostic / diameter-structure use)
# -------------------------------------------------------------------------

library(tidyverse)
library(here)

outputs_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
tmp_dir     <- here("pipeline", "01_data_prep", "01_cleaning", "tmp_files")


# ==============================================================================
# STEP 1 — Load cleaned data
# ==============================================================================

load(file.path(outputs_dir, "tree_veg_plotid_clean.rda"))
# brings in: tree_R_clean, tree_L_clean, veg_R_clean, veg_L_clean

# Present-day census records used per plot (04b_select_present_census.R):
# the documented census year, or both censuses (averaged below) for plots
# surveyed twice without a documented choice. present_year gives one
# present-day year per plot (mean of the census years for averaged plots).
present_census <- read_csv(file.path(outputs_dir, "present_census_selected.csv"),
                           show_col_types = FALSE) |>
  filter(stratum == "overstory") |>
  select(plot_id, Year = year)
present_year <- read_csv(file.path(outputs_dir, "present_census_year.csv"),
                         show_col_types = FALSE) |>
  filter(stratum == "overstory") |>
  select(plot_id, year)
n_rows_before <- nrow(tree_R_clean)
tree_R_clean <- tree_R_clean |> semi_join(present_census, by = c("plot_id", "Year"))
cat("Present overstory records kept:", nrow(tree_R_clean),
    "of", n_rows_before, "\n")

spp_lib_current <- read_csv(file.path(outputs_dir, "spp_lib_current.csv"), show_col_types = FALSE)
spp_lib_legacy  <- read_csv(file.path(outputs_dir, "spp_lib_legacy.csv"),  show_col_types = FALSE)
plot_id_list_tree <- read_csv(file.path(outputs_dir, "plot_id_list_tree.csv"), show_col_types = FALSE)

overstory_current <- spp_lib_current |> filter(canopy_broad == "overstory") |> select(prst_code, species = final_species)
overstory_legacy  <- spp_lib_legacy  |> filter(final_species %in% overstory_current$species) |> select(past_code, species = final_species)

# Same code-typo corrections as in 04_apply_spp_lib.R (repeated here
# because this script runs independently of that one).
extract_code <- function(x) str_squish(str_extract(x, "^[^-]+"))

species_name_corrections <- tribble(
  ~species_code_raw, ~corrected_code, ~reason,
  "OLPHOR",  "OPLHOR", "letters transposed (L<->P)",
  "PACMYR",  "PAXMYR", "letter C->X",
  "PAIMYR",  "PAXMYR", "letter I->X",
  "TSHUET",  "TSUHET", "letters transposed (H<->U)",
  "TSUMEN",  "TSUMER", "last letter N->R",
  "VACMEN",  "VACMEM", "last letter N->M",
  "VIOGIA",  "VIOGLA", "letters transposed (I<->L)",
  "POLMAN",  "POLMUN", "letter A->U",
  "PSUMEN",  "PSEMEN", "letters transposed (U<->E)",
  "THUPLE",  "THUPLI", "last letter E->I",
  "BELSPI",  "BLESPI", "letters transposed (E<->L)",
  "BLESPE",  "BLESPI", "last letter E->I",
  "MEMFER",  "MENFER", "letters transposed (M<->N)",
  "MENFEK",  "MENFER", "letter K->R",
  "MENFIR",  "MENFER", "letter I->E",
  "TEATRI",  "TIATRI", "letters transposed (E<->I)",
  "TIATARI", "TIATRI", "extra letter A",
  "TRIATRI", "TIATRI", "letters transposed at start",
  "ACTTRI",  "ACHTRI", "letters transposed (C<->H)",
  "CHIMAP",  "CHIMEN", "truncated/altered ending",
  "CLIUMI",  "CLIUNI", "letter M->N",
  "HEIALB",  "HIEALB", "letters transposed (E<->I)",
  "LUMMAR",  "LOMMAR", "letter U->O",
  "MUEMAC",  "MOEMAC", "letter U->O",
  "PETPAL",  "PETFRI", "different ending altogether (name updated)",
  "STRALP",  "STRAMP", "letters transposed (L<->M... P kept)",
  "STRALT",  "STRAMP", "different ending altogether",
  "LISCAU",  "NEOBAN", "genus reclassified Listera->Neottia; caurina further merged into banksiana",
  "LISCON",  "NEOCON", "genus reclassified Listera->Neottia (convallarioides)",
  "LISCOR",  "NEOCOR", "genus reclassified Listera->Neottia (cordata)"
)


# ==============================================================================
# STEP 2 — Past | Overstory | Clean and resolve species
# ==============================================================================

past_over_clean <- tree_L_clean |>
  filter(plot_id %in% plot_id_list_tree$plot_id) |>
  mutate(species_code = str_squish(Species)) |>
  left_join(species_name_corrections, by = c("species_code" = "species_code_raw")) |>
  mutate(species_code_lookup = coalesce(corrected_code, species_code)) |>
  left_join(overstory_legacy, by = c("species_code_lookup" = "past_code")) |>
  filter(!is.na(species)) |>
  rename(DBH120 = `>DBH120`)

cat("Past overstory - plots:", n_distinct(past_over_clean$plot_id), "\n")


# ==============================================================================
# STEP 3 — Past | Reshape DBH classes and bigtrees to long format
# ==============================================================================

past_over_dbh <- past_over_clean |>
  select(plot_id, species, Plot_area_m2, starts_with("DBH")) |>
  pivot_longer(cols = starts_with("DBH"), names_to = "class_size", values_to = "num_ind") |>
  mutate(
    diam_class = case_when(
      str_starts(class_size, "DBH") ~ as.numeric(str_extract(class_size, "\\d+")) + 5,
      TRUE ~ as.numeric(class_size)
    )
  )

past_over_bigtree <- past_over_clean |>
  select(plot_id, species, starts_with("bigtree")) |>
  pivot_longer(cols = starts_with("bigtree"), names_to = "bigtree_index", values_to = "class_size") |>
  filter(!is.na(class_size)) |>
  mutate(
    class_size = as.character(as.numeric(class_size)),
    num_ind    = 1,
    diam_class = as.numeric(class_size)
  ) |>
  select(-bigtree_index)

past_over_cleaner <- bind_rows(past_over_dbh, past_over_bigtree)


# ==============================================================================
# STEP 4 — Past | Basal area estimation
# ==============================================================================
# The DBH120 class is excluded to avoid double-counting with bigtrees,
# where individuals >= 120 cm are recorded with their diameter.

past_over_cov <- past_over_cleaner |>
  filter(class_size != "DBH120") |>
  mutate(num_ind = replace_na(num_ind, 0)) |>
  mutate(
    basal_area_class_m2       = (pi * (diam_class / 2)^2) / 10000,
    basal_area_class_weighted = basal_area_class_m2 * num_ind
  )

past_over_cov_spp <- past_over_cov |>
  group_by(species, plot_id) |>
  summarise(basal_area_spp_m2 = sum(basal_area_class_weighted, na.rm = TRUE), .groups = "drop")

# Plot area: the most common Plot_area_m2 value per plot (some plots have
# inconsistent values across rows in the raw data).
plot_area_fix <- past_over_clean |>
  filter(!is.na(Plot_area_m2)) |>
  count(plot_id, Plot_area_m2) |>
  group_by(plot_id) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  select(plot_id, plot_area = Plot_area_m2)

past_over_cov_spp <- past_over_cov_spp |>
  left_join(plot_area_fix, by = "plot_id") |>
  mutate(basal_area_spp_perc = (basal_area_spp_m2 / plot_area) * 100,
         period = "past", year = 1975)

past_over_count_spp <- past_over_cov |>
  group_by(plot_id, species) |>
  summarise(count = sum(num_ind, na.rm = TRUE), .groups = "drop") |>
  mutate(period = "past", year = 1975)


# ==============================================================================
# STEP 5 — Present | Overstory | Clean, resolve species, clean DBH
# ==============================================================================

prst_over_clean <- tree_R_clean |>
  filter(plot_id %in% plot_id_list_tree$plot_id) |>
  mutate(species_code = extract_code(Species)) |>
  left_join(species_name_corrections, by = c("species_code" = "species_code_raw")) |>
  mutate(species_code_lookup = coalesce(corrected_code, species_code)) |>
  left_join(overstory_current, by = c("species_code_lookup" = "prst_code")) |>
  filter(!is.na(species)) |>
  # Data-entry typo documented in the earlier cleaning script: SMO-366's
  # Abies amabilis DBH "1703" is 170.3.
  mutate(DBH = if_else(plot_id == "SMO-366" & species == "Abies amabilis" & DBH == "1703",
                       "170.3", DBH)) |>
  mutate(
    DBH_clean = str_replace_all(DBH, "[^0-9.<]", ""),
    DBH_clean = case_when(
      DBH_clean %in% c("<0.5")       ~ "0.25",
      DBH_clean %in% c("<1", "<1.0") ~ "0.5",
      DBH_clean == ""                ~ NA_character_,
      TRUE                           ~ DBH_clean
    ),
    DBH = as.numeric(DBH_clean)
  ) |>
  select(-DBH_clean) |>
  filter(!is.na(DBH))

cat("Present overstory - plots:", n_distinct(prst_over_clean$plot_id), "\n")


# ==============================================================================
# STEP 6 — Present | Bin DBH into the same 10-cm classes as the past
# survey and estimate basal area per class
# ==============================================================================

prst_over_dbh <- prst_over_clean |>
  mutate(
    class_size = case_when(
      DBH > 120 ~ as.character(DBH),
      TRUE ~ as.character(
        cut(DBH, breaks = seq(0, 130, by = 10),
            labels = paste0("DBH", seq(0, 120, by = 10)),
            include.lowest = TRUE, right = FALSE)
      )
    ),
    diam_class = case_when(
      str_starts(class_size, "DBH") ~ as.numeric(str_extract(class_size, "\\d+")) + 5,
      TRUE ~ as.numeric(class_size)
    )
  )

prst_over_cov <- prst_over_dbh |>
  group_by(plot_id, species, Year, class_size, diam_class) |>
  summarise(num_ind = n(), .groups = "drop") |>
  mutate(
    basal_area_class_m2       = (pi * (diam_class / 2)^2) / 10000,
    basal_area_class_weighted = basal_area_class_m2 * num_ind
  )

prst_over_cov_spp_year <- prst_over_cov |>
  group_by(Year, plot_id, species) |>
  summarise(basal_area_spp_m2 = sum(basal_area_class_weighted, na.rm = TRUE), .groups = "drop")

cat("Present resurvey years found:\n")
print(prst_over_cov_spp_year |> count(Year, name = "n_plot_species_rows"))

# Basal area and count averaged across the census years used for a plot
# (one year except for plots averaged in STEP 1), over the years in which the
# species was recorded.
prst_over_cov_spp <- prst_over_cov_spp_year |>
  group_by(plot_id, species) |>
  summarise(basal_area_spp_m2 = mean(basal_area_spp_m2, na.rm = TRUE), .groups = "drop") |>
  mutate(
    plot_area           = 500,  # present-day plot area (no plot-area field in the resurvey data)
    basal_area_spp_perc = (basal_area_spp_m2 / plot_area) * 100,
    period = "present"
  ) |>
  left_join(present_year, by = "plot_id")

prst_over_count_spp_year <- prst_over_cov |>
  group_by(Year, plot_id, species) |>
  summarise(count = sum(num_ind, na.rm = TRUE), .groups = "drop")

prst_over_count_spp <- prst_over_count_spp_year |>
  group_by(plot_id, species) |>
  summarise(count = mean(count, na.rm = TRUE), .groups = "drop") |>
  mutate(period = "present") |>
  left_join(present_year, by = "plot_id")


# ==============================================================================
# STEP 7 — Combine into long format
# ==============================================================================

past_over_long <- past_over_cov_spp |>
  left_join(past_over_count_spp |> select(plot_id, species, count), by = c("plot_id", "species")) |>
  mutate(occurrence = if_else(count > 0, 1, 0), density = count / plot_area) |>
  select(period, year, plot_id, plot_area, species, occurrence, count,
         basal_area_spp_m2, basal_area_spp_perc, density)

prst_over_long <- prst_over_cov_spp |>
  left_join(prst_over_count_spp |> select(plot_id, species, count), by = c("plot_id", "species")) |>
  mutate(occurrence = if_else(count > 0, 1, 0), density = count / plot_area) |>
  select(period, year, plot_id, plot_area, species, occurrence, count,
         basal_area_spp_m2, basal_area_spp_perc, density)

overstory_long <- bind_rows(past_over_long, prst_over_long)

# ==============================================================================
# STEP 7b — Restrict to the confirmed plots and species:
#   - plot_id_list_tree.csv keeps the confirmed plots with Tree data in
#     both periods, but for any tree species, not only overstory species.
#     Plots without overstory occurrence in one period are removed here.
#   - Species with no occurrence across all plots and periods are also
#     removed (not expected after pivot_wider, but checked explicitly).
# ==============================================================================

plots_with_occurrence_by_period <- overstory_long |>
  filter(occurrence == 1) |>
  distinct(plot_id, period)

# Saved for 02_explore_data - shows the real distribution of plots with
# overstory occurrence only in past, only in present, or both, BEFORE the
# both-periods restriction below collapses this to a single category.
overstory_plot_coverage <- plots_with_occurrence_by_period |>
  mutate(has = TRUE) |>
  pivot_wider(names_from = period, values_from = has, values_fill = FALSE) |>
  mutate(coverage = case_when(
    past & present ~ "both",
    past            ~ "only_past",
    present         ~ "only_present"
  )) |>
  select(plot_id, coverage)
write_csv(overstory_plot_coverage, file.path(outputs_dir, "overstory_plot_coverage.csv"))
cat("\nOverstory plot coverage by period (before both-periods filter):\n")
print(overstory_plot_coverage |> count(coverage))

plots_with_both_periods <- plots_with_occurrence_by_period |>
  count(plot_id) |>
  filter(n == 2) |>
  pull(plot_id)

species_with_any_occurrence <- overstory_long |>
  filter(occurrence == 1) |>
  distinct(species) |>
  pull(species)

n_plots_before <- n_distinct(overstory_long$plot_id)
n_species_before <- n_distinct(overstory_long$species)

overstory_long <- overstory_long |>
  filter(plot_id %in% plots_with_both_periods, species %in% species_with_any_occurrence)

cat("\nScoped to plots/species with real occurrence in both periods:\n")
cat("  Plots  :", n_plots_before, "->", n_distinct(overstory_long$plot_id), "\n")
cat("  Species:", n_species_before, "->", n_distinct(overstory_long$species), "\n")

write_csv(overstory_long, file.path(outputs_dir, "overstory_long.csv"))
cat("\nSaved overstory_long.csv:", nrow(overstory_long), "rows\n")
cat("  Plots in past long   :", n_distinct(past_over_long$plot_id), "\n")
cat("  Plots in present long:", n_distinct(prst_over_long$plot_id), "\n")
cat("  Plot sets equal      :", setequal(unique(overstory_long$plot_id[overstory_long$period == "past"]),
                                         unique(overstory_long$plot_id[overstory_long$period == "present"])), "\n")


# ==============================================================================
# STEP 8 — Wide-format species x plot matrices, one per period, then
# combined (rows = plot_id x period)
# ==============================================================================

make_matrix <- function(df, value_col) {
  df |>
    select(period, year, plot_id, species, all_of(value_col)) |>
    pivot_wider(names_from = species, values_from = all_of(value_col), values_fill = 0)
}

overstory_wide_occ      <- make_matrix(overstory_long, "occurrence")
overstory_wide_count    <- make_matrix(overstory_long, "count")
overstory_wide_ba_m2    <- make_matrix(overstory_long, "basal_area_spp_m2")
overstory_wide_ba_perc  <- make_matrix(overstory_long, "basal_area_spp_perc")

write_csv(overstory_wide_occ,     file.path(outputs_dir, "overstory_wide_occ.csv"))
write_csv(overstory_wide_count,   file.path(outputs_dir, "overstory_wide_count.csv"))
write_csv(overstory_wide_ba_m2,   file.path(outputs_dir, "overstory_wide_ba_m2.csv"))
write_csv(overstory_wide_ba_perc, file.path(outputs_dir, "overstory_wide_ba_perc.csv"))

cat("\nSaved wide-format matrices:\n")
cat("  overstory_wide_occ     :", nrow(overstory_wide_occ), "rows x", ncol(overstory_wide_occ) - 3, "species\n")
cat("  overstory_wide_ba_m2   :", nrow(overstory_wide_ba_m2), "rows x", ncol(overstory_wide_ba_m2) - 3, "species\n")


# ==============================================================================
# STEP 9 — Abundance by diameter class (diagnostic / diameter-structure use)
# ==============================================================================

past_ageclass <- past_over_cov |>
  group_by(plot_id, species, class_size) |>
  summarise(num_ind = sum(num_ind, na.rm = TRUE), .groups = "drop") |>
  mutate(period = "past", year = 1975)

prst_ageclass <- prst_over_cov |>
  group_by(plot_id, species, class_size) |>
  summarise(num_ind = sum(num_ind, na.rm = TRUE), .groups = "drop") |>
  mutate(period = "present") |>
  left_join(present_year, by = "plot_id")

overstory_ageclass <- bind_rows(past_ageclass, prst_ageclass)
write_csv(overstory_ageclass, file.path(outputs_dir, "overstory_ageclass.csv"))
cat("\nSaved overstory_ageclass.csv:", nrow(overstory_ageclass), "rows\n")