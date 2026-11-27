# -------------------------------------------------------------------------
# Author: Agostina Torres
# Project: MORAclimshift
# Script: 02_apply_plot_ids.R
# Description:
# Applies the plot_id cleaning of 01_clean_plot_ids.R to the Tree and Veg
# data (Legacy + Resurvey), for the 332 plots confirmed in both periods
# (see plot_lib.csv / plot_id_list.csv).
#
# For each dataset:
#   1. the raw plot id is replaced by the clean plot_id (via plot_lib)
#   2. rows with different raw spellings of the same plot are grouped under
#      the same clean plot_id
#   3. only plots confirmed in both periods (plot_id_list) are kept
#
# Tree and Veg only; seedlings are not part of the analyses.
# -------------------------------------------------------------------------

library(tidyverse)
library(here)

org_dir     <- here("pipeline", "01_data_prep", "01_cleaning", "data", "data_raw", "01_raw_surveys")
outputs_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")


# ==============================================================================
# STEP 1 — Load the plot_lib mapping (built by 01_clean_plot_ids.R)
# ==============================================================================

plot_lib <- read_csv(file.path(outputs_dir, "plot_lib.csv"), show_col_types = FALSE)
plot_id_list_tree <- read_csv(file.path(outputs_dir, "plot_id_list_tree.csv"), show_col_types = FALSE)
plot_id_list_veg  <- read_csv(file.path(outputs_dir, "plot_id_list_veg.csv"),  show_col_types = FALSE)

# Build simple raw -> clean lookup tables, one per dataset x period,
# straight from plot_lib's 4 org columns.
build_lookup <- function(org_col) {
  plot_lib |>
    filter(!is.na(.data[[org_col]])) |>
    distinct(plot_id_raw = .data[[org_col]], plot_id)
}

lookup_legacy_tree   <- build_lookup("plot_id_org_legacy_tree")
lookup_legacy_veg    <- build_lookup("plot_id_org_legacy_veg")
lookup_resurvey_tree <- build_lookup("plot_id_org_resurvey_tree")
lookup_resurvey_veg  <- build_lookup("plot_id_org_resurvey_veg")


# ==============================================================================
# STEP 2 — Load raw data
# ==============================================================================

tree_R <- read_csv(file.path(org_dir, "MORA_TreeData_Resurvey_Org.csv"), show_col_types = FALSE)
tree_L <- read_csv(file.path(org_dir, "MORA_TreeData_Legacy_Org.csv"),   show_col_types = FALSE)
veg_R  <- read_csv(file.path(org_dir, "MORA_VegData_Resurvey_Org.csv"),  show_col_types = FALSE)
veg_L  <- read_csv(file.path(org_dir, "MORA_VegData_Legacy_Org.csv"),    show_col_types = FALSE)

# Tree Legacy has no SiteName - it is taken from Veg Legacy (same crosswalk
# as in 01_clean_plot_ids.R; all 463 Tree Legacy Plotno values exist in
# Veg Legacy).
site_lookup <- veg_L |> distinct(Plotno, SiteName) |> filter(!is.na(Plotno), !is.na(SiteName))
tree_L <- tree_L |>
  left_join(site_lookup, by = "Plotno") |>
  mutate(plot_id_raw = str_glue("{SiteName}-{Plotno}") |> as.character())

veg_L <- veg_L |> mutate(plot_id_raw = str_glue("{SiteName}-{Plotno}") |> as.character())


# ==============================================================================
# STEP 3 — Replace raw plot id with clean plot_id, keep only confirmed plots
# ==============================================================================

replace_and_filter <- function(df, raw_col, lookup, valid_ids) {
  n_before <- nrow(df)
  result <- df |>
    left_join(lookup, by = setNames("plot_id_raw", raw_col)) |>
    filter(plot_id %in% valid_ids$plot_id)
  cat(str_glue("  {n_before} rows -> {nrow(result)} rows after replacing + filtering to confirmed plots\n"))
  result
}

cat("Tree Resurvey:\n")
tree_R_clean <- replace_and_filter(tree_R, "Plot ID", lookup_resurvey_tree, plot_id_list_tree)

cat("Tree Legacy:\n")
tree_L_clean <- replace_and_filter(tree_L, "plot_id_raw", lookup_legacy_tree, plot_id_list_tree)

cat("Veg Resurvey:\n")
veg_R_clean <- replace_and_filter(veg_R, "Plot ID", lookup_resurvey_veg, plot_id_list_veg)

cat("Veg Legacy:\n")
veg_L_clean <- replace_and_filter(veg_L, "plot_id_raw", lookup_legacy_veg, plot_id_list_veg)


# ==============================================================================
# STEP 4 — Check: number of unique clean plot_ids in each dataset
# ==============================================================================
# This step does not aggregate species measurements across rows with
# different raw spellings; it only checks that the same plot_id is used
# wherever it belongs. How values are combined per plot depends on the
# variable (sum of cover, mean DBH, etc.) and is done in the scripts that
# build each layer.

cat("\nUnique plot_id per cleaned dataframe:\n")
cat("  tree_R_clean:", n_distinct(tree_R_clean$plot_id), "\n")
cat("  tree_L_clean:", n_distinct(tree_L_clean$plot_id), "\n")
cat("  veg_R_clean :", n_distinct(veg_R_clean$plot_id), "\n")
cat("  veg_L_clean :", n_distinct(veg_L_clean$plot_id), "\n")

save(tree_R_clean, tree_L_clean, veg_R_clean, veg_L_clean,
     file = file.path(outputs_dir, "tree_veg_plotid_clean.rda"))
cat("\nSaved tree_veg_plotid_clean.rda\n")


# ==============================================================================
# STEP 5 — Check that each confirmed plot has at least one species record
# in both periods (a matching plot_id alone does not guarantee that the
# plot has species rows)
# ==============================================================================
# tree_L, tree_R, veg_R are long format (one row per species observation,
# "Species" column). veg_L is WIDE format (one row per plot, one column per
# species) - "has a species" there means at least one species column is > 0.
#
# Placeholder codes that are not a species (e.g. "NOTREE-No trees" in
# seed_R) are excluded here for tree_L/tree_R/veg_R.

placeholder_species <- c("NOTREE-No trees", "NONE", "NA-nothing found")

plots_with_species_long <- function(df, plot_col = "plot_id") {
  df |>
    filter(!is.na(Species), !(Species %in% placeholder_species)) |>
    distinct(.data[[plot_col]]) |>
    pull(1)
}

plots_with_species_wide <- function(df, plot_col = "plot_id", first_species_col, last_species_col) {
  species_cols <- names(df)[which(names(df) == first_species_col):which(names(df) == last_species_col)]
  df |>
    mutate(has_species = rowSums(across(all_of(species_cols)), na.rm = TRUE) > 0) |>
    filter(has_species) |>
    distinct(.data[[plot_col]]) |>
    pull(1)
}

tree_R_has_spp <- plots_with_species_long(tree_R_clean)
tree_L_has_spp <- plots_with_species_long(tree_L_clean)
veg_R_has_spp  <- plots_with_species_long(veg_R_clean)
veg_L_has_spp  <- plots_with_species_wide(veg_L_clean, first_species_col = "ABAM1", last_species_col = "VIOLA")

cat("\nPlots with >=1 real species record:\n")
cat("  Tree present:", length(tree_R_has_spp), "/", n_distinct(tree_R_clean$plot_id), "\n")
cat("  Tree past   :", length(tree_L_has_spp), "/", n_distinct(tree_L_clean$plot_id), "\n")
cat("  Veg present :", length(veg_R_has_spp),  "/", n_distinct(veg_R_clean$plot_id),  "\n")
cat("  Veg past    :", length(veg_L_has_spp),  "/", n_distinct(veg_L_clean$plot_id),  "\n")

# Plots confirmed in both periods and with species data in both
tree_confirmed_with_spp <- intersect(tree_R_has_spp, tree_L_has_spp)
veg_confirmed_with_spp  <- intersect(veg_R_has_spp,  veg_L_has_spp)

cat("\nTree: plots with real species data in BOTH periods:", length(tree_confirmed_with_spp), "\n")
cat("Veg : plots with real species data in BOTH periods:", length(veg_confirmed_with_spp), "\n")

# Plots that match by plot_id but have no species in one period are
# flagged for review.
tree_empty_one_side <- setdiff(union(tree_R_has_spp, tree_L_has_spp), tree_confirmed_with_spp)
veg_empty_one_side  <- setdiff(union(veg_R_has_spp,  veg_L_has_spp),  veg_confirmed_with_spp)

cat("\nTree plots confirmed by plot_id but missing real species in one period:",
    length(tree_empty_one_side), "\n")
print(tree_empty_one_side)
cat("\nVeg plots confirmed by plot_id but missing real species in one period:",
    length(veg_empty_one_side), "\n")
print(veg_empty_one_side)