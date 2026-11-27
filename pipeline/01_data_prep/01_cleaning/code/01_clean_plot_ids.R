# -------------------------------------------------------------------------
# Author: Agostina Torres
# Project: MORAclimshift
# Script: 01_clean_plot_ids.R
# Description:
# Cleans the raw plot identifiers of all survey files and compares a generic
# cleaning procedure (trim whitespace, uppercase, zero-pad) with the plot IDs
# of the earlier cleaning of the same surveys (see README), to identify where
# the generic procedure is sufficient and where individual review is needed.
#
# There is a single raw source, the _Org.csv survey files (the same files
# used in the earlier cleaning). MORA_PlotNames.csv
# (Plot_R/Plot_L) is the lookup table of the earlier cleaning, stored in the
# "Cleaned" folder; it is used only to look up the final plot ID of each raw
# value, not as a second raw source.
#
#   plot_id_raw -> raw plot ID as written in the survey data
#                   (tree/seedlings/veg x present/past)
#   plot_id_generic  -> generic cleaning of plot_id_raw
#   plot_id_ref  -> final plot ID of the earlier cleaning for the same raw
#                   value, looked up in MORA_PlotNames.csv
#
# Outputs (all in pipeline/01_data_prep/01_cleaning/outputs/):
#   plot_id_diagnostic.csv    - every unique raw plot ID, generic vs earlier cleaning
#   plot_id_final.csv         - simplified: veg_type | period | original | new | reasoning
#   plot_id_changed.csv       - plots whose ID changed
#   plot_id_needs_review.csv  - plots that could not be resolved
#   plot_lib.csv              - the 332 plots confirmed in both periods, with
#                                their original spelling in Legacy/Resurvey x
#                                Tree/Veg (one row per raw-spelling combination)
#   plot_id_list.csv          - the 332 confirmed plot_id values
#   plot_id_list_tree.csv     - the subset with Tree data in both periods
#   plot_id_list_veg.csv      - the subset with Veg data in both periods
# -------------------------------------------------------------------------

library(tidyverse)
library(here)

org_dir       <- here("pipeline", "01_data_prep", "01_cleaning", "data", "data_raw", "01_raw_surveys")
orig_cleaning_dir <- here("pipeline", "01_data_prep", "01_cleaning", "data", "data_raw", "02_earlier_cleaning")
outputs_dir   <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
tmp_dir       <- here("pipeline", "01_data_prep", "01_cleaning", "tmp_files")

dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir,     recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# STEP 1 — Gather unique raw plot IDs from all 6 raw datasets
# ==============================================================================
# tree, seedlings, veg x present (Resurvey) / past (Legacy) = 6 dataframes.
# Present uses free-text "Plot ID". Past uses "Plotno" (a plain number)
# combined with "SiteName" to get a comparable "SITE-###" spelling (e.g.
# "COUG-1"). Tree Legacy has no SiteName column of its own, so its SiteName
# is cross-referenced from Veg Legacy via the shared Plotno numbering
# (all 463 Tree Legacy Plotno values exist in Veg Legacy).

tree_R <- read_csv(file.path(org_dir, "MORA_TreeData_Resurvey_Org.csv"),  show_col_types = FALSE)
tree_L <- read_csv(file.path(org_dir, "MORA_TreeData_Legacy_Org.csv"),    show_col_types = FALSE)
seed_R <- read_csv(file.path(org_dir, "MORA_Seedlings_Resurvey_Org.csv"), show_col_types = FALSE)
seed_L <- read_csv(file.path(org_dir, "MORA_Seedlings_Legacy_Org.csv"),   show_col_types = FALSE)
veg_R  <- read_csv(file.path(org_dir, "MORA_VegData_Resurvey_Org.csv"),   show_col_types = FALSE)
veg_L  <- read_csv(file.path(org_dir, "MORA_VegData_Legacy_Org.csv"),     show_col_types = FALSE)

gather_present_ids <- function(df, dataset_label) {
  df |>
    distinct(plot_id_raw = `Plot ID`) |>
    filter(!is.na(plot_id_raw)) |>
    mutate(period = "present", source_dataset = dataset_label)
}

gather_past_ids <- function(df, dataset_label, siteName_lookup = NULL) {
  if ("SiteName" %in% names(df)) {
    # Seedlings/Veg Legacy: SiteName + Plotno, matching the "SITE-###"
    # format used elsewhere (consistent with Plot_L in MORA_PlotNames.csv,
    # e.g. "COUG-1" from SiteName = "COUG", Plotno = 1)
    df |>
      distinct(plot_id_raw = str_glue("{SiteName}-{Plotno}")) |>
      mutate(plot_id_raw = as.character(plot_id_raw)) |>
      filter(!is.na(plot_id_raw)) |>
      mutate(period = "past", source_dataset = dataset_label)
  } else {
    # Tree Legacy has no SiteName column, but its Plotno values follow the
    # same numbering as Veg/Seedlings Legacy (all 463 Tree Legacy Plotno
    # values exist in Veg Legacy). SiteName is taken from there so that Tree
    # Legacy also gets a "SITE-###" plot_id_raw instead of a bare number.
    df |>
      distinct(Plotno) |>
      left_join(siteName_lookup, by = "Plotno") |>
      mutate(plot_id_raw = str_glue("{SiteName}-{Plotno}")) |>
      mutate(plot_id_raw = as.character(plot_id_raw)) |>
      filter(!is.na(plot_id_raw)) |>
      distinct(plot_id_raw) |>
      mutate(period = "past", source_dataset = dataset_label)
  }
}

# Plotno -> SiteName crosswalk, sourced from Veg Legacy (has both columns)
siteName_lookup <- veg_L |>
  distinct(Plotno, SiteName) |>
  filter(!is.na(Plotno), !is.na(SiteName))

raw_ids <- bind_rows(
  gather_present_ids(tree_R, "tree"),
  gather_past_ids(tree_L,    "tree", siteName_lookup = siteName_lookup),
  gather_present_ids(seed_R, "seedlings"),
  gather_past_ids(seed_L,    "seedlings"),
  gather_present_ids(veg_R,  "veg"),
  gather_past_ids(veg_L,     "veg")
) |>
  distinct(plot_id_raw, period, source_dataset)

cat("Unique raw plot IDs found across all 6 datasets:", nrow(raw_ids), "\n")
cat(str_glue("  {sum(raw_ids$period == 'present')} present + {sum(raw_ids$period == 'past')} past\n"))


# ==============================================================================
# STEP 2 — Generic cleaning functions, applied to plot_id_raw -> plot_id_generic
# ==============================================================================

clean_whitespace <- function(x) str_squish(x) |> str_to_upper()

zero_pad_plot_number <- function(x) {
  site   <- str_extract(x, "^[^-]+(?=-\\d{1,2}$)")
  number <- str_extract(x, "(?<=-)\\d{1,2}$")
  padded <- str_glue("{site}-{str_pad(number, 3, pad = '0')}")
  if_else(!is.na(site) & !is.na(number), padded, x)
}

standardize_raw_plot_name <- function(x) x |> clean_whitespace() |> zero_pad_plot_number()

normalize_match_key <- function(x) {
  x |> str_to_upper() |> str_remove_all("\\s+") |> str_replace("-0+(\\d+)$", "-\\1")
}

describe_process <- function(raw) {
  trimmed <- clean_whitespace(raw)
  padded  <- zero_pad_plot_number(trimmed)
  steps <- c()
  if (str_squish(raw) != raw) steps <- c(steps, "trimmed whitespace")
  if (str_to_upper(raw) != raw) steps <- c(steps, "capitalized")
  if (padded != trimmed) steps <- c(steps, "zero-padded")
  if (length(steps) == 0) "no change needed" else paste(steps, collapse = " + ")
}

raw_ids <- raw_ids |>
  mutate(
    plot_id_generic  = standardize_raw_plot_name(plot_id_raw),
    process_AT  = map_chr(plot_id_raw, describe_process)
  )


# ==============================================================================
# STEP 3 — Pre-correction of typos in the raw values, before matching
# against MORA_PlotNames
# ==============================================================================
# 182 raw plot IDs across tree/seedlings/veg have no match in Plot_R/Plot_L,
# because in the earlier cleaning these typos were corrected on the raw
# data first and the corrected value was then matched. This table
# reproduces that pre-correction step. Named exceptions are typos that no
# generic rule can capture; consolidations are groups of raw variants that
# collapse to one plot.

named_exceptions <- tribble(
  ~plot_id_raw,       ~pre_corrected,       ~reason_KP,
  "RICK-314",         "RICK-374",           "RICK-314 does not exist; likely should be RICK-374",
  "SO-83",            "S 0 -83",            "letter O vs digit 0, not a spacing issue (the earlier-cleaning ID is not zero-padded)",
  "S O -083",         "S 0 -83",            "letter O vs digit 0, not a spacing issue (the earlier-cleaning ID is not zero-padded)",
  "MOWC-184",         "MOWC-184/MOWI-184",  "dual plot identity",
  "JUNE-478",         "JUNE-478-FIELD",     "2018-2020 seasons inconsistently added a -FIELD suffix (no clean rule)",
  "SUNR-235",         "SUNR-235-FIELD",     "same -FIELD inconsistency as JUNE-478",
  "PARA-479",         "PARA-479-FIELD",     "same -FIELD inconsistency as JUNE-478",
  "IPSU-470",         "IPSU-470-FIELD",     "same -FIELD inconsistency as JUNE-478",
  "NORT-546",         "NORT-456",           "there is not NORT-546; might not be a 546 at all",
  "ARO-219",          "AURO-219",           "misspelling",
  "NICK-314",         "NICK-318",           "there is no NICK-314; should be NICK-318",
  "RAMP-506FIELD-2",  "RAMP-506-FIELD-2",   "missing dash before FIELD-2",
  "MB-7-FILED",       "MB-7-FIELD",         "typo: FILED -> FIELD",
  "CHRI-031",         "CHRI-041",           "does not exist as CHRI-031; found via plot_id_diagnostic.csv mismatch, real correction to CHRI-041"
)

consolidations <- bind_rows(
  tibble(plot_id_raw = str_glue("RKPI-0{str_pad(1:60, 2, pad = '0')}"),
         pre_corrected = "RKPI-001",
         reason_KP = "data-entry mistake; all RKPI plots collapsed to RKPI-001"),
  tibble(plot_id_raw = c(str_glue("WHRI-00{1:9}"), str_glue("WHRI-0{10:45}")),
         pre_corrected = "WHRI-001",
         reason_KP = "data-entry mistake; all WHRI plots collapsed to WHRI-001"),
  tibble(plot_id_raw = str_glue("AURO-0{220:225}"),
         pre_corrected = "AURO-219",
         reason_KP = "data-entry mistake; collapsed to AURO-219"),
  tibble(plot_id_raw = str_glue("LOWE-{46:54}"),
         pre_corrected = "LOWE-45",
         reason_KP = "data-entry mistake; collapsed to LOWE-45"),
  tibble(plot_id_raw = c("RANG-136","RANH-136","RANG-137","RANH-137","RANH-138",
                         "RANG-139","RANH-139","RANG-140","RANH-140","RANG-141",
                         "RANH-141","RANG-142","RANH-142","RANG-143","RANH-143",
                         "RANG-144","RANH-144","RANG-145","RANH-145","RANG-146",
                         "RANH-146","RANG-147","RANH-147"),
         pre_corrected = "RANG-136",
         reason_KP = "RANG/RANH naming confusion; collapsed to RANG-136"),
  tibble(plot_id_raw = str_glue("MOW-{156:168}"),
         pre_corrected = "MOWI-156",
         reason_KP = "data-entry mistake; collapsed to MOWI-156")
) |>
  mutate(plot_id_raw = as.character(plot_id_raw), pre_corrected = as.character(pre_corrected))

pre_correction_lookup <- bind_rows(named_exceptions, consolidations) |>
  distinct(plot_id_raw, .keep_all = TRUE)

raw_ids <- raw_ids |>
  left_join(pre_correction_lookup, by = "plot_id_raw") |>
  mutate(
    # if a pre-correction exists, it is used (re-cleaned) as the lookup
    # value; otherwise the generic plot_id_generic is used
    lookup_value = coalesce(pre_corrected, plot_id_generic)
  )


# ==============================================================================
# STEP 4 — Look up the final plot ID of each raw value in MORA_PlotNames.csv
# ==============================================================================
# MORA_PlotNames.csv is the lookup table of the earlier cleaning, not a raw
# source. Plot_R is the key for present-day (Resurvey) spellings and Plot_L
# for Legacy spellings; both resolve to the same final ID, Plot_R_2. The
# lookup uses `lookup_value` (pre-corrected if needed, else plot_id_generic)
# rather than plot_id_raw (see Step 3).

plotNames <- read_csv(file.path(orig_cleaning_dir, "MORA_PlotNames.csv"), show_col_types = FALSE)

# plot_id_ref is cleaned (stray spaces removed, zero-padded) where it is
# created, so that every later use (the Step 5 diagnostic and the Step 9
# plot_info table) gets the same value. Plot_R_2 is not always formatted
# consistently: some final IDs contain a space ("L OH-275" instead of
# "LOH-275"). Without this cleanup, 7 plots (LOH-275, SPU-222, UCH-172,
# UOH-389, WFO-297, WFO-304, WFO-305) would not match the cleaned plot_id
# of plot_id_list and would be dropped from plot_info.
ref_lookup_present <- plotNames |>
  filter(!is.na(Plot_R_2), Plot_R_2 != "") |>
  distinct(Plot_R, Plot_R_2) |>
  rename(ref_lookup_key = Plot_R, plot_id_ref = Plot_R_2) |>
  mutate(
    plot_id_ref = zero_pad_plot_number(str_remove_all(plot_id_ref, "\\s+")),
    match_key  = normalize_match_key(ref_lookup_key)
  ) |>
  distinct(match_key, .keep_all = TRUE)

ref_lookup_past <- plotNames |>
  filter(!is.na(Plot_R_2), Plot_R_2 != "") |>
  distinct(Plot_L, Plot_R_2) |>
  rename(ref_lookup_key = Plot_L, plot_id_ref = Plot_R_2) |>
  mutate(
    plot_id_ref = zero_pad_plot_number(str_remove_all(plot_id_ref, "\\s+")),
    match_key  = normalize_match_key(ref_lookup_key)
  ) |>
  distinct(match_key, .keep_all = TRUE)

# Present rows are looked up in Plot_R; past rows (all three datasets have a
# SiteName-Plotno spelling) in Plot_L. Matching uses `lookup_value`
# (pre-corrected typo if one exists, else plot_id_generic); see Step 3.
present_diagnostic <- raw_ids |>
  filter(period == "present") |>
  mutate(match_key = normalize_match_key(lookup_value)) |>
  left_join(ref_lookup_present, by = "match_key")

past_diagnostic <- raw_ids |>
  filter(period == "past") |>
  mutate(match_key = normalize_match_key(lookup_value)) |>
  left_join(ref_lookup_past, by = "match_key")


# ==============================================================================
# STEP 5 — Combine, compare, and save
# ==============================================================================

plot_id_diagnostic <- bind_rows(present_diagnostic, past_diagnostic) |>
  mutate(
    # Final cleanup of the reference plot IDs: stray spaces are removed and
    # the standard zero-padding is applied, because Plot_R_2 in the lookup
    # table is not always formatted consistently (e.g. "S 0 -83" instead of
    # "S0-083"). Reference and cleaned IDs are normalised the same way so
    # that the final output is consistent.
    plot_id_ref = if_else(is.na(plot_id_ref), plot_id_ref,
                         zero_pad_plot_number(str_remove_all(plot_id_ref, "\\s+"))),
    # Error in the lookup table (not in the raw data): the raw survey data
    # of this plot use the letter "O" ("SO-83", "S O -083" in the Tree/Veg
    # Resurvey Org data), and its Name_L in MORA_PlotNames.csv ("S. O ")
    # also has a letter O, whereas Plot_R_2 has the digit "0". It is
    # corrected here because it affects the final identifier; the
    # raw-to-lookup matching still uses the literal "0" to find this row.
    plot_id_ref = if_else(plot_id_ref == "S0-083", "SO-083", plot_id_ref),
    # "match" = the final ID of the earlier cleaning was found (a bare
    # Legacy number like "45" never looks like "COUG-001" as text, so this
    # checks resolution, not textual similarity)
    match = !is.na(plot_id_ref),
    # "generic_only_match" = the generic cleaning alone (without manual
    # pre-correction) reached a resolved ID
    generic_only_match = is.na(pre_corrected) & match
  ) |>
  select(plot_id_raw, period, source_dataset, plot_id_generic, process_AT,
         pre_corrected, reason_KP, plot_id_ref, generic_only_match, match)

write_csv(plot_id_diagnostic, file.path(tmp_dir, "plot_id_diagnostic.csv"))

cat("\nSaved plot_id_diagnostic.csv:", nrow(plot_id_diagnostic), "rows\n")
cat("  Resolved total (match = TRUE)          :", sum(plot_id_diagnostic$match, na.rm = TRUE), "\n")
cat("    ...of which generic cleaning alone was enough:", sum(plot_id_diagnostic$generic_only_match, na.rm = TRUE), "\n")
cat("    ...of which needed a manual pre-correction    :",
    sum(plot_id_diagnostic$match & !plot_id_diagnostic$generic_only_match, na.rm = TRUE), "\n")
cat("  Not resolved (match = FALSE or no equivalent in the earlier cleaning):",
    sum(!plot_id_diagnostic$match | is.na(plot_id_diagnostic$match)), "\n")


# ==============================================================================
# STEP 6 — Final simplified table: veg_type | period | plot_id_original |
#          plot_id_new | reasoning
# ==============================================================================
# Final table, built from MORA_PlotNames.csv (the equivalences), the
# pre-corrections of the earlier cleaning script and their reasons, and
# the raw tree/seedlings/veg data (the original plot IDs).

plot_id_final <- plot_id_diagnostic |>
  mutate(
    veg_type = recode(source_dataset, "veg" = "understory", "seedlings" = "seedling", "tree" = "tree"),
    reasoning = case_when(
      !is.na(reason_KP) ~ reason_KP,
      match & !generic_only_match ~ "resolved via lookup, no generic change needed",
      match & generic_only_match & process_AT != "no change needed" ~ str_glue("generic cleaning: {process_AT}"),
      match ~ "no change needed",
      TRUE ~ "NOT RESOLVED - no equivalent found in MORA_PlotNames.csv, needs manual review"
    )
  ) |>
  select(veg_type, period, plot_id_original = plot_id_raw,
         plot_id_new = plot_id_ref, reasoning) |>
  mutate(
    # match: comparison of plot_id_original with plot_id_new (not a flag of
    # the earlier cleaning). TRUE = identical (no change), FALSE = the ID
    # was changed, by the generic cleaning or by the earlier cleaning.
    match = plot_id_original == plot_id_new
  )

write_csv(plot_id_final, file.path(tmp_dir, "plot_id_final.csv"))
cat("\nSaved plot_id_final.csv:", nrow(plot_id_final), "rows\n")
cat("  match = TRUE (no change) :", sum(plot_id_final$match, na.rm = TRUE), "\n")
cat("  match = FALSE (changed)  :", sum(!plot_id_final$match | is.na(plot_id_final$match)), "\n")

# Separate CSV: plots whose ID was changed
plot_id_changed <- plot_id_final |>
  filter(is.na(match) | !match)

write_csv(plot_id_changed, file.path(tmp_dir, "plot_id_changed.csv"))
cat("Saved plot_id_changed.csv:", nrow(plot_id_changed), "rows (only plots that changed)\n")


# ==============================================================================
# STEP 7 — Keep the unresolved cases for further review
# ==============================================================================

needs_review <- plot_id_diagnostic |>
  filter(is.na(match) | !match)

write_csv(needs_review, file.path(tmp_dir, "plot_id_needs_review.csv"))
cat("Saved plot_id_needs_review.csv:", nrow(needs_review), "rows to look at next\n")


# ==============================================================================
# STEP 8 — plot_lib: plots confirmed in both periods, with their original
# raw spelling in Legacy/Resurvey x Tree/Veg
# ==============================================================================
# "Confirmed in both periods" = Resurveyed == 'Y' in the lookup table, plus
# the 5 sites that the earlier cleaning script documents as mis-flagged
# "N" (they have a valid Plot_R_2 and Year despite Resurveyed == 'N', e.g.
# PARA-479, SOUT-223).
#
# When a plot has more than one raw spelling within the same dataset (e.g.
# Veg Legacy has "ARO-219" and "AURO-219" for the same plot), each variant
# gets its own row instead of being combined in one cell, so a plot can
# appear on several rows of plot_lib.

resurveyed_flag_corrections <- c("CHIN-247", "LOWE-425", "PARA-479", "SOUT-223", "U OH-389")

plots_in_both_periods <- plotNames |>
  filter(!is.na(Plot_R_2), Plot_R_2 != "") |>
  filter(str_to_upper(Resurveyed) == "Y" | Plot_R_2 %in% resurveyed_flag_corrections) |>
  mutate(plot_id = zero_pad_plot_number(str_remove_all(Plot_R_2, "\\s+"))) |>
  distinct(plot_id)

cat("\nPlots confirmed in BOTH periods:", nrow(plots_in_both_periods), "\n")

# One table of distinct variants per dataset x period (tree and veg;
# seedlings are not part of the analyses).
# plot_id_final uses `veg_type` (tree/seedling/understory), which replaced
# `source_dataset` in Step 6.
get_variants <- function(veg_type_label, per, colname) {
  plot_id_final |>
    filter(veg_type == veg_type_label, period == per, plot_id_new %in% plots_in_both_periods$plot_id) |>
    distinct(plot_id = plot_id_new, !!colname := plot_id_original)
}

legacy_tree    <- get_variants("tree",        "past",    "plot_id_org_legacy_tree")
legacy_veg     <- get_variants("understory",  "past",    "plot_id_org_legacy_veg")
resurvey_tree  <- get_variants("tree",        "present", "plot_id_org_resurvey_tree")
resurvey_veg   <- get_variants("understory",  "present", "plot_id_org_resurvey_veg")

# full_join on plot_id: when a plot has several variants in more than one
# table, this gives one row per combination (separate rows, not combined
# cells).
plot_lib <- plots_in_both_periods |>
  full_join(legacy_tree,   by = "plot_id", relationship = "many-to-many") |>
  full_join(legacy_veg,    by = "plot_id", relationship = "many-to-many") |>
  full_join(resurvey_tree, by = "plot_id", relationship = "many-to-many") |>
  full_join(resurvey_veg,  by = "plot_id", relationship = "many-to-many") |>
  arrange(plot_id)

write_csv(plot_lib, file.path(outputs_dir, "plot_lib.csv"))
cat("Saved plot_lib.csv:", nrow(plot_lib), "rows\n")

# Separate file: the unique plot_id list (all 332, regardless of dataset)
plot_id_list <- plots_in_both_periods |> arrange(plot_id)
write_csv(plot_id_list, file.path(outputs_dir, "plot_id_list.csv"))
cat("Saved plot_id_list.csv:", nrow(plot_id_list), "unique plot IDs\n")

# Tree/Veg lists: plots with data in both periods for that dataset (not
# every plot has both tree and understory measurements, e.g. above
# treeline).
plot_id_list_tree <- plot_lib |>
  filter(!is.na(plot_id_org_legacy_tree), !is.na(plot_id_org_resurvey_tree)) |>
  distinct(plot_id) |>
  arrange(plot_id)

plot_id_list_veg <- plot_lib |>
  filter(!is.na(plot_id_org_legacy_veg), !is.na(plot_id_org_resurvey_veg)) |>
  distinct(plot_id) |>
  arrange(plot_id)

write_csv(plot_id_list_tree, file.path(outputs_dir, "plot_id_list_tree.csv"))
write_csv(plot_id_list_veg,  file.path(outputs_dir, "plot_id_list_veg.csv"))
cat("Saved plot_id_list_tree.csv:", nrow(plot_id_list_tree), "plots with tree data in both periods\n")
cat("Saved plot_id_list_veg.csv:",  nrow(plot_id_list_veg),  "plots with veg data in both periods\n")


# ==============================================================================
# STEP 9 — plot_info: dedicated plot metadata (coordinates, elevation,
# aspect, slope, climate), one row per plot x period.
# ==============================================================================
# Source: MORA_Plot_Resurvey_Org.csv (present) and MORA_Plot_Legacy_Org.csv
# (past). MORA_PlotNames.csv has no environmental data, only the plot-name
# lookup used elsewhere in this script.
#
# Lat/Long (present only; Legacy has a single reading): two GPS units were
# used per plot (Lat_1/Long_1, Lat_2/Long_2). As in the earlier cleaning
# script (01_CleanMORA.R), the readings are not averaged: Lat_1/Long_1 is
# used, and Lat_2/Long_2 only when Lat_1 is missing.
#
# Data-entry corrections documented in the earlier cleaning script:
#   - one Lat value entered as text "46 .91759" -> 46.91759
#   - CHRI-097's Long_1 was 454.7823 (impossible) -> Long_2 is used
#   - ~20 plots had missing or >200 m-off coordinates -> manual corrections
#     (the list in plot_coordinate_fixes)
#
# Slope: collected in percent in 2014 and in degrees in all other years, but
# the 2014 values in the raw files are already in degrees: each one equals
# DEGREES(ATAN(%/100)) of a one-decimal percent value, in both
# MORA_Cascade_Legacy_Data_Master2018.xlsx and MORA_Plot_Resurvey_Org.csv
# (00_check_master2018_vs_org.R). Slope is therefore used as recorded for
# all years; converting 2014 again would apply the transformation twice.

plot_R_raw <- read_csv(file.path(org_dir, "MORA_Plot_Resurvey_Org.csv"), show_col_types = FALSE)
plot_L_raw <- read_csv(file.path(org_dir, "MORA_Plot_Legacy_Org.csv"),   show_col_types = FALSE)

coord_fixes <- tribble(
  ~plot_id_raw,  ~lat_fixed,    ~long_fixed,
  "CARB-126",    46.995031,     121.897284,
  "MOUN-164",    46.941648,     121.889518,
  "REFL-237",    46.77044302,   121.740122,
  "TATO-326",    46.76872599,   121.740627,
  "U CH-172",    46.86413904,   121.533261,
  "WHIT-273",    46.89943104,   121.563778,
  "KAUT-355",    46.73650897,   121.853646,
  "ABOV-418",    46.79238299,   121.910259,
  "CORI-408",    46.98505,      121.85702,
  "IPSU-141",    46.98359801,   121.839552,
  "IPSU-471",    46.95128703,   121.874041,
  "MOWI-155",    46.94241813,   121.8953349,
  "MOWI-157",    46.93520903,   121.8674933,
  "MOWI-185",    46.90764003,   121.899112,
  "MOWI-187",    46.90248298,   121.906006,
  "MOWI-472",    46.94538684,   121.8719929,
  "PINN-372",    46.76692002,   121.732751,
  "RANG-138",    46.98211,      121.85316,
  "TAHO-003",    46.74158396,   121.899559,
  "TAHO-231",    46.74215301,   121.8926694
)

plot_R_coords <- plot_R_raw |>
  mutate(
    lat_1_num  = suppressWarnings(as.numeric(if_else(`Lat 1` == "46 .91759", "46.91759", `Lat 1`))),
    long_2_num = -abs(suppressWarnings(as.numeric(`Long 2`))),
    long_1_num = -abs(suppressWarnings(as.numeric(if_else(`Plot ID` == "CHRI-097", NA_character_, `Long 1`)))),
    lat_2_num  = suppressWarnings(as.numeric(`Lat 2`)),
    lat  = coalesce(lat_1_num, lat_2_num),
    long = coalesce(long_1_num, long_2_num),
    # Slope in degrees for all years (2014 values already converted in the
    # raw files; see note above)
    slope = suppressWarnings(as.numeric(Slope))
  ) |>
  left_join(coord_fixes, by = c("Plot ID" = "plot_id_raw")) |>
  mutate(lat = coalesce(lat_fixed, lat), long = coalesce(-abs(long_fixed), long)) |>
  transmute(
    plot_id_raw = `Plot ID`, period = "present", year = Year,
    order_entered = `Order entered`,
    lat, long, elevation = NA_real_,
    aspect = suppressWarnings(as.numeric(Aspect)), slope,
    ann_t = NA_real_, wint_ppt = NA_real_, sum_ppt = NA_real_, plotarea = NA_real_
  )

# Legacy: single Lat/Long/Elevation reading, no reconciliation needed.
# Plot_ID built the same way as elsewhere in this script (SiteName-Plotno).
# No Year column exists anywhere in the Legacy datasets (Plot/Tree/Veg/
# Seedlings) - the entire Legacy survey is the single 1975 campaign, so
# year is a constant here, not a per-row value like in Resurvey.
plot_L_coords <- plot_L_raw |>
  transmute(
    plot_id_raw = str_glue("{Plotno}"),  # matched via Num_L below, not SiteName-Plotno
    period = "past", year = 1975, order_entered = NA_real_,
    lat = Latitude, long = Longitude, elevation = Elevation,
    aspect = Aspect, slope = Slope,
    ann_t = Ann_T, wint_ppt = Wint_ppt, sum_ppt = Sum_ppt, plotarea = Plotarea
  )

# Look up clean plot_id: present via Plot_R, past via Num_L (Plot_L's own
# Plotno IS Num_L - same logic as Tree Legacy's crosswalk earlier).
plot_R_info <- plot_R_coords |>
  mutate(match_key = normalize_match_key(plot_id_raw)) |>
  left_join(ref_lookup_present, by = "match_key") |>
  filter(!is.na(plot_id_ref)) |>
  transmute(plot_id = plot_id_ref, period, year, order_entered, lat, long, elevation, aspect, slope,
            ann_t, wint_ppt, sum_ppt, plotarea)

ref_lookup_num_l <- plotNames |>
  filter(!is.na(Plot_R_2), Plot_R_2 != "") |>
  distinct(Num_L, Plot_R_2) |>
  mutate(
    Plot_R_2  = zero_pad_plot_number(str_remove_all(Plot_R_2, "\\s+")),
    match_key = as.character(suppressWarnings(as.numeric(Num_L)))
  ) |>
  filter(!is.na(match_key)) |>
  distinct(match_key, .keep_all = TRUE)

plot_L_info <- plot_L_coords |>
  mutate(match_key = as.character(suppressWarnings(as.numeric(plot_id_raw)))) |>
  left_join(ref_lookup_num_l, by = "match_key") |>
  filter(!is.na(Plot_R_2)) |>
  transmute(plot_id = Plot_R_2, period, year, order_entered, lat, long, elevation, aspect, slope,
            ann_t, wint_ppt, sum_ppt, plotarea)

plot_environment_by_year <- bind_rows(plot_R_info, plot_L_info) |>
  filter(plot_id %in% plot_id_list$plot_id) |>
  arrange(plot_id, period)

write_csv(plot_environment_by_year, file.path(outputs_dir, "plot_environment_by_year.csv"))
cat("\nSaved plot_environment_by_year.csv:", nrow(plot_environment_by_year), "rows (",
    n_distinct(plot_environment_by_year$plot_id), "unique plots) - full year-level detail,\n")
cat("  kept for a FUTURE past-vs-present-by-year comparison, not analyzed here.\n")

# plot_info: one row per plot x period. Past has one row per plot (Legacy
# has no multi-year issue). Present can have more than one row when a plot
# was resurveyed in more than one year (e.g. CHIN-084). These cases are
# flagged here and resolved by the documented decisions below rather than
# by a default rule such as keeping the most recent year.

present_multi_year <- plot_environment_by_year |>
  filter(period == "present") |>
  count(plot_id) |>
  filter(n > 1)

# Per-plot decisions documented in the earlier cleaning script
# (01_CleanMORA.R) for the 11 plots where the year to keep is stated. For
# the remaining ambiguous cases, COUG-001/LONG-044 (both rows have the same
# year, 2014) and RAMP-505/RAMP-506 (surveyed twice, with no decision in the
# original script), the plot metadata follow "Order entered" (the first
# row entered). The vegetation data of these plots follow
# 04b_select_present_census.R.
documented_year_choice <- tribble(
  ~plot_id,    ~keep_year,
  "CHIN-084",  2015,
  "GOAT-468",  2016,
  "LOWE-045",  2015,
  "RAMP-093",  2015,
  "SMO-364",   2019,
  "TAHO-003",  2014,
  "WHIT-071",  2014,
  "LONG-356",  2015,
  "WHIT-273",  2014,
  "WHIT-307",  2015
  # SUNS-046 excluded: the note says "keep 14" but the plot was surveyed in
  # 2015 and 2016 (no 2014 census), so it falls back to the order entered.
)

# Documented choices written out for 04b_select_present_census.R, which keeps
# only the documented census for these plots and both censuses otherwise
write_csv(documented_year_choice |> rename(year = keep_year),
          file.path(outputs_dir, "present_year_documented.csv"))

if (nrow(present_multi_year) > 0) {
  cat("\n", nrow(present_multi_year), "plot(s) had more than one present-year record.\n")
  cat("Resolved via the documented per-plot decision where one exists;\n")
  cat("Order entered (first entered) as fallback otherwise:\n")
}

plot_info_present <- plot_environment_by_year |>
  filter(period == "present") |>
  left_join(documented_year_choice, by = "plot_id") |>
  group_by(plot_id) |>
  mutate(n_rows = n()) |>
  filter(
    n_rows == 1 |
      (!is.na(keep_year) & year == keep_year) |
      (is.na(keep_year) & order_entered == min(order_entered, na.rm = TRUE))
  ) |>
  slice(1) |>  # safety net: if a rule still leaves >1 row, the first is kept
  ungroup() |>
  select(-keep_year, -n_rows)

if (nrow(present_multi_year) > 0) {
  print(
    plot_info_present |>
      filter(plot_id %in% present_multi_year$plot_id) |>
      select(plot_id, year, order_entered) |>
      arrange(plot_id)
  )
}

plot_info <- bind_rows(
  plot_environment_by_year |> filter(period == "past"),
  plot_info_present
) |>
  select(-order_entered)

# Complete: every plot in plot_id_list gets one row per period, even when
# the metadata sheet lacks that plot for one period (e.g. SUNR-235 has no
# row in MORA_Plot_Legacy_Org.csv, although the plot has 1975 Tree/Veg
# data in another file). Missing combinations get an explicit NA row.
plot_info <- plot_id_list |>
  expand_grid(period = c("past", "present")) |>
  left_join(plot_info, by = c("plot_id", "period")) |>
  arrange(plot_id, period)

missing_metadata <- plot_info |> filter(is.na(year))
if (nrow(missing_metadata) > 0) {
  cat("\nPlots with a period entirely missing from the source metadata (now an explicit NA row):\n")
  print(missing_metadata |> select(plot_id, period))
}

write_csv(plot_info, file.path(outputs_dir, "plot_info.csv"))
cat("\nSaved plot_info.csv:", nrow(plot_info), "rows (", n_distinct(plot_info$plot_id), "unique plots) -",
    "one row per plot x period\n")
cat("  Elevation available (past only - present has no source yet, e.g. a DEM):",
    sum(!is.na(plot_info$elevation)), "\n")
