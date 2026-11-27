# -------------------------------------------------------------------------
# Author: Agostina Torres
# Project: MORAclimshift
# Script: 06_build_understory.R
# Description:
# Builds long-format understory (Veg) tables, past (1975) and present
# combined. Present-day data were collected in 4 x 1-m2 subplots (QN/QE/QS/
# QW) plus a nested 50-m2 census of the plot ("ALL" records); past data are
# one % cover value per species per plot (variable plot size - see Plotarea
# in plot_info.csv). Cover is used as recorded in each period, without
# normalization to a common area; plot-size effects are addressed in the
# robustness module.
#
# PRESENT-DAY CENSUS: some plots were surveyed in two years. Where the data
# provider documented which census to keep, only that census is used; for the
# remaining plots both censuses are used (present_census_selected.csv,
# 04b_select_present_census.R) and present-day cover is the mean over all
# subplot x census combinations (absent = 0), i.e. the mean of the per-census
# subplot means; occurrence is 1 if the species was recorded in either
# census. The present-day year is the census year, or the mean of the two
# years for averaged plots (present_census_year.csv).
#
# SPECIES COMPLEXES: congeners not identified consistently between surveys
# are merged into one taxon (species_complexes, STEP 1), in both periods,
# before any aggregation.
#
# Occurrence: present-day occurrence is based on detection in a Q-subplot
# (QN/QE/QS/QW) or in "ALL", the nested 50-m2 presence/absence census;
# omitting "ALL" would understate presence for species recorded there but
# not in the four 1-m2 subplots. This applies to every present-day output of
# this script, so that occurrence is the same whatever the later use.
# Cover is based on the Q-subplots only: "ALL" has no % cover, so a species
# detected only via "ALL" gets cover = NA (not 0, which would mean
# "measured and absent"; see STEP 3a). Cover-weighted analyses (CWM,
# cover-weighted CTI/CMI) exclude NA-cover species from the weighting
# rather than imputing a value; occurrence-based metrics (richness,
# presence/absence) treat the species as present.
#
# Five outputs, because the present-day rows used differ by purpose:
#   understory_long.csv          - main comparable table: present
#                                   occurrence including "ALL" (see above);
#                                   present cover averaged across the 4
#                                   Q-subplots, NA for ALL-only detections;
#                                   the equivalent of overstory_long.csv.
#   understory_long_subplot.csv  - present-day subplot-level detail (one
#                                   row per QN/QE/QS/QW, plus "ALL" as its
#                                   own replicate row).
#   understory_presence.csv      - presence/absence only (no cover value;
#                                   detections only, not a full 0/1
#                                   matrix); present day includes "ALL"
#                                   together with the 4 subplots.
#   understory_occurrence_incl_all.csv / understory_wide_occ_incl_all.csv -
#                                   full 0/1 occurrence matrix, both
#                                   periods, not restricted to plots with
#                                   both periods (STEP 7); used for HMSC.
#   understory_cover_incl_all.csv / understory_wide_cov_incl_all.csv -
#                                   cover matrix, both periods, same scope
#                                   as above, NA for ALL-only present-day
#                                   detections.
#
# Cover-value conventions:
#   "T"/"t" (trace)        -> 0.5% cover. Occurs almost only in the
#                              Q-subplot records (QN/QE/QS/QW, ~1,860
#                              rows), so it enters the Q-subplot cover
#                              average of understory_long.csv.
#   "present"/"Present"    -> 0.5% cover (same value as trace). Occurs
#                              almost only in "ALL" records (~1,150 rows):
#                              "T" and "present" are the two conventions
#                              for a trace record, used in different
#                              replicate types.
#   "0"                    -> explicit zero cover (absence)
#   "4m circle"/blank/etc. -> not a subplot code, excluded from all outputs
#                              (not QN/QE/QS/QW/ALL)
#
# Outputs (pipeline/01_data_prep/01_cleaning/outputs/):
#   understory_long.csv, understory_long_subplot.csv, understory_presence.csv
#   understory_wide_cov.csv, understory_wide_occ.csv (from understory_long)
#
# Plus, in outputs/final/:
#   understory_occurrence_incl_all.csv, understory_wide_occ_incl_all.csv
#   understory_cover_incl_all.csv, understory_wide_cov_incl_all.csv
# -------------------------------------------------------------------------

library(tidyverse)
library(here)

outputs_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
tmp_dir     <- here("pipeline", "01_data_prep", "01_cleaning", "tmp_files")
final_dir   <- file.path(outputs_dir, "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)


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
  filter(stratum == "understory") |>
  select(plot_id, Year = year)
present_year <- read_csv(file.path(outputs_dir, "present_census_year.csv"),
                         show_col_types = FALSE) |>
  filter(stratum == "understory") |>
  select(plot_id, year)
n_rows_before <- nrow(veg_R_clean)
veg_R_clean <- veg_R_clean |> semi_join(present_census, by = c("plot_id", "Year"))
cat("Present understory records kept:", nrow(veg_R_clean),
    "of", n_rows_before, "\n")

# Species complexes: congeners whose records were not identified
# consistently between surveys (one species recorded in the first survey and
# the other in the second, in the same plots, far more often than the
# reverse). Records are merged into the retained taxon in both periods.
#   V. alaskaense -> V. ovalifolium : 22 plots ovalifolium (past) ->
#                                     alaskaense (present), 2 the reverse
#   V. orbiculata -> V. sempervirens: 38 plots sempervirens (past) ->
#                                     orbiculata (present), 0 the reverse
species_complexes <- c("Vaccinium alaskaense" = "Vaccinium ovalifolium",
                       "Viola orbiculata"     = "Viola sempervirens")

spp_lib_current <- read_csv(file.path(outputs_dir, "spp_lib_current.csv"), show_col_types = FALSE)
spp_lib_legacy  <- read_csv(file.path(outputs_dir, "spp_lib_legacy.csv"),  show_col_types = FALSE)
plot_id_list_veg <- read_csv(file.path(outputs_dir, "plot_id_list_veg.csv"), show_col_types = FALSE)
plot_info <- read_csv(file.path(outputs_dir, "plot_info.csv"), show_col_types = FALSE)

understory_current <- spp_lib_current |> filter(canopy_broad == "understory") |> select(prst_code, species = final_species)
understory_legacy  <- spp_lib_legacy  |> filter(final_species %in% understory_current$species) |> select(past_code, species = final_species)

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

# Cover-text -> numeric conversion, shared by the three present-day
# resolutions below.
clean_cover_value <- function(x) {
  x_clean <- str_squish(x)
  case_when(
    str_to_upper(x_clean) %in% c("T")                    ~ 0.5,
    str_to_lower(x_clean) %in% c("present")               ~ 0.5,
    TRUE ~ suppressWarnings(as.numeric(x_clean))
  )
}


# ==============================================================================
# STEP 2 — Past | Understory | Clean and resolve species (WIDE source)
# ==============================================================================

non_species_cols <- c("SiteName", "Plotno", "plot_id_raw", "plot_id")
veg_L_species_cols <- setdiff(names(veg_L_clean), non_species_cols)

past_under_long_raw <- veg_L_clean |>
  filter(plot_id %in% plot_id_list_veg$plot_id) |>
  select(plot_id, all_of(veg_L_species_cols)) |>
  pivot_longer(cols = -plot_id, names_to = "species_code_raw", values_to = "cover") |>
  mutate(species_code_lookup = str_squish(species_code_raw)) |>
  left_join(understory_legacy, by = c("species_code_lookup" = "past_code")) |>
  filter(!is.na(species)) |>
  mutate(species = recode(species, !!!species_complexes))

past_plot_area <- plot_info |>
  filter(period == "past") |>
  select(plot_id, plot_area = plotarea)

cat("Past understory - plots:", n_distinct(past_under_long_raw$plot_id), "\n")

past_under_for_long <- past_under_long_raw |>
  left_join(past_plot_area, by = "plot_id") |>
  # Several past_codes can resolve to the same final_species (including
  # merged species complexes); their cover is summed so that this object
  # has one row per plot_id x species.
  group_by(plot_id, plot_area, species) |>
  summarise(cover = sum(cover, na.rm = TRUE), .groups = "drop") |>
  mutate(period = "past", year = 1975,
         occurrence = if_else(!is.na(cover) & cover > 0, 1, 0))


# ==============================================================================
# STEP 3 — Present | Understory | Clean, resolve species, convert cover
# ==============================================================================

prst_under_clean <- veg_R_clean |>
  filter(plot_id %in% plot_id_list_veg$plot_id) |>
  mutate(species_code = extract_code(Species)) |>
  left_join(species_name_corrections, by = c("species_code" = "species_code_raw")) |>
  mutate(species_code_lookup = coalesce(corrected_code, species_code)) |>
  left_join(understory_current, by = c("species_code_lookup" = "prst_code")) |>
  filter(!is.na(species)) |>
  mutate(species = recode(species, !!!species_complexes)) |>
  mutate(Replicate = str_squish(Replicate)) |>
  mutate(cover = clean_cover_value(`Percent Cover`)) |>
  # "ALL" is presence/absence only and has no cover quantity, so a
  # trace-text record there ("Present"/"present", the only trace form in
  # ALL) is recoded from 0.5 (the trace value of clean_cover_value, meant
  # for % cover in the Q-subplots) to 1, a detection flag. This does not
  # change occurrence (both are > 0 in the all_detected check below).
  mutate(cover = if_else(Replicate == "ALL" & cover == 0.5, 1, cover))

cat("Present understory - plots:", n_distinct(prst_under_clean$plot_id), "\n")
cat("Replicate values found:\n")
print(prst_under_clean |> count(Replicate))

subplot_codes <- c("QN", "QE", "QS", "QW")

# Rows that are neither a subplot code nor "ALL" (e.g. "4m circle", blank)
# are not used in any output and are dropped.
prst_under_other <- prst_under_clean |> filter(!Replicate %in% c(subplot_codes, "ALL"))
if (nrow(prst_under_other) > 0) {
  cat("\nDropped", nrow(prst_under_other), "present rows with a non-subplot/non-ALL Replicate value:\n")
  print(prst_under_other |> count(Replicate))
}

prst_under_subplot_only <- prst_under_clean |> filter(Replicate %in% subplot_codes)


# ==============================================================================
# STEP 3a — Present-day occurrence (including "ALL") and cover (Q-subplots
# only, NA for ALL-only detections). This occurrence/cover table is used by
# every present-day output below (STEP 4 onward).
#
# Occurrence: 1 if detected in a Q-subplot or in "ALL" (the nested 50-m2
# presence/absence census; a species recorded only there is present in the
# plot).
# Cover: the Q-subplot mean over a full grid of the 4 Q-subplots (subplots
# where the species was absent count as 0), except when the species was
# detected only via "ALL", in which case cover is NA (no % cover was
# measured there; 0 would mean "measured and absent").
# ==============================================================================

species_per_plot_incl_all <- prst_under_clean |>
  filter(Replicate %in% c(subplot_codes, "ALL")) |>
  distinct(plot_id, species)

# Subplots are identified by census year as well, so that plots with two
# censuses contribute each census's subplots separately to the mean
subplots_sampled_per_plot <- prst_under_subplot_only |> distinct(plot_id, Year, Replicate)

q_true_zero_grid <- species_per_plot_incl_all |>
  inner_join(subplots_sampled_per_plot, by = "plot_id", relationship = "many-to-many") |>
  left_join(
    prst_under_subplot_only |> group_by(plot_id, Year, Replicate, species) |> summarise(cover = sum(cover, na.rm = TRUE), .groups = "drop"),
    by = c("plot_id", "Year", "Replicate", "species")
  ) |>
  mutate(cover = coalesce(cover, 0))

q_summary <- q_true_zero_grid |>
  group_by(plot_id, species) |>
  summarise(
    occurrence_Q = as.integer(any(cover > 0)),
    cover_Q_mean = mean(cover, na.rm = TRUE),
    .groups = "drop"
  )

all_detected <- prst_under_clean |>
  filter(Replicate == "ALL", !is.na(cover), cover > 0) |>
  distinct(plot_id, species) |>
  mutate(occurrence_ALL = 1L)

present_occ_cov <- species_per_plot_incl_all |>
  left_join(q_summary, by = c("plot_id", "species")) |>
  left_join(all_detected, by = c("plot_id", "species")) |>
  mutate(
    occurrence_ALL = coalesce(occurrence_ALL, 0L),
    occurrence = if_else(occurrence_Q == 1 | occurrence_ALL == 1, 1L, 0L),
    cover = if_else(occurrence_Q == 0 & occurrence_ALL == 1, NA_real_, cover_Q_mean),
    period = "present", plot_area = 500
  ) |>
  left_join(present_year, by = "plot_id") |>
  select(period, year, plot_id, plot_area, species, cover, occurrence, occurrence_Q, occurrence_ALL)

cat("\nPresent-day occurrence/cover incl. ALL:", nrow(present_occ_cov), "plot x species rows\n")
cat("  Detected via ALL only (cover forced to NA, not 0):",
    sum(present_occ_cov$occurrence_Q == 0 & present_occ_cov$occurrence_ALL == 1, na.rm = TRUE), "\n")

# occurrence_Q/cover_Q_mean are NA (not 0) for a plot with no Q-subplot
# rows (sampled only via "ALL"): q_true_zero_grid has no grid for that plot,
# so it has no q_summary row. This is the expected state for cover (no
# subplot measurement), and it is reported separately as a plot-level gap,
# distinct from the species-level ALL-only case above.
plots_with_no_subplot_sampling <- setdiff(species_per_plot_incl_all$plot_id, subplots_sampled_per_plot$plot_id)
cat("  Plots with ALL records but zero Q-subplot rows (cover NA for every",
    "species there):", length(plots_with_no_subplot_sampling), "\n")


# ==============================================================================
# STEP 4 — OUTPUT 1: understory_long.csv - main comparable table. Present
# occurrence/cover from STEP 3a ("ALL" included for occurrence, Q-subplots
# only for cover). Equivalent of overstory_long.csv.
# ==============================================================================

understory_long <- bind_rows(
  past_under_for_long |> select(period, year, plot_id, plot_area, species, cover, occurrence),
  present_occ_cov     |> select(period, year, plot_id, plot_area, species, cover, occurrence)
)

write_csv(understory_long, file.path(outputs_dir, "understory_long.csv"))
cat("\nSaved understory_long.csv:", nrow(understory_long), "rows\n")
cat("  Plots in past   :", n_distinct(understory_long$plot_id[understory_long$period == "past"]), "\n")
cat("  Plots in present:", n_distinct(understory_long$plot_id[understory_long$period == "present"]), "\n")


# ==============================================================================
# STEP 5 — OUTPUT 2: understory_long_subplot.csv - present-day detail per
# replicate, one row per QN/QE/QS/QW (1 m2 each) and "ALL" (50 m2, the
# nested presence/absence unit). The past data are added as they are (no
# subplots, plot_area from plot_info.csv) so that both periods are in one
# file.
# ==============================================================================

# area per replicate - Q-codes are 1 m2 subplots, "ALL" is the 50 m2 nested unit
replicate_area <- function(replicate) {
  if_else(replicate == "ALL", 50, 1)
}

prst_under_subplot_detail <- prst_under_clean |>
  filter(Replicate %in% c(subplot_codes, "ALL")) |>
  left_join(present_year, by = "plot_id") |>
  transmute(period = "present", year, plot_id,
            plot_area = replicate_area(Replicate),
            replicate = Replicate, species, cover,
            occurrence = if_else(!is.na(cover) & cover > 0, 1, 0))

past_under_subplot_detail <- past_under_for_long |>
  transmute(period = "past", year = 1975, plot_id, plot_area,
            replicate = NA_character_,
            species, cover, occurrence)

understory_long_subplot <- bind_rows(past_under_subplot_detail, prst_under_subplot_detail)

write_csv(understory_long_subplot, file.path(outputs_dir, "understory_long_subplot.csv"))
cat("\nSaved understory_long_subplot.csv:", nrow(understory_long_subplot), "rows\n")
cat("  Present replicates included:\n")
print(understory_long_subplot |> filter(period == "present") |> count(replicate, plot_area))


# ==============================================================================
# STEP 6 — OUTPUT 3: understory_presence.csv - presence/absence only;
# present day includes "ALL" together with the 4 subplots. This table has
# detections only (occurrence == 1 rows), so it is not a full 0/1
# occurrence matrix (see understory_long.csv /
# understory_occurrence_incl_all.csv for that).
# ==============================================================================

prst_under_presence <- prst_under_clean |>
  filter(Replicate %in% c(subplot_codes, "ALL")) |>
  filter(!is.na(cover), cover > 0) |>
  distinct(plot_id, species) |>
  left_join(present_year, by = "plot_id") |>
  mutate(period = "present", occurrence = 1)

past_under_presence <- past_under_long_raw |>
  filter(!is.na(cover), cover > 0) |>
  distinct(plot_id, species) |>
  mutate(period = "past", year = 1975, occurrence = 1)

understory_presence <- bind_rows(past_under_presence, prst_under_presence)

write_csv(understory_presence, file.path(outputs_dir, "understory_presence.csv"))
cat("\nSaved understory_presence.csv:", nrow(understory_presence), "rows\n")


# ==============================================================================
# STEP 7 — Restrict understory_long to plots/species with occurrence in
# both periods (same logic as overstory Step 7b), excluding plots without
# present-day cover subplots. Occurrence includes "ALL" (STEP 3a), so a
# species or plot recorded in both periods thanks to an ALL-only
# present-day detection is retained.
# ==============================================================================

# Saved for 02_explore_data: number of plots with understory occurrence
# only in the past, only in the present, or in both, before the
# both-periods restriction below.
understory_plot_coverage <- understory_long |>
  filter(occurrence == 1) |>
  distinct(plot_id, period) |>
  mutate(has = TRUE) |>
  pivot_wider(names_from = period, values_from = has, values_fill = FALSE) |>
  mutate(coverage = case_when(
    past & present ~ "both",
    past            ~ "only_past",
    present         ~ "only_present"
  )) |>
  select(plot_id, coverage)
write_csv(understory_plot_coverage, file.path(outputs_dir, "understory_plot_coverage.csv"))
cat("\nUnderstory plot coverage by period (before both-periods filter):\n")
print(understory_plot_coverage |> count(coverage))

plots_with_both_periods <- understory_long |>
  filter(occurrence == 1) |>
  distinct(plot_id, period) |>
  count(plot_id) |>
  filter(n == 2) |>
  pull(plot_id)

species_with_any_occurrence <- understory_long |>
  filter(occurrence == 1) |>
  distinct(species) |>
  pull(species)

n_plots_before   <- n_distinct(understory_long$plot_id)
n_species_before <- n_distinct(understory_long$species)

# Plots without present-day cover subplots (presence recorded only in the
# 50-m2 "ALL" census; STEP 3a) are excluded from the analysis tables, so
# that every retained plot has comparable occurrence and cover data in both
# surveys.
cat("\nPlots excluded for lacking present-day cover subplots:",
    paste(plots_with_no_subplot_sampling, collapse = ", "), "\n")
plots_with_both_periods <- setdiff(plots_with_both_periods, plots_with_no_subplot_sampling)

understory_long <- understory_long |>
  filter(plot_id %in% plots_with_both_periods, species %in% species_with_any_occurrence)

cat("\nScoped to plots/species with real occurrence in both periods:\n")
cat("  Plots  :", n_plots_before, "->", n_distinct(understory_long$plot_id), "\n")
cat("  Species:", n_species_before, "->", n_distinct(understory_long$species), "\n")

write_csv(understory_long, file.path(outputs_dir, "understory_long.csv"))


# ==============================================================================
# STEP 8 — Wide-format species x plot matrices (cover + occurrence), from
# the Step 7 understory_long. A species with cover == NA (ALL-only
# present-day detection) keeps that NA through pivot_wider; only absent
# plot x species combinations are filled with 0.
# ==============================================================================

understory_wide_cov <- understory_long |>
  select(period, year, plot_id, species, cover) |>
  pivot_wider(names_from = species, values_from = cover, values_fill = 0)

understory_wide_occ <- understory_long |>
  select(period, year, plot_id, species, occurrence) |>
  pivot_wider(names_from = species, values_from = occurrence, values_fill = 0)

write_csv(understory_wide_cov, file.path(outputs_dir, "understory_wide_cov.csv"))
write_csv(understory_wide_occ, file.path(outputs_dir, "understory_wide_occ.csv"))

cat("\nSaved wide-format matrices:\n")
cat("  understory_wide_cov:", nrow(understory_wide_cov), "rows x", ncol(understory_wide_cov) - 3, "species\n")
cat("  understory_wide_occ:", nrow(understory_wide_occ), "rows x", ncol(understory_wide_occ) - 3, "species\n")


# ==============================================================================
# STEP 9 — HMSC outputs: occurrence and cover given presence, not
# restricted to plots with both periods (STEP 7), so that HMSC uses the
# full present/past data rather than the thermophilization subset.
# Uses present_occ_cov from STEP 3a (same occurrence/cover definition as
# understory_long.csv; the only difference from STEP 4-8 is the absence of
# the both-periods restriction).
# ==============================================================================

past_occ_cov <- past_under_for_long |>
  mutate(occurrence_Q = NA_integer_, occurrence_ALL = NA_integer_) |>
  select(period, year, plot_id, species, occurrence, cover, occurrence_Q, occurrence_ALL)

present_occ_cov_full <- present_occ_cov |>
  select(period, year, plot_id, species, occurrence, cover, occurrence_Q, occurrence_ALL)

understory_occ_cov_incl_all <- bind_rows(past_occ_cov, present_occ_cov_full)

write_csv(understory_occ_cov_incl_all |> select(period, year, plot_id, species, occurrence),
          file.path(final_dir, "understory_occurrence_incl_all.csv"))
write_csv(understory_occ_cov_incl_all |> select(period, year, plot_id, species, cover),
          file.path(final_dir, "understory_cover_incl_all.csv"))

cat("\nSaved outputs/final/understory_occurrence_incl_all.csv:", nrow(understory_occ_cov_incl_all), "rows\n")
cat("Saved outputs/final/understory_cover_incl_all.csv:", nrow(understory_occ_cov_incl_all), "rows,",
    sum(is.na(understory_occ_cov_incl_all$cover)), "with NA cover (present-day, ALL-only detections)\n")

# Wide matrices (species x plot), the HMSC Y-input shape. Missing plot x
# species combinations are filled with 0 for both occurrence and cover
# (the species does not occur in that plot), as distinct from the NA cover
# of an ALL-only detection in a plot where the species does occur.
understory_wide_occ_incl_all <- understory_occ_cov_incl_all |>
  select(period, year, plot_id, species, occurrence) |>
  pivot_wider(names_from = species, values_from = occurrence, values_fill = 0)

understory_wide_cov_incl_all <- understory_occ_cov_incl_all |>
  select(period, year, plot_id, species, cover) |>
  pivot_wider(names_from = species, values_from = cover, values_fill = 0)

write_csv(understory_wide_occ_incl_all, file.path(final_dir, "understory_wide_occ_incl_all.csv"))
write_csv(understory_wide_cov_incl_all, file.path(final_dir, "understory_wide_cov_incl_all.csv"))

cat("\nSaved wide incl.-ALL matrices (outputs/final/):\n")
cat("  understory_wide_occ_incl_all:", nrow(understory_wide_occ_incl_all), "rows x", ncol(understory_wide_occ_incl_all) - 3, "species\n")
cat("  understory_wide_cov_incl_all:", nrow(understory_wide_cov_incl_all), "rows x", ncol(understory_wide_cov_incl_all) - 3, "species\n")