# -------------------------------------------------------------------------
# Author: Agostina Torres
# Project: MORAclimshift
# Script: 03_clean_spp_names.R
# Description:
# Standardizes the species list of the MORA dataset and builds the lookup
# tables used to join species codes to the vegetation data. A single lookup
# table cannot serve both joins, because prst_code (current 6-letter code)
# and past_code (1975 4-letter code) have different cardinality:
#
#   - present-day data (2014+) join by prst_code -> one row per prst_code
#   - 1975 data join by past_code -> a past_code maps to only one
#     prst_code, but a prst_code can have several past_codes (a species
#     was recorded under different 4-letter codes across surveys)
#
# About 51 species have several past_codes, so two separate lookup tables
# are built, both keyed to the scientific name (final_species) rather than
# the code.
#
# Source data:
# - MORA_PlantList_AllSurvey_Org.csv (original survey species list)
#   Species codes from 2023 (SpeciesCode2023) are the current codes.
#   Four-letter codes from 1975 (1975 4 Letter Code) are the main legacy
#   survey codes. 45 species have no value there but have one in the
#   alternate "1975/2014 4 letter code" column, which is used as a fallback
#   (e.g. POTR -> Populus trichocarpa).
#
# Outputs (pipeline/01_data_prep/01_cleaning/outputs/):
#   spp_lib_current.csv  : one row per prst_code, to join present-day
#                           tree/veg data. Includes the canopy_broad
#                           classification.
#   spp_lib_legacy.csv    : one row per past_code, to join legacy
#                            tree/veg data
#
# Exclusions:
# Two legacy codes (MITEL, PYEL) correspond to several current species that
# cannot be distinguished in the 1975 data. The species associated with
# these codes are excluded from both outputs:
#   MITEL -> Brewerimitella breweri, Ozomelis trifida, Pectiantia pentandra
#   PYEL  -> Pyrola elliptica, Pyrola minor
# -------------------------------------------------------------------------

library(tidyverse)
library(here)
library(TNRS)

org_dir     <- here("pipeline", "01_data_prep", "01_cleaning", "data", "data_raw", "01_raw_surveys")
outputs_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# STEP 1 — Load and format raw species list
# ==============================================================================

sppKV <- read_csv(file.path(org_dir, "MORA_PlantList_AllSurvey_Org.csv"), show_col_types = FALSE)

spp_raw <- sppKV |>
  select(
    prst_code = SpeciesCode2023,
    prst_name = SpeciesNames2023,
    past_code_primary  = `1975 4 Letter Code`,
    past_code_fallback = `1975/2014 4 letter code`,
    past_name = `1975 name`,
    id_level  = `Id Level`,
    family    = Family,
    life_form = PlantType,
    tree_type = TreeType
  ) |>
  # "not in list" is a placeholder meaning "no code", not a real code
  mutate(past_code_fallback = na_if(past_code_fallback, "not in list")) |>
  # Fallback: 45 species have no value in the primary 1975 column but
  # have one in the alternate "1975/2014" column (e.g. POTR -> Populus
  # trichocarpa); it is used when the primary is empty.
  mutate(past_code = coalesce(past_code_primary, past_code_fallback)) |>
  select(-past_code_primary, -past_code_fallback) |>
  # Correction of a legacy code typo
  mutate(past_code = if_else(past_code == "TSH4", "TSHE", past_code)) |>
  # "species" written in full is not a specific epithet but a placeholder
  # for an unidentified species of the genus, like "sp." or "spp.". These
  # names are reclassified here, before Step 2, so that they are handled
  # consistently by the species-level filter and enter the genus-level
  # catalog of Step 4 instead of being dropped from both.
  mutate(id_level = if_else(str_detect(prst_name, "(?i)\\bspecies$"), "Genus", id_level))

cat("Raw rows loaded            :", nrow(sppKV), "\n")
cat("Unique prst_codes at start :", n_distinct(spp_raw$prst_code), "\n")


# ==============================================================================
# STEP 2 — Clean species names and filter to species-level identifications
# ==============================================================================

spp_species_level <- spp_raw |>
  distinct(prst_code, .keep_all = TRUE) |>
  select(-past_name, -past_code) |>
  mutate(clean_name = str_squish(prst_name)) |>
  mutate(clean_name = str_remove_all(clean_name, "\\s?-\\s?type$")) |>
  mutate(id_level = case_when(
    str_detect(prst_name,  "\\b(spp?|spp-|sp)\\.?$")   ~ "Genus",
    str_detect(clean_name, "\\b(spp?|spp-|sp.|sp)\\.?$") ~ "Genus",
    TRUE ~ id_level
  )) |>
  filter(id_level == "Species") |>
  mutate(clean_name = str_to_sentence(str_to_lower(clean_name)))

cat("Unique prst_codes after species-level filter :", n_distinct(spp_species_level$prst_code), "\n")


# ==============================================================================
# STEP 3 — Standardize names via TNRS
# ==============================================================================

tnrs_res <- TNRS(spp_species_level$clean_name)

cat("Names resolved by TNRS :",
    sum(!is.na(tnrs_res$Accepted_name) & tnrs_res$Accepted_name != ""), "\n")
cat("Names kept as original :",
    sum( is.na(tnrs_res$Accepted_name) | tnrs_res$Accepted_name == ""), "\n")

spp_catalog_current <- spp_species_level |>
  left_join(
    select(tnrs_res, Name_submitted, Accepted_name, Accepted_family),
    by = c("clean_name" = "Name_submitted")
  ) |>
  mutate(
    final_species = if_else(!is.na(Accepted_name) & Accepted_name != "", Accepted_name, clean_name),
    final_family  = if_else(!is.na(Accepted_family) & Accepted_family != "", Accepted_family, family)
  ) |>
  select(prst_code, final_species, final_family, life_form, tree_type)

# Manual overrides from the review of TNRS discrepancies:
#   Group A - COLHET/FESIDA/HEUGLA were identified only to genus in the
#     field ("Collomia spp"/"Festuca spp"/"Heuchera spp"). A TNRS match on
#     the full name would assign them to a species, so they are kept at
#     genus level (id_level overridden) and excluded from spp_lib_current
#     like any other genus-only identification.
#   Group B - Aruncus dioicus var. acuminatus / Alnus alnobetula subsp.
#     sinuata: the parent-species name is used.
#   Group C - Mitella caulescens: the TNRS synonym is adopted
#     (Mitellastra caulescens -> Mitella caulescens), as for the
#     Listera -> Neottia names.
# Groups D/E (Spiraea splendens, Boykinia richardsonii) are not changed:
# the TNRS matches carry low-confidence warnings ("better higher taxonomic
# match available" / "ambiguous match"), so the original
# names are kept.
manual_species_overrides <- tribble(
  ~prst_code, ~override_final_species,
  "ARUDIO",   "Aruncus dioicus",
  "ALNVIR",   "Alnus alnobetula",
  "MITCAU",   "Mitella caulescens"
)
genus_only_overrides <- c("COLHET", "FESIDA", "HEUGLA")

spp_catalog_current <- spp_catalog_current |>
  left_join(manual_species_overrides, by = "prst_code") |>
  mutate(final_species = coalesce(override_final_species, final_species)) |>
  select(-override_final_species)

cat("Unique prst_codes after TNRS   :", n_distinct(spp_catalog_current$prst_code), "\n")
cat("Unique accepted species names  :", n_distinct(spp_catalog_current$final_species), "\n")

# ==============================================================================
# STEP 4 — Add missing present-day species (the 1975 ambiguous codes are
# resolved in Step 6, where past_code_lookup is built)
# ==============================================================================
# Checked against the project's species reference table (columns
# ProposedGroup/Code6/Code4/...), which is more complete than
# MORA_PlantList_AllSurvey_Org.csv:
#
#   - COCA, PYEL, MADI: not ambiguous. MORA_PlantList has an extra or wrong
#     destination for each that the reference table does not support:
#     COCA -> Cornus canadensis should be Cornus unalaschkensis; PYEL has
#     two destinations (Pyrola elliptica + Pyrola minor) but the reference
#     table supports only Pyrola minor; MADI has two destinations (Madia
#     gracilis + Maianthemum dilatatum) but Madia gracilis has no current
#     code in the reference table (never recorded in the present-day
#     survey), so MADI resolves to Maianthemum dilatatum. The three are
#     corrected in Step 6, where past_code_lookup is built.
#   - MITEL, VIOL: ambiguous in the 1975 data (3 species for MITEL, 2 for
#     VIOL share the same 1975 code in the reference table), but the
#     present-day data distinguish them with their own codes: VIOADU/VIOORB
#     are in MORA_PlantList_AllSurvey_Org.csv, and MITBRE/MITPEN/MITTRI,
#     missing from that list but used in the raw Veg Resurvey data (49
#     occurrences combined), are added below. The ambiguity is therefore
#     only in the past: the species are kept in the catalog and only the
#     1975 "MITEL"/"VIOL" codes are removed from past_code_lookup (Step 6).

# The 3 missing present-day species (codes used in the raw Veg Resurvey
# data but absent from MORA_PlantList_AllSurvey_Org.csv) are added and run
# through TNRS like the other species.
missing_present_species <- tribble(
  ~prst_code, ~raw_name,
  "MITBRE",   "Mitella breweri",
  "MITPEN",   "Mitella pentandra",
  "MITTRI",   "Mitella trifida"
)

missing_tnrs <- TNRS(missing_present_species$raw_name)

missing_present_species <- missing_present_species |>
  left_join(select(missing_tnrs, Name_submitted, Accepted_name, Accepted_family),
            by = c("raw_name" = "Name_submitted")) |>
  transmute(
    prst_code,
    final_species = if_else(!is.na(Accepted_name) & Accepted_name != "", Accepted_name, raw_name),
    final_family  = Accepted_family,
    life_form = "Forb", tree_type = NA_character_  # matches the other Mitella-group entries
  )

cat("\nAdded missing present-day species (not in MORA_PlantList_AllSurvey_Org.csv):\n")
print(missing_present_species |> select(prst_code, final_species))

# Genus-level and other entries (removed by the species-level filter of
# Step 2) are also part of the full catalog, tagged as such. They are
# built here from spp_raw without TNRS, since they are already "Genus spp."
# names rather than binomials.
spp_catalog_genus <- spp_raw |>
  distinct(prst_code, .keep_all = TRUE) |>
  filter(id_level != "Species") |>
  transmute(prst_code, final_species = prst_name, final_family = family,
            life_form, tree_type, catalog_status = "genus_level")

spp_catalog_current <- spp_catalog_current |>
  mutate(catalog_status = if_else(prst_code %in% genus_only_overrides,
                                  "genus_level", "species_level")) |>
  bind_rows(spp_catalog_genus) |>
  bind_rows(missing_present_species |> mutate(catalog_status = "species_level"))


write_csv(spp_catalog_current, file.path(outputs_dir, "spp_catalog_current.csv"))
cat("\nSaved spp_catalog_current.csv:", nrow(spp_catalog_current), "rows (FULL catalog, nothing dropped)\n")
cat("Catalog status breakdown:\n")
print(spp_catalog_current |> count(catalog_status))

spp_with_final_name <- spp_catalog_current |>
  filter(catalog_status == "species_level") |>
  select(-catalog_status)


# ==============================================================================
# STEP 5 — spp_lib_current: one row per prst_code, for joining PRESENT data
# ==============================================================================

spp_lib_current <- spp_with_final_name |>
  distinct(prst_code, .keep_all = TRUE)

cat("\nspp_lib_current rows (= unique prst_codes) :", nrow(spp_lib_current), "\n")


# ==============================================================================
# STEP 6 — spp_lib_legacy: one row per past_code, for joining LEGACY data
# ==============================================================================
# A past_code maps to only one prst_code/final_species, but a prst_code can
# have several past_codes across surveys. Building this table separately,
# keyed by past_code, gives each join a table with the right cardinality.

past_code_lookup <- spp_raw |>
  select(prst_code, past_code) |>
  filter(!is.na(past_code), past_code != "") |>
  # Data-entry errors in MORA_PlantList_AllSurvey_Org.csv, checked against
  # the project's species reference table: each of these 3 codes has an
  # extra or wrong destination that the reference table does not support,
  # so each resolves to a single species:
  #   COCA -> CORCAN (Cornus canadensis) corrected to CORUNA (Cornus
  #           unalaschkensis), a field-identification correction; the two
  #           species remain distinct, neither excluded nor merged
  #   PYEL -> the reference table supports only PYRMIN (Pyrola minor), not PYRELI
  #   MADI -> Madia gracilis has no current code in the reference table
  #           (never recorded in the present-day survey), so the destination
  #           is MAIDIL (Maianthemum dilatatum)
  mutate(prst_code = case_when(
    past_code == "COCA" ~ "CORUNA",
    past_code == "PYEL" ~ "PYRMIN",
    past_code == "MADI" ~ "MAIDIL",
    TRUE ~ prst_code
  )) |>
  # MITEL/VIOL: ambiguous in 1975 (3 species share MITEL and 2 share VIOL in
  # the reference table), but the present-day data distinguish them with
  # their own codes (MITBRE/MITPEN/MITTRI added in Step 4; VIOADU/VIOORB
  # already in the catalog). The ambiguity is only in the past, so these two
  # past_codes are removed from the lookup rather than assigned to one
  # species; present-day data for these species are unaffected.
  filter(!past_code %in% c("MITEL", "VIOL"))

# Full past-side catalog (against the full current-side catalog, including
# genus-level entries), built from past_code_lookup, which resolves
# COCA/PYEL/MADI to a single destination and drops MITEL/VIOL (see above).
spp_catalog_legacy <- past_code_lookup |>
  inner_join(spp_catalog_current, by = "prst_code") |>
  distinct(past_code, final_species, .keep_all = TRUE) |>
  select(past_code, final_species, final_family, life_form, tree_type, catalog_status)

write_csv(spp_catalog_legacy, file.path(outputs_dir, "spp_catalog_legacy.csv"))
cat("\nSaved spp_catalog_legacy.csv:", nrow(spp_catalog_legacy), "rows (FULL catalog, nothing dropped)\n")

spp_lib_legacy <- past_code_lookup |>
  inner_join(spp_with_final_name, by = "prst_code") |>
  distinct(past_code, .keep_all = TRUE) |>
  select(past_code, final_species, final_family, life_form, tree_type)

cat("spp_lib_legacy rows (= unique past_codes)  :", nrow(spp_lib_legacy), "\n")

# Check: no past_code should map to more than one final_species, since
# spp_lib_legacy has one row per past_code.
ambiguous_check <- past_code_lookup |>
  inner_join(spp_with_final_name, by = "prst_code") |>
  distinct(past_code, final_species) |>
  count(past_code) |>
  filter(n > 1)

if (nrow(ambiguous_check) > 0) {
  cat("\nWARNING: the following past_codes still map to multiple final_species\n")
  cat("(should be empty after the MITEL/PYEL exclusion in Step 4):\n")
  print(ambiguous_check)
} else {
  cat("Confirmed: every past_code maps to exactly one final_species.\n")
}


# ==============================================================================
# STEP 7 — canopy_broad: overstory / understory (binary botanical
# classification: whether the species can be a canopy tree at MORA)
# ==============================================================================
# Whether a species has Tree data in both periods (tree_data_coverage) is
# computed from the Tree observations in 04_apply_spp_lib.R rather
# than from a fixed species list, so that the
# classification follows the data.

species_overstory_full <- c(
  "Abies amabilis", "Abies grandis", "Abies lasiocarpa", "Abies procera",
  "Acer circinatum", "Acer glabrum", "Acer macrophyllum",
  "Alnus alnobetula", "Alnus rubra", "Betula papyrifera",
  "Callitropsis nootkatensis", "Malus fusca", "Picea engelmannii",
  "Picea sitchensis", "Pinus albicaulis", "Pinus contorta", "Pinus monticola",
  "Pinus ponderosa", "Populus balsamifera", "Prunus emarginata",
  "Pseudotsuga menziesii", "Sorbus sitchensis", "Taxus brevifolia",
  "Thuja plicata", "Tsuga heterophylla", "Tsuga mertensiana"
)

classify_canopy_broad <- function(species) {
  if_else(species %in% species_overstory_full, "overstory", "understory")
}

# Attached to spp_lib_current (not a separate spp_list file) so that the
# canopy classification travels with the species catalog;
# 04_apply_spp_lib.R uses it to build spp_list_final.
spp_lib_current <- spp_lib_current |>
  mutate(canopy_broad = classify_canopy_broad(final_species))

cat("\nCanopy broad breakdown (unique species):\n")
print(spp_lib_current |> distinct(final_species, canopy_broad) |> count(canopy_broad))


# ==============================================================================
# STEP 8 — spp_info: species metadata, separate from the code-matching
# tables (spp_lib_current/legacy resolve codes; spp_info describes species
# attributes - family, life form, canopy - once resolved).
#
# life_form (from PlantType): Forb / Tree / Shrub / Grass / Fern / Sedge /
#   Rush - NA for ~66 species not categorized in the source list.
# tree_type (from TreeType): Conifer / Non-coniferous - only populated for
#   life_form == "Tree", NA otherwise.
# ==============================================================================

spp_info <- spp_lib_current |>
  distinct(species = final_species, family = final_family, life_form, tree_type,
           canopy_broad) |>
  arrange(species)

cat("\nspp_info rows (= unique species) :", nrow(spp_info), "\n")
cat("life_form breakdown:\n")
print(spp_info |> count(life_form))


# ==============================================================================
# STEP 9 — Save
# ==============================================================================

write_csv(spp_lib_current, file.path(outputs_dir, "spp_lib_current.csv"))
write_csv(spp_lib_legacy,  file.path(outputs_dir, "spp_lib_legacy.csv"))
write_csv(spp_info,        file.path(outputs_dir, "spp_info.csv"))

cat("\nSaved spp_lib_current.csv, spp_lib_legacy.csv, spp_info.csv\n")
