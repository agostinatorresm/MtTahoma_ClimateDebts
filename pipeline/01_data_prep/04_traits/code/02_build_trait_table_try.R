# ------------------------------------------------------------------------------
# 02_build_trait_table_try.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: extract the six traits available from TRY (SLA, seed dry mass,
# plant lifespan/longevity, wood density, growth form, woodiness) from the
# raw TRY exports, reduce each to one value per species, and match them to
# the project species list (overstory + understory combined) by direct
# match with a WFO synonym fallback.
#
# TRY only. GIFT is added next (03_build_trait_table_gift.R) and genus-level
# imputation after that (05_impute_trait_by_genus.R); the output of this
# script is an intermediate table, not the final per-layer trait table.
#
# Species with no TRY record under their accepted name get a synonym
# fallback: if TRY has the trait under one of the species' WFO synonyms
# (from 01_build_species_synonyms.R), that record is used. Each value
# carries a "<trait>_src" tag ("direct", "synonym:<name used>", or
# "unresolved"), so that the coverage gained through synonyms can be traced.
#
# Numeric traits (SLA, seed mass, longevity, wood density) are reduced to the
# mean of TRY's standardized StdValue across all records for that species.
# Categorical traits (growth form, woodiness) do not populate StdValue, so
# OrigValueStr is used instead. TRY pools many original datasets with
# inconsistent coding for these two traits, so each raw value is first
# harmonized to a small set of target categories (section 2b) before
# taking the per-species mode:
#   - growth_form -> tree / shrub / vine / grass (incl. sedges) / herb
#   - woodiness    -> woody / non-woody
# A raw value matching none or more than one category (unrecognized
# spelling, or a mixed report like "shrub/tree") is excluded.
#
# TRAIT SOURCE FILES (raw TRY exports, one species x trait record per row):
#   - 38172 (dispersal request):        TraitID 26   Seed dry mass
#   - 40000 (leaves & longevity request): TraitID 59  Plant lifespan (longevity)
#                                         TraitIDs 3115/3116/3117 SLA (leaf
#                                         area per leaf dry mass -- pooled,
#                                         same units, differ only in whether
#                                         the petiole is included)
#   - 52607 (wood density request):      TraitID 4    Stem specific density
#     (SSD, stem dry mass per unit stem fresh volume), stem tissue density,
#     wood density (TraitName field of this file).
#   - 52732 (growth form & woodiness request): TraitID 42 Plant growth form
#                                               TraitID 38 Plant woodiness
# (38363, phenology/flowering time, is not used here; it is not one of the
# HMSC traits.)
#
# SPECIES LIST: distinct species of overstory_long.csv + understory_long.csv
# with at least one occurrence (Q-subplot or ALL), restricted (both layers)
# to species_pool_broad.csv, the project species pool
# (pipeline/01_data_prep/02_species_pool; species present in both surveys,
# no plot-count threshold), used for the main trait and HMSC analyses. The
# pool is defined once upstream so that every module uses the same species
# list.
#
# VALUE USED: TRY's StdValue (already unit-standardized across original
# studies) for numeric traits, mode of OrigValueStr for categorical traits.
#
# Output:
#   - tmp_files/try_species_trait_summary.csv: one row per species, direct
#     TRY match only (no synonym fallback yet), columns SLA/seed_mass/
#     longevity/wood_density (mean StdValue) + growth_form/woodiness (mode of
#     the harmonized category) + n_records per trait
#   - tmp_files/try_species_trait_summary_resolved.csv: full project pool,
#     after direct match + WFO synonym fallback, with "<trait>_src" provenance
#   - tmp_files/try_trait_coverage_summary.csv: % of the project's species
#     with at least one TRY record, per trait, split overstory/understory
#   - tmp_files/try_species_missing_<trait>.csv: species with zero TRY
#     records (direct or synonym), per trait -- candidates for GIFT / genus
#     imputation downstream
#   - tmp_files/trait_matching_intermediate.rds: raw per-trait record tables +
#     trait_species_resolved + species pools + synonym table, loaded by the
#     next scripts instead of re-reading the raw TRY exports (several
#     hundred MB, ~2.5M rows for growth_form alone)
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(data.table)

# 0a. Module directories --------------------------------------------------------
this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)
cleaning_dir     <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)

# outputs/ holds only the final per-layer trait tables, written at the end
# of 05_impute_trait_by_genus.R; everything here -- raw per-species TRY
# reductions, the coverage summary, the per-trait missing-species lists -- is
# intermediate/diagnostic and goes to tmp_files/.
for (d in c("outputs", "tmp_files")) {
  dir.create(this_module_dir(d), recursive = TRUE, showWarnings = FALSE)
}

# 0b. Raw TRY export paths, relative to this module's data_raw ------------------
# Flat in data_raw/ (no per-request subfolders).
TRY_FILES <- c(
  dispersal        = this_module_dir("data_raw", "38172.txt"),
  leaves_longevity = this_module_dir("data_raw", "40000.txt"),
  wood_density     = this_module_dir("data_raw", "52607.txt"),
  lifeform         = this_module_dir("data_raw", "52732.txt")
)

# 1. Load species list: overstory + understory, with the species names of the
#    vegetation tables (the same names as the species pool and the WFO synonym
#    table), restricted to species with at least one real occurrence
#    (Q-subplot or ALL) anywhere in the long table. Names are not changed to
#    other taxonomic treatments (e.g. Berberis -> Mahonia) before matching:
#    the direct match uses the vegetation name and the synonym fallback covers
#    every WFO-listed alternative name, so TRY records filed under either
#    treatment are found ---------------------------------------------------------
species_pool_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)

overstory_long_raw  <- read_csv(cleaning_dir("outputs", "overstory_long.csv"),  show_col_types = FALSE)
understory_long_raw <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE)

# Project species pool (built in 02_species_pool from the same two long
# tables). species_pool_broad.csv (species present in both surveys, no
# plot-count threshold) is used for both layers; species_pool_strict.csv
# (additional >5 plots per period threshold) is used for HMSC and is not
# needed for trait matching.
species_pool_broad <- read_csv(species_pool_dir("data_clean", "species_pool_broad.csv"), show_col_types = FALSE)

# Occurrence filter: a species with occurrence == 0 in every row (never
# detected, in either period or sampling grain) is a catalog entry rather
# than a community member, and is dropped from the trait species pool for
# both layers.
drop_zero_occurrence_species <- function(long_df, stratum_label) {
  occ_by_species <- long_df %>%
    group_by(species) %>%
    summarise(total_occ = sum(occurrence, na.rm = TRUE), .groups = "drop")
  
  n_phantom <- occ_by_species %>% filter(total_occ == 0) %>% nrow()
  cat(stratum_label, ": dropped", n_phantom,
      "species with zero real occurrence anywhere (template artifacts).\n")
  
  occ_by_species %>% filter(total_occ > 0) %>% pull(species)
}

overstory_real_species  <- drop_zero_occurrence_species(overstory_long_raw,  "Overstory")
understory_real_species <- drop_zero_occurrence_species(understory_long_raw, "Understory")

# Both layers are restricted to species_pool_broad.csv -- the centralized
# species pool (see header note above), keyed on the species names of
# overstory_long.csv / understory_long.csv.
overstory_species <- tibble(species = overstory_real_species) %>%
  inner_join(species_pool_broad %>% filter(stratum == "overstory") %>% select(species),
             by = "species") %>%
  mutate(stratum = "overstory")

cat("Overstory: ", length(overstory_real_species), "species with real occurrence ->",
    nrow(overstory_species), "in species_pool_broad.csv.\n")

understory_species <- tibble(species = understory_real_species) %>%
  inner_join(species_pool_broad %>% filter(stratum == "understory") %>% select(species),
             by = "species") %>%
  distinct(species) %>%
  mutate(stratum = "understory")

cat("Understory: ", length(unique(understory_real_species)), "species with real occurrence ->",
    nrow(understory_species), "in species_pool_broad.csv.\n")

project_species <- bind_rows(overstory_species, understory_species) %>%
  distinct(species, .keep_all = TRUE)  # a species in both layers counts once for coverage

cat("Project species (overstory + understory, deduplicated, real occurrences only):",
    nrow(project_species), "\n")

# 1b. Species synonym table (WFO), restricted to the project species
#     pool; used below as a fallback when a species has no TRY record
#     under its accepted name but TRY has it under a synonym. Built by
#     01_build_species_synonyms.R from the WFO Dec-2025 COLDP release; it
#     covers the full species catalog (spp_lib_current.csv, including names
#     never recorded) and is filtered here to the species recorded in the
#     vegetation data.
spp_synonyms_all <- read_csv(this_module_dir("outputs", "spp_synonyms.csv"), show_col_types = FALSE)

# Alternative names not resolved through WFO's exact-name match: names under
# another genus treatment and orthographic variants of the epithet. Added to
# the synonym table so that TRY / GIFT records filed under them are found.
extra_synonyms <- tribble(
  ~final_species,                ~synonym,
  "Berberis aquifolium",         "Mahonia aquifolium",
  "Berberis nervosa",            "Mahonia nervosa",
  "Monotropa hypopitys",         "Hypopitys monotropa",
  "Lactuca muralis",             "Mycelis muralis",
  "Trautvetteria caroliniensis", "Trautvetteria carolinensis",
  "Mitella trifida",             "Ozomelis trifida"
)

spp_synonyms <- spp_synonyms_all %>%
  select(final_species, synonym) %>%
  bind_rows(extra_synonyms) %>%
  distinct(final_species, synonym) %>%
  semi_join(project_species, by = c("final_species" = "species"))

cat("\nSpecies-synonym fallback table:", n_distinct(spp_synonyms$final_species),
    "of", nrow(project_species), "project species have >=1 WFO-listed synonym",
    "available as a TRY matching fallback.\n")

# 2. Extract each trait from its raw TRY file ------------------------------------
read_try_trait <- function(path, trait_id) {
  if (!file.exists(path)) {
    cat("FILE NOT FOUND:", path, "-- skipping this trait.\n")
    return(tibble(species = character(), value = double()))
  }
  
  dt <- fread(path, sep = "\t", header = TRUE, quote = "", encoding = "Latin-1", showProgress = FALSE)
  dt_trait <- dt[TraitID %in% trait_id] %>% as_tibble()
  
  # Drop likely-erroneous records: TRY's own ErrorRisk column flags values that
  # deviate more than N standard deviations from the species/trait mean (data
  # entry or unit errors); the standard threshold used across TRY workflows is
  # to exclude ErrorRisk > 4. Records with no ErrorRisk (NA) are kept.
  n_before <- nrow(dt_trait)
  if ("ErrorRisk" %in% names(dt_trait)) {
    dt_trait <- dt_trait %>% filter(is.na(ErrorRisk) | ErrorRisk <= 4)
    cat("  ErrorRisk filter (>4 dropped):", n_before - nrow(dt_trait), "of", n_before, "records removed\n")
  } else {
    cat("  No ErrorRisk column found in this file -- skipped that filter.\n")
  }
  
  dt_trait %>%
    filter(!is.na(StdValue), !is.na(AccSpeciesName)) %>%
    transmute(species = AccSpeciesName, value = StdValue)
}

seed_mass_raw  <- read_try_trait(TRY_FILES["dispersal"],        26)   # TraitID 26: Seed dry mass (StdValue unit: mg)
longevity_raw  <- read_try_trait(TRY_FILES["leaves_longevity"], 59)   # TraitID 59: Plant lifespan
# SLA: pool all three TRY variants of the same measurement (same units,
# mm2/mg) that differ only in whether the petiole is included in the leaf
# area used.
sla_raw        <- read_try_trait(TRY_FILES["leaves_longevity"], c(3115, 3116, 3117))
wood_dens_raw  <- read_try_trait(TRY_FILES["wood_density"],     4)    # TraitID 4: Wood density (SSD)

cat("\nRaw records read:\n")
cat(" seed_mass:    ", nrow(seed_mass_raw), "\n")
cat(" longevity:    ", nrow(longevity_raw), "\n")
cat(" SLA:          ", nrow(sla_raw), "\n")
cat(" wood_density: ", nrow(wood_dens_raw), "\n")

# 2a. Categorical traits (growth form, woodiness) --------------------------------
# These don't populate StdValue, so they're extracted from OrigValueStr instead
# and reduced to a per-species mode (most frequent value) rather than a mean.
# Same ErrorRisk > 4 filter as the numeric traits.
read_try_categorical_trait <- function(path, trait_id) {
  if (!file.exists(path)) {
    cat("FILE NOT FOUND:", path, "-- skipping this trait.\n")
    return(tibble(species = character(), value = character()))
  }
  
  dt <- fread(path, sep = "\t", header = TRUE, quote = "", encoding = "Latin-1", showProgress = FALSE)
  dt_trait <- dt[TraitID == trait_id] %>% as_tibble()
  
  n_before <- nrow(dt_trait)
  if ("ErrorRisk" %in% names(dt_trait)) {
    dt_trait <- dt_trait %>% filter(is.na(ErrorRisk) | ErrorRisk <= 4)
    cat("  ErrorRisk filter (>4 dropped):", n_before - nrow(dt_trait), "of", n_before, "records removed\n")
  } else {
    cat("  No ErrorRisk column found in this file -- skipped that filter.\n")
  }
  
  dt_trait %>%
    filter(!is.na(OrigValueStr), OrigValueStr != "", !is.na(AccSpeciesName)) %>%
    transmute(species = AccSpeciesName, value = str_trim(OrigValueStr))
}

growth_form_raw_all <- read_try_categorical_trait(TRY_FILES["lifeform"], 42)  # TraitID 42: Plant growth form
woodiness_raw_all    <- read_try_categorical_trait(TRY_FILES["lifeform"], 38)  # TraitID 38: Plant woodiness

cat(" growth_form (raw, before harmonization):", nrow(growth_form_raw_all), "\n")
cat(" woodiness (raw, before harmonization):  ", nrow(woodiness_raw_all), "\n")

# 2b. Harmonize raw categories down to a small, unambiguous target set -----------
# TRY pools many original datasets, each with its own coding (full words,
# abbreviations, mixed-form strings like "shrub/tree", undocumented numeric
# codes). Rather than trust a mode taken over ~1700 raw growth-form spellings
# (or ~40 raw woodiness spellings) as-is, each raw value is classified into a
# fixed target category by keyword match; anything matching zero or more than
# one target category (an unrecognized spelling, or a genuinely mixed-form
# report like "shrub/tree" or "non-woody/woody") is left unmapped (NA) rather
# than guessed. This is deliberately conservative -- coverage below reflects
# only the records giving unambiguous evidence for one category.
#
# growth_form target categories: tree, shrub, vine, herb, grass
# (graminoid/sedge folded into grass). Single-letter codes are only resolved
# where the TRY convention is well established (the T/S/H growth-form coding
# seen across several datasets, e.g. LEDA): T = tree, S = shrub, H = herb.
harmonize_growth_form <- function(raw_value) {
  v <- str_to_lower(str_trim(raw_value))
  is_tree  <- str_detect(v, "tree") | v == "t"
  is_shrub <- str_detect(v, "shrub") | v == "s"
  is_vine  <- str_detect(v, "vine") | str_detect(v, "liana") | str_detect(v, "climb")
  is_grass <- str_detect(v, "grass") | str_detect(v, "graminoid") | str_detect(v, "sedge")
  is_herb  <- str_detect(v, "herb") | str_detect(v, "\\bforb\\b") | v == "h"
  
  n_matches <- is_tree + is_shrub + is_vine + is_grass + is_herb
  
  case_when(
    n_matches != 1 ~ NA_character_,  # unrecognized, or a mixed-form report -- excluded, not guessed
    is_tree  ~ "tree",
    is_shrub ~ "shrub",
    is_vine  ~ "vine",
    is_grass ~ "grass",
    is_herb  ~ "herb",
    TRUE     ~ NA_character_
  )
}

# woodiness target categories: woody, non-woody. Single-letter
# codes resolved only where well established: W = woody, H = herbaceous.
# A mixed-form report (e.g. "non-woody/woody") is excluded, not guessed.
# "Herb"/"herb" alone (distinct from the growth-form trait's own "herb") is
# also treated as non-woody -- under the woodiness trait specifically, it can
# only mean herbaceous.
harmonize_woodiness <- function(raw_value) {
  v <- str_to_lower(str_trim(raw_value))
  is_mixed    <- str_detect(v, "/") | str_detect(v, "\\band\\b")
  is_nonwoody <- str_detect(v, "non-?\\s*woody") | str_detect(v, "herbaceous") | str_detect(v, "\\bherb\\b") | v == "h"
  is_woody    <- (str_detect(v, "woody") & !is_nonwoody) | v == "w"
  
  case_when(
    is_mixed    ~ NA_character_,
    is_nonwoody ~ "non-woody",
    is_woody    ~ "woody",
    TRUE        ~ NA_character_
  )
}

growth_form_raw <- growth_form_raw_all %>%
  mutate(value = harmonize_growth_form(value)) %>%
  filter(!is.na(value))

woodiness_raw <- woodiness_raw_all %>%
  mutate(value = harmonize_woodiness(value)) %>%
  filter(!is.na(value))

cat("\n growth_form -- mapped:", nrow(growth_form_raw), "of", nrow(growth_form_raw_all),
    "raw records (", round(100 * nrow(growth_form_raw) / nrow(growth_form_raw_all), 1), "% )\n")
cat(" growth_form harmonized categories:\n")
print(count(growth_form_raw, value, sort = TRUE))

cat("\n woodiness -- mapped:", nrow(woodiness_raw), "of", nrow(woodiness_raw_all),
    "raw records (", round(100 * nrow(woodiness_raw) / nrow(woodiness_raw_all), 1), "% )\n")
cat(" woodiness harmonized categories:\n")
print(count(woodiness_raw, value, sort = TRUE))

cat("\n growth_form -- top unmapped raw values (check none of these should have been\n",
    " caught by the keyword match -- e.g. a language variant or synonym missed):\n")
print(growth_form_raw_all %>% filter(is.na(harmonize_growth_form(value))) %>% count(value, sort = TRUE) %>% slice_head(n = 15))

cat("\n woodiness -- top unmapped raw values:\n")
print(woodiness_raw_all %>% filter(is.na(harmonize_woodiness(value))) %>% count(value, sort = TRUE) %>% slice_head(n = 15))

# 3. Reduce each trait to one value per species ----------------------------------
# Numeric traits: mean of StdValue. Categorical traits: mode of OrigValueStr
# (most frequent category across records for that species; ties broken by
# first-encountered category).
summarise_trait <- function(df, trait_label) {
  df %>%
    group_by(species) %>%
    summarise(value = mean(value, na.rm = TRUE), n_records = n(), .groups = "drop") %>%
    rename_with(~ paste0(trait_label, "_", .x), c(value, n_records))
}

summarise_categorical_trait <- function(df, trait_label) {
  df %>%
    count(species, value, name = "n") %>%
    group_by(species) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(df %>% count(species, name = "n_records"), by = "species") %>%
    select(species, value, n_records) %>%
    rename_with(~ paste0(trait_label, "_", .x), c(value, n_records))
}

trait_species <- summarise_trait(seed_mass_raw, "seed_mass") %>%
  full_join(summarise_trait(longevity_raw, "longevity"), by = "species") %>%
  full_join(summarise_trait(sla_raw, "SLA"), by = "species") %>%
  full_join(summarise_trait(wood_dens_raw, "wood_density"), by = "species") %>%
  full_join(summarise_categorical_trait(growth_form_raw, "growth_form"), by = "species") %>%
  full_join(summarise_categorical_trait(woodiness_raw, "woodiness"), by = "species")

write_csv(trait_species, this_module_dir("tmp_files", "try_species_trait_summary.csv"))

# 3b. Resolve each trait for the project species, with a synonym fallback -------
# Canonical match first (species found directly under its own accepted name
# in trait_species); for anything still missing, check whether TRY has the
# trait recorded under one of that species' WFO-listed synonyms. When more
# than one synonym has TRY data for a trait, the one backed by the most
# records is used (same most-supported-wins logic as the direct per-species
# reduction above). Provenance per trait is kept in a "<trait>_src" column.
resolve_with_synonym_fallback <- function(trait_species, trait_label, project_species, spp_synonyms) {
  value_col <- paste0(trait_label, "_value")
  n_col     <- paste0(trait_label, "_n_records")
  src_col   <- paste0(trait_label, "_src")
  
  base <- project_species %>%
    select(species) %>%
    left_join(trait_species %>% select(species, all_of(value_col), all_of(n_col)), by = "species")
  
  direct <- base %>%
    filter(!is.na(.data[[value_col]])) %>%
    mutate(!!src_col := "direct")
  
  fallback <- base %>%
    filter(is.na(.data[[value_col]])) %>%
    select(species) %>%
    inner_join(spp_synonyms, by = c("species" = "final_species")) %>%
    inner_join(
      trait_species %>% select(synonym = species, all_of(value_col), all_of(n_col)),
      by = "synonym"
    ) %>%
    filter(!is.na(.data[[value_col]])) %>%
    group_by(species) %>%
    slice_max(.data[[n_col]], n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(!!src_col := paste0("synonym:", synonym)) %>%
    select(species, all_of(value_col), all_of(n_col), all_of(src_col))
  
  unresolved <- project_species %>%
    select(species) %>%
    anti_join(direct, by = "species") %>%
    anti_join(fallback, by = "species") %>%
    mutate(!!value_col := NA, !!n_col := NA_integer_, !!src_col := "unresolved")
  
  bind_rows(direct, fallback, unresolved) %>% arrange(species)
}

trait_labels <- c("seed_mass", "longevity", "SLA", "wood_density", "growth_form", "woodiness")

trait_species_resolved <- project_species %>% select(species)
for (tr in trait_labels) {
  resolved_tr <- resolve_with_synonym_fallback(trait_species, tr, project_species, spp_synonyms)
  trait_species_resolved <- trait_species_resolved %>% left_join(resolved_tr, by = "species")
}

write_csv(trait_species_resolved, this_module_dir("tmp_files", "try_species_trait_summary_resolved.csv"))

cat("\nSynonym-fallback gains, per trait (species with no direct TRY match that",
    "were resolved via a WFO synonym instead):\n")
for (tr in trait_labels) {
  src_col <- paste0(tr, "_src")
  n_fallback <- sum(str_starts(trait_species_resolved[[src_col]], "synonym:"), na.rm = TRUE)
  cat(" ", tr, ":", n_fallback, "species gained via synonym fallback\n")
}

# 4. Coverage against the project species list, per trait x layer -------------
# Uses trait_species_resolved (direct match + synonym fallback), not the raw
# trait_species table, so the coverage numbers below reflect the fallback.
# wood_density is a woody-species trait -- it only applies to overstory (trees);
# understory herbs/shrubs don't have a meaningful wood density, so that
# combination is dropped rather than reported as a coverage gap to fill.
# growth_form and woodiness apply to both layers (every species, woody or not,
# has a legitimate growth form and woodiness value), so no layer filter for
# those two.
coverage <- project_species %>%
  left_join(trait_species_resolved, by = "species") %>%
  mutate(across(c(seed_mass_value, longevity_value, SLA_value, wood_density_value), as.character)) %>%
  pivot_longer(
    cols = c(seed_mass_value, longevity_value, SLA_value, wood_density_value,
             growth_form_value, woodiness_value),
    names_to = "trait", values_to = "value"
  ) %>%
  mutate(trait = str_remove(trait, "_value$"), has_data = !is.na(value)) %>%
  filter(!(trait == "wood_density" & stratum == "understory")) %>%
  group_by(stratum, trait) %>%
  summarise(
    n_species = n(),
    n_with_data = sum(has_data),
    pct_covered = round(100 * n_with_data / n_species, 1),
    .groups = "drop"
  ) %>%
  arrange(trait, stratum)

write_csv(coverage, this_module_dir("tmp_files", "try_trait_coverage_summary.csv"))

cat("\nCoverage (% of species with at least one TRY record), by trait x stratum:\n")
print(coverage)

# 5. Species missing each trait -- the GIFT / genus-imputation candidates -------
# After direct match + synonym fallback (section 3b): species with no TRY
# record under any name.
for (tr in trait_labels) {
  missing_sp <- project_species %>%
    left_join(trait_species_resolved, by = "species") %>%
    filter(is.na(.data[[paste0(tr, "_value")]])) %>%
    select(species, stratum)
  
  write_csv(missing_sp, this_module_dir("tmp_files", paste0("try_species_missing_", tr, ".csv")))
  cat(tr, "-- missing for", nrow(missing_sp), "species (written to tmp_files/try_species_missing_", tr, ".csv)\n", sep = "")
}

# 6. Summary ----------------------------------------------------------------------
cat("\n========================================\n")
cat("TRY TRAIT x SPECIES COVERAGE (direct + synonym) -- COMPLETE\n")
cat("========================================\n")
cat("No final outputs written here -- this is an intermediate resolution step.\n")
cat("Next: 03_build_trait_table_gift.R folds in GIFT; 04_impute_trait_by_genus.R\n")
cat("adds genus-mean/mode imputation and writes the final outputs/overstory_traits.csv\n")
cat("and outputs/understory_traits.csv.\n")
cat("\nIntermediate/diagnostic (tmp_files/):\n")
cat(" - try_species_trait_summary.csv (direct TRY match only, no fallback)\n")
cat(" - try_species_trait_summary_resolved.csv (direct + synonym fallback, full project pool)\n")
cat(" - try_trait_coverage_summary.csv (% covered per trait x stratum)\n")
cat(" - try_species_missing_<trait>.csv (x6, genuinely unresolved after fallback)\n")
cat("\nNOTE: growth_form (tree/shrub/vine/grass/herb) and woodiness\n",
    "(woody/non-woody) were harmonized from TRY's raw coding by keyword match\n",
    "(section 2b) -- check the two 'top unmapped raw values' lists printed above\n",
    "for anything that should actually have been caught (a missed synonym or\n",
    "language variant); everything else unmapped was a genuinely ambiguous or\n",
    "unrecognized code and was correctly excluded rather than guessed.\n")

# 7. Save intermediate objects for downstream scripts ----------------------------
# (03_build_trait_table_gift.R, 05_impute_trait_by_genus.R,
# 06_build_phylogeny.R) -- avoids re-parsing the raw TRY exports
# (several hundred MB, ~2.5M rows for growth_form alone) every time a
# downstream script runs.
saveRDS(
  list(
    seed_mass_raw = seed_mass_raw, longevity_raw = longevity_raw,
    sla_raw = sla_raw, wood_dens_raw = wood_dens_raw,
    growth_form_raw = growth_form_raw, woodiness_raw = woodiness_raw,
    trait_species_resolved = trait_species_resolved,
    project_species = project_species,
    overstory_species = overstory_species,
    understory_species = understory_species,
    spp_synonyms = spp_synonyms
  ),
  this_module_dir("tmp_files", "trait_matching_intermediate.rds")
)
cat("\nSaved tmp_files/trait_matching_intermediate.rds for downstream scripts.\n")
