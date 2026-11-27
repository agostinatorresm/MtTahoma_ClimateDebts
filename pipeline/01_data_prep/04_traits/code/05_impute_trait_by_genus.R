# ------------------------------------------------------------------------------
# 05_impute_trait_by_genus.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: genus-level imputation for species x trait cells still unresolved
# after the species-level sources (TRY 02, GIFT 03, BIEN 03b).
#
# Genus used for a species: the genus of its vegetation name; only when that
# genus has no records for the trait, the genera of its synonyms (e.g. Mahonia
# for Berberis nervosa). Synonym genera are a fallback and are not pooled with
# the accepted genus, because WFO synonyms include old combinations in other
# genera (e.g. Abies menziesii for Picea sitchensis and Pseudotsuga
# menziesii), which would mix unrelated genera. The genus value is computed
# from all records of the genus in each source's full download (any species,
# not only project species), so it is available even when no project congener
# has data.
#
# Evidence pooled per trait (the same pool is used in the LOO validation, 09):
#   seed_mass, SLA : TRY + GIFT + BIEN raw records
#   growth_form    : TRY + GIFT
#   longevity, wood_density, woodiness : TRY
# Numeric traits: mean of all records of the genus. Categorical traits: mode
# of all records of the genus.
#
# Provenance: values resolved upstream keep their <trait>_src tag. Values
# imputed here are tagged "genus_mean_<sources>:<genera>" or
# "genus_mode_<sources>:<genus>" (e.g. "genus_mean_try+gift+bien:Tiarella"),
# so the evidence behind each value is visible in the table. A cell
# whose genus set has no records stays "unresolved".
#
# Literature values for single species are added afterwards in
# 08_integration.R, only to cells still unresolved here.
#
# Inputs : tmp_files/trait_matching_intermediate.rds (02, 03, 03b)
# Outputs: tmp_files/try_species_trait_summary_genus_imputed.csv  all project
#            species, after species-level sources and genus imputation
#          tmp_files/try_species_still_missing_<trait>.csv  species still
#            unresolved per trait
#          outputs/overstory_traits.csv, outputs/understory_traits.csv
#            per-layer trait tables (finalised by 08_integration.R)
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

obj <- readRDS(this_module_dir("tmp_files", "trait_matching_intermediate.rds"))
list2env(obj, envir = environment())

trait_labels       <- c("seed_mass", "longevity", "SLA", "wood_density", "growth_form", "woodiness")
categorical_traits <- c("growth_form", "woodiness")

# Starting point: species-level values after TRY + GIFT + BIEN
# (trait_species_resolved_bien, 03b); TRY + GIFT only if 03b has not been run.
start_table <- if (exists("trait_species_resolved_bien")) {
  trait_species_resolved_bien
} else {
  cat("NOTE: trait_species_resolved_bien not found; run 03b_build_trait_table_bien.R.",
      "Proceeding from TRY + GIFT.\n")
  trait_species_resolved_try_gift
}

raw_by_trait <- list(
  seed_mass = seed_mass_raw, longevity = longevity_raw, SLA = sla_raw,
  wood_density = wood_dens_raw, growth_form = growth_form_raw, woodiness = woodiness_raw
)
pool_sources <- setNames(as.list(rep("try", length(trait_labels))), trait_labels)

add_source <- function(label) {
  if (!exists(label, inherits = TRUE)) {
    cat("NOTE:", label, "not found in the .rds; genus summaries proceed without it.\n")
    return(invisible(NULL))
  }
  extra <- get(label)
  for (tr in intersect(names(extra), trait_labels)) {
    raw_by_trait[[tr]] <<- bind_rows(raw_by_trait[[tr]], extra[[tr]] %>% select(species, value))
    pool_sources[[tr]] <<- c(pool_sources[[tr]], sub("_raw_by_trait", "", label))
  }
}
add_source("gift_raw_by_trait")
add_source("bien_raw_by_trait")
cat("Genus evidence pooled per trait:\n")
for (tr in trait_labels) cat(" ", tr, ":", paste(pool_sources[[tr]], collapse = " + "), "\n")

extract_genus <- function(species_binomial) word(species_binomial, 1)

# Candidate genera of each project species: genus of the vegetation name
# (priority 1) and genera of its synonyms (priority 2, fallback)
species_genera <- bind_rows(
  project_species %>% transmute(species, genus = extract_genus(species), priority = 1L),
  spp_synonyms %>% filter(!is.na(synonym)) %>%
    transmute(species = final_species, genus = extract_genus(synonym), priority = 2L)
) %>%
  semi_join(project_species, by = "species") %>%
  group_by(species, genus) %>% summarise(priority = min(priority), .groups = "drop")

# Records available to each project species: those of its own genus if any,
# otherwise those of its synonyms' genera
genus_records <- function(raw_df) {
  species_genera %>%
    inner_join(raw_df %>% mutate(genus = extract_genus(species)) %>% select(genus, rec_species = species, value),
               by = "genus", relationship = "many-to-many") %>%
    group_by(species) %>%
    slice_min(priority, with_ties = TRUE) %>%
    ungroup()
}

genus_summary_numeric <- function(raw_df) {
  genus_records(raw_df) %>%
    group_by(species) %>%
    summarise(genus_value = mean(value, na.rm = TRUE),
              genus_n_records = n(),
              genus_n_species = n_distinct(rec_species),
              genera = paste(sort(unique(genus)), collapse = "/"),
              .groups = "drop")
}

genus_summary_categorical <- function(raw_df) {
  pooled <- genus_records(raw_df)
  mode_tab <- pooled %>%
    count(species, value, name = "n") %>%
    group_by(species) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(species, genus_value = value)
  support_tab <- pooled %>% group_by(species) %>%
    summarise(genus_n_records = n(), genera = paste(sort(unique(genus)), collapse = "/"), .groups = "drop")
  mode_tab %>% left_join(support_tab, by = "species") %>% mutate(genus_n_species = NA_integer_)
}

# For one trait: fill "unresolved" rows from the genus table, keep everything
# already resolved (direct/synonym/GIFT) untouched.
impute_genus <- function(resolved_tbl, trait_label, raw_df, is_categorical) {
  value_col <- paste0(trait_label, "_value")
  n_col     <- paste0(trait_label, "_n_records")
  src_col   <- paste0(trait_label, "_src")
  method    <- paste0(if (is_categorical) "genus_mode_" else "genus_mean_",
                      paste(pool_sources[[trait_label]], collapse = "+"))
  
  genus_tab <- if (is_categorical) genus_summary_categorical(raw_df) else genus_summary_numeric(raw_df)
  
  fill <- resolved_tbl %>%
    select(species, !!src_col) %>%
    filter(.data[[src_col]] == "unresolved") %>%
    left_join(genus_tab, by = "species") %>%
    filter(!is.na(genus_value)) %>%
    select(species, genus_value, genus_n_records, genera)
  
  cat(" ", trait_label, ": genus imputation resolves", nrow(fill), "additional species\n")
  
  out <- resolved_tbl %>%
    left_join(fill, by = "species")
  if (!is_categorical) out <- out %>% mutate(!!value_col := as.numeric(.data[[value_col]]))
  out %>%
    mutate(
      !!src_col   := if_else(!is.na(genus_value), paste0(method, ":", genera), .data[[src_col]]),
      !!n_col     := if_else(!is.na(genus_value), as.integer(genus_n_records), as.integer(.data[[n_col]])),
      !!value_col := coalesce(.data[[value_col]], genus_value)
    ) %>%
    select(-genus_value, -genus_n_records, -genera)
}

cat("Genus-level imputation gains, per trait:\n")
trait_species_genus_imputed <- start_table
for (tr in trait_labels) {
  trait_species_genus_imputed <- impute_genus(
    trait_species_genus_imputed, tr, raw_by_trait[[tr]], tr %in% categorical_traits
  )
}

write_csv(trait_species_genus_imputed, this_module_dir("tmp_files", "try_species_trait_summary_genus_imputed.csv"))

# Species still unresolved after genus imputation --------------------------------
cat("\nStill unresolved after TRY + GIFT + BIEN + genus imputation:\n")
for (tr in trait_labels) {
  src_col <- paste0(tr, "_src")
  still_missing_sp <- trait_species_genus_imputed %>%
    filter(.data[[src_col]] == "unresolved") %>%
    select(species)
  write_csv(still_missing_sp, this_module_dir("tmp_files", paste0("try_species_still_missing_", tr, ".csv")))
  cat(" ", tr, ": ", nrow(still_missing_sp), " species\n", sep = "")
}

# Final per-layer trait tables -------------------------------------------------
# wood_density is dropped from the understory table (not a structurally
# meaningful trait for non-woody/sub-canopy species); every other trait
# applies to both layers.
clean_final_table <- function(species_list, resolved_table, drop_wood_density) {
  out <- species_list %>%
    select(species) %>%
    left_join(resolved_table, by = "species") %>%
    rename_with(~ str_remove(.x, "_value$"), ends_with("_value"))
  if (drop_wood_density) out <- out %>% select(-wood_density, -wood_density_n_records, -wood_density_src)
  out %>% arrange(species)
}

overstory_traits  <- clean_final_table(overstory_species,  trait_species_genus_imputed, FALSE)
understory_traits <- clean_final_table(understory_species, trait_species_genus_imputed, TRUE)

write_csv(overstory_traits,  this_module_dir("outputs", "overstory_traits.csv"))
write_csv(understory_traits, this_module_dir("outputs", "understory_traits.csv"))

cat("\n========================================\n")
cat("TRAIT TABLES (TRY + GIFT + BIEN + genus imputation)\n")
cat("========================================\n")
cat("outputs/overstory_traits.csv: ", nrow(overstory_traits), "species x all 6 traits\n")
cat("outputs/understory_traits.csv:", nrow(understory_traits), "species x 5 traits (no wood_density)\n")

# Number of distinct species missing at least one applicable trait (the
# per-trait counts above can count a species more than once).
count_missing_any <- function(traits_df) {
  src_cols <- names(traits_df) %>% str_subset("_src$")
  traits_df %>%
    mutate(any_missing = rowSums(across(all_of(src_cols), ~ .x == "unresolved")) > 0)
}

overstory_missing_any  <- count_missing_any(overstory_traits)
understory_missing_any <- count_missing_any(understory_traits)

cat("\nDistinct species missing at least one applicable trait",
    "(after species-level sources and genus imputation):\n")
cat(" overstory: ",  sum(overstory_missing_any$any_missing),  "of", nrow(overstory_traits),  "\n")
cat(" understory:", sum(understory_missing_any$any_missing), "of", nrow(understory_traits), "\n")


# Traits missing for each species (species x trait detail), e.g. "Carex sp.
# lacks longevity and SLA", rather than only the number of species with
# some missing trait.
missing_trait_detail <- function(traits_df) {
  src_cols <- names(traits_df) %>% str_subset("_src$")
  traits_df %>%
    select(species, all_of(src_cols)) %>%
    rename_with(~ str_remove(.x, "_src$"), all_of(src_cols)) %>%
    pivot_longer(-species, names_to = "trait", values_to = "src") %>%
    filter(src == "unresolved") %>%
    group_by(species) %>%
    summarise(missing_traits = paste(trait, collapse = ", "), n_missing = n(), .groups = "drop") %>%
    arrange(desc(n_missing), species)
}

overstory_missing_detail  <- missing_trait_detail(overstory_traits)
understory_missing_detail <- missing_trait_detail(understory_traits)

cat("\nOverstory -- species x which trait(s) still missing:\n")
print(overstory_missing_detail, n = Inf)

cat("\nUnderstory -- species x which trait(s) still missing:\n")
print(understory_missing_detail, n = Inf)

# missing_trait_detail (species + exactly which trait(s), written below) is
# strictly more informative than a bare species list, so that's the only
# missing-species file written from here -- no separate species-only file.
write_csv(overstory_missing_detail,  this_module_dir("tmp_files", "overstory_missing_trait_detail.csv"))
write_csv(understory_missing_detail, this_module_dir("tmp_files", "understory_missing_trait_detail.csv"))

# Number of understory species with all three traits of a chosen set
# (SLA, seed_mass, woodiness); missing longevity/wood_density/growth_form
# does not count here, since those traits are not in the set.
check_traits <- c("SLA", "seed_mass", "woodiness")
n_complete_set <- understory_traits %>%
  filter(if_all(all_of(paste0(check_traits, "_src")), ~ .x != "unresolved")) %>%
  nrow()

cat("\nUnderstory species with", paste(check_traits, collapse = " + "), "all resolved:",
    n_complete_set, "of", nrow(understory_traits), "\n")
