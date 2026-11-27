# ------------------------------------------------------------------------------
# 01_build_species_synonyms.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: build a species x synonym table for the project's full resolved
# species catalog (spp_lib_current.csv's final_species, the TNRS-accepted
# names already used everywhere downstream), using WFO (World Flora Online)
# as the synonym source. TNRS itself only resolves a submitted name to its
# accepted name -- it does not expose the list of OTHER names (synonyms)
# that also point to that accepted name. WFO is one of TNRS's own backbone
# sources, so matching the already-TNRS-accepted names against WFO's
# accepted taxa and pulling every synonym pointing back at that same taxon
# gives a synonym list consistent with the taxonomic treatment already
# adopted in 03_clean_spp_names.R, rather than introducing a second,
# potentially conflicting backbone.
#
# This runs FIRST in the trait/phylogeny pipeline (before any TRY or GIFT
# matching) because both of those matching cascades (direct match, then WFO
# synonym fallback) depend on this table already existing.
#
# Data source: WFO Plant List, Catalogue of Life Data Package (COLDP) format,
# 2025-12 release (files: taxon.tsv, name.tsv, synonym.tsv). This differs
# from the "classic" WFO backbone expected by WorldFlora::WFO.match() (which
# needs classification.csv) -- COLDP exposes the accepted/synonym
# relationship directly via taxonID/nameID keys, so no fuzzy-matching
# package is needed for this step; matching is done on exact scientificName
# strings (scientificName excludes authorship in this format, consistent
# with the binomial-only format of final_species).
#
# COLDP table roles used here:
#   - name.tsv:    one row per scientific name string (ID, scientificName).
#                  Both accepted names and synonym names live in this table.
#   - taxon.tsv:   one row per ACCEPTED taxon (ID = taxonID, nameID = the
#                  name.tsv ID of its accepted name).
#   - synonym.tsv: one row per synonym (taxonID = the accepted taxon it
#                  belongs to, nameID = the name.tsv ID of the synonym's
#                  OWN name string).
#
# Kept in 04_traits (not 01_cleaning) -- its only consumer right now is
# the trait-matching cascade (TRY and GIFT), so it lives next to that rather
# than in the species-identity module; move it later if a second consumer
# shows up.
#
# Output:
#   - outputs/spp_synonyms.csv: one row per (final_species, synonym) pair --
#     final_species is the project's TNRS-accepted name (the join key used
#     everywhere else), synonym is one WFO-listed alternate name for it.
#     A species with zero WFO-listed synonyms does not appear in this table
#     at all (nothing to add over the canonical name already in spp_lib).
#   - outputs/spp_synonyms_match_summary.csv: one row per final_species --
#     whether/how it matched a WFO taxon, and how many synonyms (if any)
#     WFO lists for it. Species with zero synonyms gain nothing from the
#     synonym-fallback cascade regardless of how well everything else works.
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)
cleaning_dir     <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)

# Path to the extracted WFO release.
wfo_dir <- this_module_dir("data_raw", "WFO_Backbone", "WFO_Backbone")
dir.create(this_module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)

# --- Load the three COLDP tables, essential columns only ----------------------
wfo_name <- read_tsv(
  file.path(wfo_dir, "name.tsv"),
  col_select = c(ID, scientificName, authorship, rank, genus, specificEpithet),
  show_col_types = FALSE
)

wfo_taxon <- read_tsv(
  file.path(wfo_dir, "taxon.tsv"),
  col_select = c(ID, nameID),
  show_col_types = FALSE
) %>%
  rename(taxonID = ID, accepted_nameID = nameID)

wfo_synonym <- read_tsv(
  file.path(wfo_dir, "synonym.tsv"),
  col_select = c(taxonID, nameID),
  show_col_types = FALSE
) %>%
  rename(synonym_nameID = nameID)

# --- Project species to resolve ------------------------------------------------
spp_lib_current <- read_csv(cleaning_dir("outputs", "spp_lib_current.csv"), show_col_types = FALSE)

project_species <- spp_lib_current %>%
  distinct(final_species) %>%
  filter(!is.na(final_species), final_species != "") %>%
  mutate(final_species = str_trim(final_species))

cat("Project species to match against WFO:", nrow(project_species), "\n")

# --- Step 1: resolve each project species to a WFO accepted taxonID -----------
# Path A: the project species name matches a name used as ACCEPTED in taxon.tsv
match_as_accepted <- project_species %>%
  inner_join(wfo_name, by = c("final_species" = "scientificName")) %>%
  inner_join(wfo_taxon, by = c("ID" = "accepted_nameID")) %>%
  distinct(final_species, taxonID) %>%
  mutate(match_type = "accepted")

# Path B: for species not resolved above, check whether the name is itself
# filed as a SYNONYM in WFO -- if so, follow it to its accepted taxonID. This
# covers cases where TNRS and WFO disagree on which name is currently accepted.
unresolved_1 <- project_species %>% anti_join(match_as_accepted, by = "final_species")

match_as_synonym <- unresolved_1 %>%
  inner_join(wfo_name, by = c("final_species" = "scientificName")) %>%
  inner_join(wfo_synonym, by = c("ID" = "synonym_nameID")) %>%
  distinct(final_species, taxonID) %>%
  mutate(match_type = "resolved_via_synonym")

resolved <- bind_rows(match_as_accepted, match_as_synonym)

# Homonym guard: a name string matching more than one WFO taxon is ambiguous --
# excluded rather than guessed, flagged for manual review.
ambiguous_species <- resolved %>% count(final_species) %>% filter(n > 1) %>% pull(final_species)
if (length(ambiguous_species) > 0) {
  cat("WARNING:", length(ambiguous_species), "species matched more than one WFO taxon",
      "(homonym/ambiguous) -- excluded from synonym resolution:\n")
  print(ambiguous_species)
}
resolved <- resolved %>% filter(!final_species %in% ambiguous_species)

unresolved_final <- project_species %>% anti_join(resolved, by = "final_species")

cat(nrow(resolved), "of", nrow(project_species), "project species resolved to a WFO taxon;",
    nrow(unresolved_final), "not found in WFO under this exact name.\n")

# --- Step 2: pull every synonym listed for each resolved taxonID --------------
spp_synonyms <- resolved %>%
  select(final_species, taxonID) %>%
  # many-to-many is expected here: one taxonID can have several synonym
  # rows, and all of them are kept.
  inner_join(wfo_synonym, by = "taxonID", relationship = "many-to-many") %>%
  left_join(wfo_name, by = c("synonym_nameID" = "ID")) %>%
  filter(!is.na(scientificName), scientificName != final_species) %>%
  distinct(final_species, synonym = scientificName) %>%
  arrange(final_species, synonym)

# --- Match / coverage summary --------------------------------------------------
spp_synonyms_match_summary <- project_species %>%
  left_join(resolved %>% select(final_species, match_type), by = "final_species") %>%
  mutate(match_type = replace_na(match_type, "unresolved")) %>%
  left_join(
    spp_synonyms %>% count(final_species, name = "n_synonyms"),
    by = "final_species"
  ) %>%
  mutate(n_synonyms = replace_na(n_synonyms, 0L))

cat("\nMatch type breakdown:\n")
print(count(spp_synonyms_match_summary, match_type))
cat("\nSpecies gaining >=1 synonym from this table:",
    sum(spp_synonyms_match_summary$n_synonyms > 0), "of", nrow(project_species), "\n")

# --- Save -----------------------------------------------------------------------
write_csv(spp_synonyms, this_module_dir("outputs", "spp_synonyms.csv"))
write_csv(spp_synonyms_match_summary, this_module_dir("outputs", "spp_synonyms_match_summary.csv"))
