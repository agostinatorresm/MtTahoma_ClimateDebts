# ------------------------------------------------------------------------------
# 01_gbif_download_occurrences.R
# Module: pipeline/01_data_prep/05_species_distribution/01_gbif_occurrences
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: download GBIF occurrence records for every species in the MORA
# species library (over + understory combined), restricted to records at or
# before 1985. This predates the recent-warming period and keeps the derived
# species thermal niches (computed downstream by extracting CHELSA-bioclim
# 1981-2010 values at these occurrence points) from being biased by species
# having already shifted their realized distributions in response to more
# recent climate change.
#
# SCOPE: kept as a parameter for provenance, but the study uses
# SCOPE <- "north_america" only -- an initial unrestricted ("global")
# download was checked with a diagnostic continent breakdown and found that
# 83% of all records, driven by the most common species, came from Europe,
# where citizen-science reporting density is much higher for the same
# Eurasian-native species that are also common in North America. Left
# unrestricted, this would make each species' climatic affinity (mean
# annual temperature/precipitation across its occurrences) reflect European
# climate rather than the conditions under which the species actually grows
# in North America, which is what is relevant for comparison against the
# MORA plots. A convex-hull range-size trait (computed downstream) also only
# makes ecological sense within a single continent, reinforcing the
# North-America-only decision.
#   - SCOPE <- "north_america": adds pred("continent", "NORTH_AMERICA")
#     (GBIF's continent field, spanning the US, Canada, Mexico, Central
#     America and the Caribbean).
#   - SCOPE <- "global": no geographic restriction (kept in the code for
#     reference; not used in the final analysis).
#
# METHOD: uses rgbif::occ_download() (not occ_search()) because it is the
# only GBIF access method that returns a citable DOI for the exact query,
# which is needed for the methods section of a submission. occ_download()
# requires a GBIF account, with credentials (GBIF_USER, GBIF_PWD, GBIF_EMAIL)
# stored as environment variables (usually in .Renviron) -- this script stops
# early with instructions if they are not set, rather than failing deep
# inside the download call.
#
# Species names are resolved to GBIF backbone taxon keys first (via
# name_backbone_checklist()), restricted to exact/fuzzy matches at species
# rank, and any name that fails to resolve cleanly is written out separately
# for manual review rather than silently dropped or matched to the wrong
# taxon (e.g. falling back to genus level, which would pull in occurrences
# for congeners). This taxon resolution is scope-independent, so it is
# cached (data_raw/gbif_taxon_key_matches.csv) and reused across both runs
# instead of re-querying the GBIF backbone twice.
#
# Output (all filenames tagged with SCOPE, e.g. "..._global.csv" /
# "..._north_america.csv"):
#   - data_raw/gbif_taxon_key_matches.csv: species name -> resolved GBIF
#     taxonKey, match type and confidence (scope-independent, not tagged)
#   - data_raw/gbif_taxon_key_unmatched.csv: species that did not resolve
#     cleanly (scope-independent, not tagged)
#   - data_raw/gbif_download_doi_<SCOPE>.txt: the citable DOI for this exact
#     query
#   - data_raw/gbif_occurrences_raw_<SCOPE>.csv: the downloaded occurrence
#     records (species, coordinates, year, and other core Darwin Core
#     fields), joined back to canopy_broad (overstory/understory) from
#     spp_lib
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
# install.packages("rgbif") if not already installed
library(rgbif)

# 0a. Scope for this run -----------------------------------------------------------
# Set to "global" for one run and "north_america" for the other -- run the
# whole script once per value.
SCOPE <- "north_america"  # "global" or "north_america"
stopifnot(SCOPE %in% c("global", "north_america"))

cat("SCOPE =", SCOPE, "\n")

# 0b. Check GBIF credentials are available before doing anything else ----------
required_env <- c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL")
missing_env <- required_env[Sys.getenv(required_env) == ""]
if (length(missing_env) > 0) {
  stop(
    "Missing GBIF environment variables: ", paste(missing_env, collapse = ", "),
    ". Add them to .Renviron (usethis::edit_r_environ()) as:\n",
    "  GBIF_USER=<user>\n  GBIF_PWD=<password>\n  GBIF_EMAIL=<email>\n",
    "and restart R before running this script again."
  )
}

dir.create(here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw"),
           recursive = TRUE, showWarnings = FALSE)

# 0c. Use an existing unscoped download, if present -------------------------------
# A download saved as gbif_occurrences_raw.csv (unrestricted scope) is used as
# the "global" download when the scoped file does not exist, instead of
# requesting the same data from GBIF again.
legacy_raw <- here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw", "gbif_occurrences_raw.csv")
global_raw <- here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw", "gbif_occurrences_raw_global.csv")
if (SCOPE == "global" && file.exists(legacy_raw) && !file.exists(global_raw)) {
  cat("Found gbif_occurrences_raw.csv (earlier global download, no scope tag) -- copied to",
      basename(global_raw), "instead of downloading again.\n")
  file.copy(legacy_raw, global_raw)
  legacy_doi <- here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw", "gbif_download_doi.txt")
  global_doi <- here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw", "gbif_download_doi_global.txt")
  if (file.exists(legacy_doi) && !file.exists(global_doi)) file.copy(legacy_doi, global_doi)
  cat("Done -- the rest of this script is not needed for SCOPE = \"global\".\n")
  cat("(To download again from GBIF, delete", basename(global_raw), "and run again.)\n")
} else {
  
  # 1. Load species list -----------------------------------------------------------
  spp_lib <- read_csv(
    here("pipeline", "01_data_prep", "01_cleaning", "outputs", "spp_lib_current.csv"),
    show_col_types = FALSE
  )
  
  cat("Species in spp_lib_current.csv:", nrow(spp_lib), "\n")
  cat("  overstory:", sum(spp_lib$canopy_broad == "overstory"), "\n")
  cat("  understory:", sum(spp_lib$canopy_broad == "understory"), "\n")
  
  species_names <- spp_lib %>% distinct(final_species) %>% pull(final_species)
  cat("Unique species names to resolve in GBIF:", length(species_names), "\n")
  
  # 2. Resolve species names to GBIF backbone taxon keys (scope-independent, cached) --
  matches_path   <- here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw", "gbif_taxon_key_matches.csv")
  unmatched_path <- here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw", "gbif_taxon_key_unmatched.csv")
  
  if (file.exists(matches_path)) {
    cat("\nUsing the saved name resolution in gbif_taxon_key_matches.csv (independent of SCOPE).\n")
    backbone_matches <- read_csv(matches_path, show_col_types = FALSE)
  } else {
    cat("\nResolving names against the GBIF backbone...\n")
    backbone_matches <- name_backbone_checklist(species_names)
    write_csv(backbone_matches, matches_path)
  }
  
  clean_matches <- backbone_matches %>%
    filter(matchType %in% c("EXACT", "FUZZY"), rank == "SPECIES")
  
  unmatched <- backbone_matches %>%
    filter(!(matchType %in% c("EXACT", "FUZZY") & rank == "SPECIES"))
  
  if (nrow(unmatched) > 0 && !file.exists(unmatched_path)) {
    write_csv(unmatched, unmatched_path)
    cat("\n", nrow(unmatched), "species did not resolve to an exact species-level match --\n")
    cat("see gbif_taxon_key_unmatched.csv to check the completeness of the download.\n")
  }
  
  taxon_keys <- unique(clean_matches$usageKey)
  cat("\nTaxon keys resolved:", length(taxon_keys), "of", length(species_names), "species\n")
  
  # 3. Submit the GBIF occurrence download ------------------------------------------
  cat("\nSubmitting the GBIF download request (occ_download), SCOPE =", SCOPE, "...\n")
  
  predicates <- list(
    pred_in("taxonKey", taxon_keys),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred_lte("year", 1985)
  )
  if (SCOPE == "north_america") {
    predicates <- c(predicates, list(pred("continent", "NORTH_AMERICA")))
  }
  
  download_key <- do.call(occ_download, c(predicates, format = "SIMPLE_CSV"))
  
  cat("Download submitted, key:", download_key, "\n")
  cat("Waiting for GBIF to prepare the file (may take several minutes)...\n")
  
  occ_download_wait(download_key)
  
  # 4. Save the DOI for citation in the manuscript -----------------------------------
  download_meta <- occ_download_meta(download_key)
  doi_line <- paste0("https://doi.org/", download_meta$doi)
  
  filter_desc <- if (SCOPE == "north_america") {
    "year <= 1985, continent = NORTH_AMERICA, hasCoordinate = TRUE, hasGeospatialIssue = FALSE"
  } else {
    "year <= 1985, hasCoordinate = TRUE, hasGeospatialIssue = FALSE (no geographic restriction)"
  }
  
  writeLines(
    c(
      paste("GBIF occurrence download DOI:", doi_line),
      paste("Download key:", download_key),
      paste("Scope:", SCOPE),
      paste("Query date:", Sys.Date()),
      paste("Filter:", filter_desc),
      paste("Taxon keys:", length(taxon_keys), "species (from spp_lib_current.csv)")
    ),
    here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw", paste0("gbif_download_doi_", SCOPE, ".txt"))
  )
  
  cat("\nDownload DOI (", SCOPE, ") saved:", doi_line, "\n")
  
  # 5. Download and import the occurrence data ---------------------------------------
  cat("\nDownloading the occurrence file...\n")
  zip_path <- occ_download_get(
    download_key,
    path = here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw"),
    overwrite = TRUE
  )
  
  occurrences_raw <- occ_download_import(zip_path)
  
  cat("Occurrences downloaded (", SCOPE, "):", nrow(occurrences_raw), "\n")
  
  # 6. Join back to canopy_broad (overstory/understory) and save ---------------------
  # spp_lib can have more than one row per final_species (e.g. two prst_code
  # entries that both resolved to the same valid name), which would otherwise
  # turn this into a many-to-many join and silently duplicate occurrence rows
  # for the affected species. Deduplicating on (final_species, canopy_broad,
  # life_form) first collapses those cases as long as canopy_broad/life_form
  # agree across the duplicate rows -- if they don't (a genuine inconsistency
  # in spp_lib_current.csv), this stops so it gets resolved by hand instead of
  # silently picking one.
  spp_lib_dedup_check <- spp_lib %>%
    distinct(final_species, canopy_broad, life_form) %>%
    count(final_species) %>%
    filter(n > 1)
  
  if (nrow(spp_lib_dedup_check) > 0) {
    stop(
      "These species have inconsistent canopy_broad/life_form across rows of spp_lib_current.csv: ",
      paste(spp_lib_dedup_check$final_species, collapse = ", "),
      ". To be resolved before continuing."
    )
  }
  
  spp_lib_join <- spp_lib %>% distinct(final_species, canopy_broad, life_form)
  
  occurrences <- occurrences_raw %>%
    left_join(spp_lib_join, by = c("species" = "final_species"))
  
  stopifnot(nrow(occurrences) == nrow(occurrences_raw))
  
  write_csv(occurrences,
            here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences",
                 "data_raw", paste0("gbif_occurrences_raw_", SCOPE, ".csv")))
  
  cat("\n========================================\n")
  cat("GBIF DOWNLOAD COMPLETE (SCOPE =", SCOPE, ")\n")
  cat("========================================\n")
  cat("Occurrences saved in pipeline/01_data_prep/05_species_distribution/01_gbif_occurrences/data_raw/gbif_occurrences_raw_",
      SCOPE, ".csv\n", sep = "")
  cat("DOI to cite in the manuscript:", doi_line, "\n")
  if (SCOPE == "global") {
    cat("Next step: run this script with SCOPE <- \"north_america\".\n")
  } else {
    cat("Next step: 02_gbif_clean_occurrences.R (with SCOPE = \"north_america\").\n")
  }
  
}  # end of else block (actual download path)
