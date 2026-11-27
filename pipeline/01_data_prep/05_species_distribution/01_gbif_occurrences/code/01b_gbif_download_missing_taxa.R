# ------------------------------------------------------------------------------
# 01b_gbif_download_missing_taxa.R
# Module: pipeline/01_data_prep/05_species_distribution/01_gbif_occurrences
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Supplementary GBIF download for project species that 01_gbif_download_
# occurrences.R could not resolve to a species-rank taxon in the GBIF
# backbone (gbif_taxon_key_unmatched.csv):
#   - species whose project name is not in the backbone, which GBIF matched
#     to the genus only (e.g. Epilobium angustifolium -> genus Chamaenerion):
#     they are queried under the name GBIF accepts (alt_names below);
#   - varieties (e.g. Petasites frigidus var. palmatus): matched exactly but
#     excluded by the species-rank filter of 01; they are queried with the
#     variety's own taxon key, so that only records of that variety are used.
# Genus-level project names (e.g. Amelanchier) are not queried.
#
# The download uses the same filters as 01 (coordinates present, no
# geospatial issue, year <= 1985, North America). The species column of the
# downloaded records is replaced by the project species name, so that the
# records join the main download under the vegetation name in
# 02_gbif_clean_occurrences.R.
#
# Inputs : data_raw/gbif_taxon_key_unmatched.csv (01)
#          01_cleaning/outputs/spp_lib_current.csv
# Outputs (data_raw/):
#   gbif_taxon_key_matches_supplement.csv   project name -> query name -> key
#   gbif_occurrences_raw_<SCOPE>_supplement.csv
#   gbif_download_doi_<SCOPE>_supplement.txt
# Requires GBIF_USER, GBIF_PWD and GBIF_EMAIL in .Renviron (as for 01).
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(rgbif)

SCOPE <- "north_america"

module_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution",
                                 "01_gbif_occurrences", ...)

# Names accepted by the GBIF backbone for project species matched only to genus
alt_names <- tribble(
  ~project_species,          ~query_name,
  "Epilobium angustifolium", "Chamaenerion angustifolium",
  "Epilobium latifolium",    "Chamaenerion latifolium",
  "Blechnum spicant",        "Struthiopteris spicant",
  "Anemonoides grayi",       "Anemone grayi"
)

# ------------------------------------------------------------------------------
# 1. Taxa to query
# ------------------------------------------------------------------------------
unmatched <- read_csv(module_dir("data_raw", "gbif_taxon_key_unmatched.csv"), show_col_types = FALSE)

# Varieties / subspecies matched exactly: their own taxon key
# (speciesKey of the parent species is kept to assign records filed under a
# synonym of the variety, which carry the parent species key)
infraspecific <- unmatched %>%
  filter(rank %in% c("VARIETY", "SUBSPECIES"), matchType %in% c("EXACT", "FUZZY")) %>%
  transmute(project_species = verbatim_name, query_name = verbatim_name,
            usageKey, parent_species_key = as.numeric(speciesKey), rank, matchType, scientificName)

# Genus-only matches with a known alternative name: resolve the alternative
alt_resolved <- name_backbone_checklist(alt_names$query_name)
# acceptedUsageKey is only returned when at least one name is a synonym
if (!"acceptedUsageKey" %in% names(alt_resolved)) alt_resolved$acceptedUsageKey <- NA_real_
alt_resolved <- alt_resolved %>%
  select(verbatim_name, usageKey, rank, matchType, scientificName, status, acceptedUsageKey) %>%
  rename(query_name = verbatim_name) %>%
  right_join(alt_names, by = "query_name")
print(alt_resolved, width = Inf)

alt_ok <- alt_resolved %>%
  filter(rank %in% c("SPECIES", "VARIETY", "SUBSPECIES"), matchType %in% c("EXACT", "FUZZY")) %>%
  # synonyms in the backbone: query the accepted taxon
  mutate(usageKey = coalesce(as.numeric(acceptedUsageKey), as.numeric(usageKey))) %>%
  select(project_species, query_name, usageKey, rank, matchType, scientificName)

not_resolved <- setdiff(alt_names$project_species, alt_ok$project_species)
if (length(not_resolved) > 0) {
  warning("Alternative names not resolved in the GBIF backbone (not downloaded): ",
          paste(not_resolved, collapse = ", "))
}

supplement_keys <- bind_rows(infraspecific, alt_ok) %>% distinct(project_species, .keep_all = TRUE)
write_csv(supplement_keys, module_dir("data_raw", "gbif_taxon_key_matches_supplement.csv"))
cat("\nTaxa in the supplementary download:\n")
print(supplement_keys, width = Inf)

skipped <- setdiff(unmatched$verbatim_name, supplement_keys$project_species)
cat("\nNot queried (genus-level names or unresolved):", paste(skipped, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# 2. Download (same filters as 01)
# ------------------------------------------------------------------------------
# A completed download for these taxa is reused if present
doi_file <- module_dir("data_raw", paste0("gbif_download_doi_", SCOPE, "_supplement.txt"))
download_key <- NULL
if (file.exists(doi_file)) {
  key_line <- grep("^Download key:", readLines(doi_file), value = TRUE)
  if (length(key_line) == 1) {
    previous_key <- str_trim(str_remove(key_line, "^Download key:"))
    if (file.exists(module_dir("data_raw", paste0(previous_key, ".zip")))) download_key <- previous_key
  }
}
if (is.null(download_key)) {
  download_key <- occ_download(
    pred_in("taxonKey", unique(supplement_keys$usageKey)),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred_lte("year", 1985),
    pred("continent", "NORTH_AMERICA"),
    format = "SIMPLE_CSV"
  )
  occ_download_wait(download_key)
} else {
  cat("Reusing download", download_key, "\n")
}
cat("Download key:", download_key, "\n")

meta <- occ_download_meta(download_key)
writeLines(c(paste("GBIF occurrence download DOI:", paste0("https://doi.org/", meta$doi)),
             paste("Download key:", download_key),
             paste("Scope:", SCOPE),
             paste("Query date:", Sys.Date()),
             "Filter: year <= 1985, continent = NORTH_AMERICA, hasCoordinate = TRUE, hasGeospatialIssue = FALSE",
             paste("Taxa:", paste(supplement_keys$query_name, collapse = "; "))),
           doi_file)

zip_path <- occ_download_get(download_key, path = module_dir("data_raw"), overwrite = TRUE)
occ_supp <- occ_download_import(zip_path)

# ------------------------------------------------------------------------------
# 3. Assign records to project species
# ------------------------------------------------------------------------------
# A record belongs to a queried taxon when its taxonKey, speciesKey or
# acceptedTaxonKey equals the queried key (records of a variety carry the
# variety key in taxonKey; records of a species filed under a synonym carry
# the species key in speciesKey). Records filed under a synonym of a variety
# carry only the parent species key; since the download contains only the
# queried taxa, these are assigned to the variety through that key.
taxon_keys <- bind_rows(
  supplement_keys %>% transmute(project_species, key = as.numeric(usageKey)),
  supplement_keys %>% filter(!is.na(parent_species_key)) %>%
    transmute(project_species, key = parent_species_key)
) %>% distinct()
key_cols <- intersect(c("taxonKey", "acceptedTaxonKey", "speciesKey"), names(occ_supp))
n_downloaded <- nrow(occ_supp)
occ_supp <- occ_supp %>%
  mutate(.row = row_number()) %>%
  pivot_longer(all_of(key_cols), names_to = "key_col", values_to = "key", values_drop_na = TRUE) %>%
  mutate(key = as.numeric(key)) %>%
  inner_join(taxon_keys, by = "key", relationship = "many-to-many") %>%
  distinct(.row, project_species) %>%
  left_join(occ_supp %>% mutate(.row = row_number()), by = ".row") %>%
  mutate(species = project_species) %>%
  select(-.row, -project_species)

spp_lib <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", "spp_lib_current.csv"),
                    show_col_types = FALSE) %>%
  distinct(final_species, canopy_broad, life_form)
occ_supp <- occ_supp %>% left_join(spp_lib, by = c("species" = "final_species"))

cat("\nRecords downloaded:", n_downloaded, "| assigned to a project species:", nrow(occ_supp), "\n")
cat("Records per project species:\n")
print(count(occ_supp, species))

write_csv(occ_supp, module_dir("data_raw", paste0("gbif_occurrences_raw_", SCOPE, "_supplement.csv")))
