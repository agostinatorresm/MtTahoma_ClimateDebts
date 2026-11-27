# ------------------------------------------------------------------------------
# 01_combine_gbif_chelsa.R
# Module: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: combine the cleaned GBIF occurrences (North America, year <= 1985;
# gbif_occurrences module) with CHELSA-bioclim 1981-2010 (bio1 = mean annual
# temperature, bio12 = annual precipitation; downloaded by
# 00_download_chelsa_bioclim.R in this same module) to compute
# each species' climatic affinity: a single truncated-mean value for use as
# an HMSC trait, and a set of 100 bootstrap resamples for propagating
# uncertainty into the community-level thermophilization/CTI analysis
# (Fadrique et al., 2018 approach).
#
# SCOPE: North America only (see gbif_occurrences/02_gbif_clean_occurrences.R).
# A convex hull spanning several continents has no ecological meaning as a
# species' range, so the geographic scope is kept continental throughout.
#
# SYNONYM HARMONIZATION: GBIF occurrences are filed under the name assigned
# by the GBIF taxonomic backbone, which does not always match the project's
# species name (e.g. occurrences of Callitropsis nootkatensis are filed
# under "Xanthocyparis nootkatensis"). Occurrences are therefore relabelled
# to the project species name (via the WFO synonym table of
# 04_traits/02_build_trait_table_try.R) before computing any per-species
# affinity, so that occurrences filed under a synonym count toward the
# project species in every downstream script. Requires the spp_synonyms
# table saved by 04_traits/02_build_trait_table_try.R (section 1b).
#
# METHOD:
#   1. Extract bio1 and bio12 at every cleaned occurrence point (local
#      GeoTIFFs from data_raw/chelsa_bioclim/, downloaded by
#      00_download_chelsa_bioclim.R).
#   2. For each species, exclude values outside the 5th-95th percentile of
#      its own distribution (reduces influence of outlying/erroneous
#      records).
#   3. Single-value affinity (for HMSC traits): mean of the truncated
#      distribution.
#   4. Uncertainty propagation (for thermophilization/CTI): 100 values
#      resampled with replacement from the truncated distribution, kept in
#      long format so downstream scripts can join on (species, rep).
#
# TEMPERATURE UNITS: checked empirically (if raw bio1 > 200, assumed Kelvin
# and 273.15 is subtracted) rather than assumed, since CHELSA's monthly and
# bioclim products are not guaranteed to share the same embedded offset
# convention.
#
# Output:
#   - data_clean/species_climatic_affinity.csv: one row per species --
#     n_occurrences (after truncation), MAT_mean, MAP_mean (truncated means,
#     for use as HMSC traits)
#   - data_clean/species_climatic_affinity_resamples.csv: long format,
#     (species, rep 1-100, MAT_resample, MAP_resample) for CTI bootstrap
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(terra)

# 0a. Module directories ------------------------------------------------------------
# The CHELSA-bioclim download (00_download_chelsa_bioclim.R) lives in this
# same module, since it's downloaded specifically to be combined with GBIF
# here -- so there's a single module directory, not a separate one per data
# source.
occ_module_dir  <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", ...)
this_module_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
# Cross-module read, for the WFO synonym table only (see section 1b).
traits_dir   <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

N_RESAMPLES <- 100
set.seed(42)  # reproducible bootstrap resampling

# 1. Load cleaned North America occurrences ----------------------------------------
occ <- read_csv(
  occ_module_dir("data_clean", "gbif_occurrences_clean_north_america.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

cat("Clean occurrences (North America):", nrow(occ), "\n")
cat("Species (GBIF name, before synonym harmonization):", n_distinct(occ$species), "\n")

# 1b. Harmonize occurrence species names to the project's canonical name -----------
# spp_synonyms (columns final_species, synonym) is the same WFO synonym table
# used for TRY matching (02_build_trait_table_try.R) and for the LOO-
# validation alias fix (09_trait_imputation.R) -- reused here rather than
# rebuilt, so there's a single synonym source of truth for the whole project.
synonyms_rds_path <- traits_dir("tmp_files", "trait_matching_intermediate.rds")

if (file.exists(synonyms_rds_path)) {
  spp_synonyms <- readRDS(synonyms_rds_path)$spp_synonyms
  
  synonym_lookup <- spp_synonyms %>%
    distinct(synonym, final_species) %>%
    # A synonym mapped to more than one canonical species would make the
    # relabel ambiguous -- drop those (logged below) rather than guess.
    add_count(synonym, name = "n_targets") %>%
    filter(n_targets == 1) %>%
    select(synonym, final_species)
  
  n_ambiguous <- spp_synonyms %>% distinct(synonym, final_species) %>%
    add_count(synonym, name = "n_targets") %>% filter(n_targets > 1) %>%
    distinct(synonym) %>% nrow()
  if (n_ambiguous > 0) {
    cat("NOTE:", n_ambiguous, "synonym(s) map to more than one canonical species",
        "-- left un-relabeled (ambiguous), see spp_synonyms for detail.\n")
  }
  
  occ <- occ %>%
    left_join(synonym_lookup, by = c("species" = "synonym")) %>%
    mutate(species_gbif_raw = species,
           species = coalesce(final_species, species)) %>%
    select(-final_species)
  
  n_relabeled_records <- sum(occ$species != occ$species_gbif_raw)
  relabeled_species <- occ %>% filter(species != species_gbif_raw) %>%
    distinct(species_gbif_raw, species)
  
  cat("Synonyms harmonized:", n_relabeled_records, "occurrences relabelled,",
      nrow(relabeled_species), "species affected:\n")
  if (nrow(relabeled_species) > 0) print(relabeled_species, n = Inf)
  
  occ <- occ %>% select(-species_gbif_raw)
} else {
  cat("NOTE:", synonyms_rds_path, "not found -- run",
      "04_traits/02_build_trait_table_try.R first. Without this file,",
      "occurrences keep their GBIF names (no synonym",
      "harmonization).\n")
}

cat("Species (after synonym harmonization):", n_distinct(occ$species), "\n")

# 2. Extract CHELSA-bioclim bio1 (MAT) and bio12 (MAP) at occurrence points ---------
bio1_path  <- this_module_dir("data_raw", "chelsa_bioclim", "CHELSA_bio1_1981-2010_V.2.1.tif")
bio12_path <- this_module_dir("data_raw", "chelsa_bioclim", "CHELSA_bio12_1981-2010_V.2.1.tif")

if (!file.exists(bio1_path) || !file.exists(bio12_path)) {
  stop("CHELSA-bioclim GeoTIFFs not found -- run 00_download_chelsa_bioclim.R first.")
}

pts_vect <- vect(occ, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

cat("\nExtracting bio1 (MAT)...\n")
r_bio1 <- rast(bio1_path)
bio1_raw <- terra::extract(r_bio1, pts_vect, ID = FALSE)[[1]]

# Empirical Kelvin check.
cat("raw bio1 -- range:", round(range(bio1_raw, na.rm = TRUE), 2), "\n")
if (mean(bio1_raw, na.rm = TRUE) > 200) {
  cat("Range indicates Kelvin -- subtracting 273.15\n")
  occ$MAT <- bio1_raw - 273.15
} else {
  cat("Range indicates Celsius -- no conversion\n")
  occ$MAT <- bio1_raw
}
cat("MAT -- range:", round(range(occ$MAT, na.rm = TRUE), 2), "\n")

cat("\nExtracting bio12 (annual precipitation)...\n")
r_bio12 <- rast(bio12_path)
occ$MAP <- terra::extract(r_bio12, pts_vect, ID = FALSE)[[1]]
cat("MAP -- range:", round(range(occ$MAP, na.rm = TRUE), 2), "\n")

dir.create(this_module_dir("data_clean"), recursive = TRUE, showWarnings = FALSE)
write_csv(occ, this_module_dir("data_clean", "gbif_occurrences_with_bioclim.csv"))

# 3. Truncate to 5th-95th percentile per species -------------------------------------
truncate_5_95 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(x)
  q <- quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)
  x[x >= q[1] & x <= q[2]]
}

species_list <- occ %>% distinct(species) %>% pull(species)

affinity_summary <- vector("list", length(species_list))
resamples_list   <- vector("list", length(species_list))

for (i in seq_along(species_list)) {
  sp <- species_list[i]
  sp_data <- occ %>% filter(species == sp)
  
  mat_trunc <- truncate_5_95(sp_data$MAT)
  map_trunc <- truncate_5_95(sp_data$MAP)
  
  n_occ <- nrow(sp_data)
  
  affinity_summary[[i]] <- tibble(
    species = sp,
    n_occurrences = n_occ,
    n_truncated_MAT = length(mat_trunc),
    n_truncated_MAP = length(map_trunc),
    MAT_mean = if (length(mat_trunc) > 0) mean(mat_trunc) else NA_real_,
    MAP_mean = if (length(map_trunc) > 0) mean(map_trunc) else NA_real_
  )
  
  if (length(mat_trunc) > 0 && length(map_trunc) > 0) {
    resamples_list[[i]] <- tibble(
      species = sp,
      rep = 1:N_RESAMPLES,
      MAT_resample = sample(mat_trunc, N_RESAMPLES, replace = TRUE),
      MAP_resample = sample(map_trunc, N_RESAMPLES, replace = TRUE)
    )
  }
}

species_climatic_affinity <- bind_rows(affinity_summary)
species_climatic_affinity_resamples <- bind_rows(resamples_list)

n_no_affinity <- sum(is.na(species_climatic_affinity$MAT_mean))
cat("\n", n_no_affinity, "species without enough occurrences to compute a climatic affinity\n")

write_csv(species_climatic_affinity, this_module_dir("data_clean", "species_climatic_affinity.csv"))
write_csv(species_climatic_affinity_resamples,
          this_module_dir("data_clean", "species_climatic_affinity_resamples.csv"))

cat("\n========================================\n")
cat("GBIF + CHELSA INTEGRATION COMPLETE\n")
cat("========================================\n")
cat("Species processed:", nrow(species_climatic_affinity), "\n")
cat("Saved in pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/data_clean/\n")
cat("  - species_climatic_affinity.csv (single value per species, for HMSC)\n")
cat("  - species_climatic_affinity_resamples.csv (100 resamples per species, for CTI/CMI)\n")
