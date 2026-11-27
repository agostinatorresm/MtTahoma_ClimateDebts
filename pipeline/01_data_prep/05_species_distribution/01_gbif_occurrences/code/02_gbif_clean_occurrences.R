# ------------------------------------------------------------------------------
# 02_gbif_clean_occurrences.R
# Module: pipeline/01_data_prep/05_species_distribution/01_gbif_occurrences
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: clean the raw GBIF occurrence records (downloaded by
# 01_gbif_download_occurrences.R; year <= 1985) before they are used to
# compute species' climatic affinities and range size, following the
# spatial-cleaning approach in Lynn et al. (2021):
#   1. Drop records with reported spatial uncertainty >= 100 km (records with
#      no reported uncertainty are kept, not dropped -- see note below).
#   2. Spatially thin to one occurrence per species per 0.0083-degree grid
#      cell (~1 km, matching CHELSA/ClimateNA resolution), to reduce bias
#      from spatially clustered records (e.g. citizen-science hotspots or
#      herbarium collecting localities revisited repeatedly).
#   3. Compute each species' range size as the area of the convex hull of its
#      cleaned occurrences, in a Lambert azimuthal equal-area projection
#      centered on North America (so area is not distorted by using raw
#      lon/lat degrees) -- a candidate range-size trait, kept separate from
#      the thermal/moisture climatic affinities computed downstream.
#
# SCOPE: "global" or "north_america", matching the raw download(s) of
# 01_gbif_download_occurrences.R to clean. Records from the supplementary
# download (01b_gbif_download_missing_taxa.R) are added when present.
#
# Missing coordinateUncertaintyInMeters: many GBIF records, especially
# older, pre-GPS records (all records here, since the dataset is restricted
# to year <= 1985), have no reported uncertainty rather than a large one.
# Treating "not reported" as "imprecise" would discard a large fraction of
# the historical records, so only records with an explicit uncertainty
# >= 100 km are dropped; records with NA uncertainty are kept and counted
# separately.
#
# Output (tagged with SCOPE):
#   - data_clean/gbif_occurrences_clean_<SCOPE>.csv: cleaned, spatially
#     thinned occurrences (one row per species per grid cell)
#   - outputs/gbif_species_range_size_<SCOPE>.csv: per species, convex hull
#     area (km2) from the cleaned occurrences, plus the point count it was
#     computed from (NA for species with < 3 points, which cannot form a
#     hull)
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(sf)

# 0a. Scope for this run -----------------------------------------------------------
SCOPE <- "north_america"  # "global" or "north_america" -- match the raw file to clean
stopifnot(SCOPE %in% c("global", "north_america"))
cat("SCOPE =", SCOPE, "\n")

# 1. Load raw occurrences ---------------------------------------------------------
raw_path <- here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw",
                 paste0("gbif_occurrences_raw_", SCOPE, ".csv"))

if (!file.exists(raw_path)) {
  stop(basename(raw_path), " not found -- run 01_gbif_download_occurrences.R first with SCOPE = \"",
       SCOPE, "\".")
}

occurrences_raw <- read_csv(raw_path, show_col_types = FALSE)

# Supplementary download (01b_gbif_download_missing_taxa.R): project species
# not resolved at species rank by 01, already labelled with the project name
supp_path <- sub("\\.csv$", "_supplement.csv", raw_path)
if (file.exists(supp_path)) {
  occurrences_supp <- read_csv(supp_path, col_types = cols(.default = col_character())) %>%
    select(any_of(names(occurrences_raw))) %>%
    type_convert(col_types = spec(occurrences_raw)) %>%
    mutate(across(any_of("gbifID"), as.character))
  occurrences_raw <- occurrences_raw %>%
    mutate(across(any_of("gbifID"), as.character)) %>%
    filter(!species %in% unique(occurrences_supp$species)) %>%
    bind_rows(occurrences_supp)
  cat("Supplementary records added:", nrow(occurrences_supp), "for",
      n_distinct(occurrences_supp$species), "species\n")
}

cat("Raw occurrences (", SCOPE, "):", nrow(occurrences_raw), "\n")

# 2. Spatial uncertainty filter ----------------------------------------------------
UNCERTAINTY_THRESHOLD_M <- 100000  # 100 km, in meters (matches coordinateUncertaintyInMeters units)

n_high_uncertainty <- sum(occurrences_raw$coordinateUncertaintyInMeters >= UNCERTAINTY_THRESHOLD_M,
                          na.rm = TRUE)
n_missing_uncertainty <- sum(is.na(occurrences_raw$coordinateUncertaintyInMeters))

cat(n_high_uncertainty, "records with reported spatial uncertainty >= 100 km (removed)\n")
cat(n_missing_uncertainty, "records without reported uncertainty (",
    round(100 * n_missing_uncertainty / nrow(occurrences_raw), 1),
    "% of total -- kept, see the script header)\n")

occurrences_precise <- occurrences_raw %>%
  filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < UNCERTAINTY_THRESHOLD_M)

cat("Occurrences after the spatial-uncertainty filter:", nrow(occurrences_precise), "\n")

# 3. Spatial thinning: one occurrence per species per 0.0083-degree grid cell -------
GRID_SIZE_DEG <- 0.0083  # ~1 km, matches CHELSA/ClimateNA resolution

set.seed(42)  # slice_sample below is otherwise non-reproducible across runs

occurrences_thinned <- occurrences_precise %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  mutate(
    grid_lon = floor(decimalLongitude / GRID_SIZE_DEG),
    grid_lat = floor(decimalLatitude / GRID_SIZE_DEG)
  ) %>%
  group_by(species, grid_lon, grid_lat) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  select(-grid_lon, -grid_lat)

cat("Occurrences after spatial thinning (1 per species per 0.0083° cell):",
    nrow(occurrences_thinned), "\n")

dir.create(here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_clean"),
           recursive = TRUE, showWarnings = FALSE)
write_csv(occurrences_thinned,
          here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_clean",
               paste0("gbif_occurrences_clean_", SCOPE, ".csv")))

# 4. Range size: convex hull area per species, in a North America LAEA projection --
# ESRI:102008 (North America Lambert Azimuthal Equal-Area, centered 45N/100W)
# equivalent proj4 string, used directly so this doesn't depend on the EPSG/
# ESRI code being present in the local PROJ database. Used for both scopes so
# range-size values stay comparable -- a projection centered on North America
# is not distorted for North American points; for the global scope it is
# only used as a common, science-standard equal-area projection, not a claim
# that the resulting hull is itself restricted to North America.
NA_LAEA_CRS <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

occ_sf <- occurrences_thinned %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(NA_LAEA_CRS)

compute_hull_area_km2 <- function(pts_sf) {
  n <- nrow(pts_sf)
  if (n < 3) return(NA_real_)
  hull <- st_convex_hull(st_union(pts_sf))
  as.numeric(st_area(hull)) / 1e6  # m2 -> km2
}

range_size <- occ_sf %>%
  group_by(species) %>%
  group_modify(~ tibble(
    n_points_used = nrow(.x),
    range_size_km2 = compute_hull_area_km2(.x)
  )) %>%
  ungroup() %>%
  st_drop_geometry()

n_no_hull <- sum(is.na(range_size$range_size_km2))
cat("\n", n_no_hull, "species with < 3 points after cleaning -- no convex hull (range_size_km2 = NA)\n")

dir.create(here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "outputs"),
           recursive = TRUE, showWarnings = FALSE)
write_csv(range_size,
          here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "outputs",
               paste0("gbif_species_range_size_", SCOPE, ".csv")))

cat("\n========================================\n")
cat("CLEANING COMPLETE (SCOPE =", SCOPE, ")\n")
cat("========================================\n")
cat("Clean occurrences saved in data_clean/gbif_occurrences_clean_", SCOPE, ".csv\n", sep = "")
cat("Range size (convex hull) per species saved in outputs/gbif_species_range_size_", SCOPE, ".csv\n", sep = "")
if (SCOPE == "global") {
  cat("Next step: run this script with SCOPE <- \"north_america\".\n")
} else {
  cat("Next step: 02b_gbif_visual_check.R (visual check).\n")
}
