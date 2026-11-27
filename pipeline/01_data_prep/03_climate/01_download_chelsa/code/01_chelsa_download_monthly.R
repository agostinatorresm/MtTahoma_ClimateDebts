# ------------------------------------------------------------------------------
# 01_chelsa_download_monthly.R
# Module: pipeline/01_data_prep/03_climate/01_download_chelsa
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: extract CHELSA monthly climate data (mean temperature, minimum
# temperature, precipitation) at MORA plot coordinates, for the FULL
# continuous monthly record available from CHELSA V2.1 (1979-01 to 2021-12),
# not just the past/present survey windows -- this feeds both:
#   - a full time-series description of climate change at the site
#     (02_explore_climate), and
#   - the past (~1975 survey) vs. present (~2015 survey) window summaries
#     used as candidate HMSC covariates, computed downstream from the
#     annual series this script's output feeds into.
#
# METHOD: reads CHELSA's monthly GeoTIFFs directly via GDAL's /vsicurl/
# streaming (terra::rast() on the file URL, then terra::extract() at the plot
# coordinates), matching the monthly archive path pattern already verified
# against CHELSA's live cloud storage listing:
#   https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly/<var>/<year>/CHELSA_<var>_<MM>_<year>_V.2.1.tif
#
# SCALE/OFFSET: CHELSA's GeoTIFFs carry scale (0.1) and offset (0) as embedded
# raster metadata, and terra::rast()/extract() applies that automatically on
# read (terra's scoff() mechanism: value <- value * scale + offset). So the
# values extract() returns are already in physical units (Kelvin for
# tas/tasmin, mm for pr) -- no manual "* 0.1" is needed. The only unit
# conversion actually required is Kelvin -> Celsius for the temperature
# variables.
#
# COVERAGE: CHELSA's monthly archive spans 1979-01 to 2021-12 (confirmed
# against the CHELSA V2.1 dataset description, chelsa-climate.org/models/chelsa).
# tas/tasmin have no data for 1979 itself (pr for 1979 is available) -- this
# script attempts every year in the requested range and records failures
# per (var, year) rather than assuming any single year is complete.
#
# EFFICIENCY / RESUME: each month's raster is read and extracted ONCE for all
# plots together (terra::extract() takes the whole set of plot coordinates in
# a single call per file), so the cost is 12 months x n_years x 3 variables
# file reads, not that multiplied by n_plots. Each (variable, year) is cached
# to its own CSV and skipped on re-run if already present, since the full
# range (1979-2021, 43 years x 3 variables x 12 months = ~1,550 file reads)
# is a long-running extraction -- safe to stop and resume.
#
# Output: one CSV per variable per year in
# data_raw/chelsa_monthly/<variable>/<variable>_<year>.csv, in long format
# (plot_id, year, month, value) -- already at monthly resolution, so
# 02_explore_climate/code/02_compute_annual_variables.R can go
# directly to annual/seasonal/PAS calculations without any further
# daily-to-monthly aggregation step.
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(terra)

# 0a. Module directory --------------------------------------------------------
module_dir <- function(...) here("pipeline", "01_data_prep", "03_climate", "01_download_chelsa", ...)

# 1. Load plot coordinates (and elevation, if present) ---------------------------
plot_info <- read_csv(
  here("pipeline", "01_data_prep", "01_cleaning", "outputs", "plot_info.csv"),
  show_col_types = FALSE
)

# Column names are matched case-insensitively against common variants rather
# than assumed exactly. lon/lat are required (script stops if missing);
# elevation is optional (script warns and continues without it if no matching
# column is found, since it isn't needed for the extraction itself, only for
# the downstream description/PCA).
find_col <- function(df, candidates, required = TRUE) {
  hit <- names(df)[match(tolower(candidates), tolower(names(df)))]
  hit <- hit[!is.na(hit)][1]
  if (is.na(hit)) {
    if (required) {
      stop("None of these columns was found: ", paste(candidates, collapse = ", "),
           ". Available columns: ", paste(names(df), collapse = ", "))
    }
    return(NA_character_)
  }
  hit
}

col_plot_id <- find_col(plot_info, c("plot_id", "plotid", "plot"))
col_lon     <- find_col(plot_info, c("lon", "long", "longitude"))
col_lat     <- find_col(plot_info, c("lat", "latitude"))
col_elev    <- find_col(plot_info, c("elevation", "elev", "elev_m", "elevation_m", "altitude", "alt"),
                        required = FALSE)

if (is.na(col_elev)) {
  cat("WARNING: no elevation column found in plot_info.csv -- coords_lookup.csv will not include it.\n")
  cat("Available columns:", paste(names(plot_info), collapse = ", "), "\n")
} else {
  cat("Elevation taken from column:", col_elev, "\n")
}

# plot_info.csv has one row per plot PER SURVEY PERIOD (a "period" column
# distinguishes "past"/"present"), not one row per plot -- the same physical
# plot's past and present coordinates can differ by GPS noise between the two
# surveys, so deduplicating by exact (lat, lon) would treat many single plots
# as two separate locations and query CHELSA twice per plot at two nearly-
# identical points. Deduplicating by plot_id instead, keeping the "past"
# (1975 establishment) row's coordinates as the plot's canonical location,
# gives one query point per real plot.
if ("period" %in% names(plot_info)) {
  plot_info <- plot_info %>%
    arrange(plot_id, period != "past")  # "past" row sorts first per plot_id
}

select_cols <- c("plot_id", "lon", "lat", if (!is.na(col_elev)) "elevation")
coords_lookup <- plot_info %>%
  select(
    plot_id = all_of(col_plot_id),
    lon = all_of(col_lon),
    lat = all_of(col_lat),
    !!!(if (!is.na(col_elev)) rlang::set_names(col_elev, "elevation") else NULL)
  ) %>%
  distinct(plot_id, .keep_all = TRUE)

cat("Unique plots:", nrow(coords_lookup), "\n")

# Sanity check: confirm the past/present coordinate difference for each plot
# is small (GPS noise between the two surveys), not a real difference in
# location -- which would mean this dedup step is discarding a genuinely
# different point rather than a duplicate.
if ("period" %in% names(plot_info)) {
  coord_spread <- plot_info %>%
    select(plot_id = all_of(col_plot_id), lon = all_of(col_lon), lat = all_of(col_lat)) %>%
    group_by(plot_id) %>%
    filter(n() > 1) %>%
    summarise(lon_range = max(lon, na.rm = TRUE) - min(lon, na.rm = TRUE),
              lat_range = max(lat, na.rm = TRUE) - min(lat, na.rm = TRUE),
              any_missing_coord = anyNA(lon) || anyNA(lat),
              .groups = "drop")
  
  n_missing <- sum(coord_spread$any_missing_coord)
  if (n_missing > 0) {
    cat(n_missing, "plot(s) with missing lat/lon in some row:\n")
    print(coord_spread$plot_id[coord_spread$any_missing_coord])
  }
  
  spread_ok <- coord_spread %>% filter(!any_missing_coord)
  if (nrow(spread_ok) > 0) {
    cat("Maximum past vs. present coordinate difference (degrees):\n")
    cat("  lon:", max(spread_ok$lon_range, na.rm = TRUE),
        " | lat:", max(spread_ok$lat_range, na.rm = TRUE), "\n")
    cat("(for reference, 0.01 degrees ~ 1 km at this latitude)\n")
  }
}

dir.create(module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
write_csv(coords_lookup, module_dir("outputs", "coords_lookup.csv"))

# terra::vect() preserves row order, so extract() results line up positionally
# with coords_lookup$plot_id.
pts_vect <- vect(coords_lookup, geom = c("lon", "lat"), crs = "EPSG:4326")

# 2. Full monthly record available from CHELSA V2.1 -------------------------------
all_years <- 1979:2021

cat("Extracting the full monthly series available in CHELSA V2.1:",
    min(all_years), "-", max(all_years), "\n")
cat("(tas/tasmin have no data for 1979 -- recorded as a failed year, the script continues)\n")

# 3. URL builder and extraction helper -------------------------------------------
chelsa_monthly_url <- function(var, year, month) {
  mm <- sprintf("%02d", month)
  paste0(
    "/vsicurl/https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/monthly/",
    var, "/", year, "/CHELSA_", var, "_", mm, "_", year, "_V.2.1.tif"
  )
}

extract_monthly <- function(pts_vect, var, year, month) {
  r <- rast(chelsa_monthly_url(var, year, month))
  vals <- terra::extract(r, pts_vect, ID = FALSE)[[1]]
  if (var %in% c("tas", "tasmin")) vals - 273.15 else vals
}

# 4. Extraction function (per variable, per year, with resume support) -----------
extract_chelsa_var_year <- function(var_code, years, coords_lookup, pts_vect, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (year in years) {
    out_file <- file.path(out_dir, paste0(var_code, "_", year, ".csv"))
    
    if (file.exists(out_file)) {
      cat("  [", var_code, year, "] already extracted, skipping\n")
      next
    }
    
    cat(">>> Extracting", var_code, year, "(12 months)...\n")
    start_time <- Sys.time()
    
    result <- tryCatch({
      monthly_vals <- map_dfr(1:12, function(m) {
        vals <- extract_monthly(pts_vect, var_code, year, m)
        tibble(plot_id = coords_lookup$plot_id, year = year, month = m, value = vals)
      })
      monthly_vals
    }, error = function(e) {
      cat("    ERROR in", var_code, year, ":", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(result)) {
      write_csv(result, out_file)
      elapsed <- round(as.numeric(Sys.time() - start_time, units = "secs"), 1)
      cat("    OK (", elapsed, "s )\n")
    }
    
    Sys.sleep(1)  # be polite to the CHELSA server
  }
}

# 5. Run extraction for the three variables needed --------------------------------
# tas    = monthly mean temperature   -> annual/seasonal MAT, PAS
# pr     = monthly precipitation      -> annual/seasonal MAP, PAS
# tasmin = monthly minimum temperature -> kept for consistency with the HMSC
#          covariate table (chelsa_tasmin); not required for PAS itself.

extract_chelsa_var_year("tas",    all_years, coords_lookup, pts_vect,
                        module_dir("data_raw", "chelsa_monthly", "tas"))
extract_chelsa_var_year("pr",     all_years, coords_lookup, pts_vect,
                        module_dir("data_raw", "chelsa_monthly", "pr"))
extract_chelsa_var_year("tasmin", all_years, coords_lookup, pts_vect,
                        module_dir("data_raw", "chelsa_monthly", "tasmin"))

cat("\n========================================\n")
cat("EXTRACTION COMPLETE\n")
cat("========================================\n")
cat("Files saved in pipeline/01_data_prep/03_climate/01_download_chelsa/data_raw/chelsa_monthly/{tas,pr,tasmin}/\n")
cat("Next step: 01_download_chelsa/code/02_compute_annual_variables.R\n")
