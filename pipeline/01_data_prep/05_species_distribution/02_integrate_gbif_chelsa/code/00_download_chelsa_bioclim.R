# ------------------------------------------------------------------------------
# 00_download_chelsa_bioclim.R
# Module: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: download the two CHELSA-bioclim 1981-2010 GeoTIFFs needed for the
# species climatic-affinity extraction (bio1 = mean annual temperature,
# bio12 = annual precipitation), cached locally in this module's data_raw/.
# Independent of the GBIF occurrence download -- can be run at any time,
# including while a GBIF occ_download() is still processing, since it only
# needs network access to CHELSA's own cloud storage.
#
# Each file is global, ~1km resolution, so this can take a while and the
# files are a few hundred MB each -- run once, the result is reused by
# 01_combine_gbif_chelsa.R as long as the files stay in
# data_raw/chelsa_bioclim/.
#
# Output:
#   - data_raw/chelsa_bioclim/CHELSA_bio1_1981-2010_V.2.1.tif
#   - data_raw/chelsa_bioclim/CHELSA_bio12_1981-2010_V.2.1.tif
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(here)

# 0a. R's default per-download timeout (60s) is too short for these files
# (bio12 alone is ~625 MB) -- raise it for this session.
options(timeout = 1200)  # 20 minutes

# 1. Setup ---------------------------------------------------------------------
this_module_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)

chelsa_bioclim_dir <- this_module_dir("data_raw", "chelsa_bioclim")
dir.create(chelsa_bioclim_dir, recursive = TRUE, showWarnings = FALSE)

download_bioclim <- function(var) {
  url <- paste0(
    "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/",
    "CHELSA_", var, "_1981-2010_V.2.1.tif"
  )
  dest <- file.path(chelsa_bioclim_dir, paste0("CHELSA_", var, "_1981-2010_V.2.1.tif"))
  if (file.exists(dest)) {
    cat(var, "already downloaded in", dest, "\n")
    return(invisible(NULL))
  }
  cat("Downloading", var, "(global file, several hundred MB)...\n")
  result <- tryCatch({
    download.file(url, dest, mode = "wb")
    TRUE
  }, error = function(e) {
    cat("  ERROR downloading", var, ":", conditionMessage(e), "\n")
    # Delete a partial/corrupt file so the next run doesn't mistake it for a
    # complete download (file.exists() above would otherwise skip it).
    if (file.exists(dest)) {
      file.remove(dest)
      cat("  Partial file deleted -- run this script again to retry.\n")
    }
    FALSE
  })
  if (isTRUE(result)) cat("  OK --", dest, "\n")
}

# 2. Download bio1 (MAT) and bio12 (annual precipitation) --------------------------
download_bioclim("bio1")
download_bioclim("bio12")

cat("\n========================================\n")
cat("CHELSA-BIOCLIM DOWNLOAD COMPLETE\n")
cat("========================================\n")
cat("Files saved in pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/data_raw/chelsa_bioclim/\n")
cat("Next step: 01_combine_gbif_chelsa.R, once the clean GBIF occurrences are available\n")
cat("(it reuses these files instead of downloading them again).\n")
