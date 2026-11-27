# ------------------------------------------------------------------------------
# 02_export_raw_surveys.R
# Module: pipeline/07_data_archive
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: copy the original survey files read by the cleaning scripts
# (pipeline/01_data_prep/01_cleaning/code/) into the data archive, so that the
# cleaning can be reproduced from the archive. Columns that are not used by the
# cleaning and should not be published are removed: observer names and
# initials, directions to the plots, photo file names, and free-text notes.
# Empty unnamed columns are also removed. All other values are kept exactly as
# recorded (every column is read and written as text, without trimming).
#
# Two steps:
#   1. INSPECT <- TRUE  : lists, for every file, the columns kept and removed.
#   2. INSPECT <- FALSE : writes the files to
#                         pipeline/07_data_archive/archive/raw_surveys/.
# Coordinates are rounded to COORD_DIGITS decimal places (and UTM coordinates
# to the nearest km) if COORD_DIGITS is not NA.
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

INSPECT      <- TRUE
COORD_DIGITS <- NA    # NA = as recorded; e.g. 2 to round to ~1 km

# Columns removed (exact names); any column without a name is also removed
DROP_COLUMNS <- c(
  # people
  "Observers", "RLB Notes", "EHL Notes",
  # plot location and photos
  "Directions to plot", "Camera", "FacingS", "FacingW", "FacingN", "FacingE",
  "Photos labeled", "2015 picture name",
  # free-text notes
  "Notes", "Notes 1", "Notes 2", "EnviromentNotes", "SiteDescription",
  "Representative Homogeneous notes", "Disturbance notes", "Comments",
  "T0CombineNotes", "ID notes", "Collection Notes", "Addition notes"
)
COORD_COLUMNS <- c("Latitude", "Longitude", "Lat 1", "Long 1", "Lat 2", "Long 2",
                   "Lat1", "Long1", "Long1_2", "Lat2", "Long2", "Long2_2")
UTM_COLUMNS   <- c("UTME", "UTMN")

raw_dir <- here("pipeline", "01_data_prep", "01_cleaning", "data", "data_raw")
files <- c(
  list.files(file.path(raw_dir, "01_raw_surveys"), pattern = "\\.csv$", full.names = TRUE),
  file.path(raw_dir, "02_earlier_cleaning", "MORA_PlotNames.csv")
)
files <- files[file.exists(files)]
out_dir <- here("pipeline", "07_data_archive", "archive", "raw_surveys")

round_text <- function(x, digits) {
  num <- suppressWarnings(as.numeric(x))
  ifelse(is.na(num), x, as.character(round(num, digits)))
}

for (f in files) {
  d <- read_csv(f, col_types = cols(.default = col_character()), na = character(),
                trim_ws = FALSE, name_repair = "minimal", progress = FALSE)
  nm <- names(d)
  drop <- which(nm %in% DROP_COLUMNS | trimws(nm) == "")
  keep <- setdiff(seq_along(nm), drop)
  if (INSPECT) {
    cat("\n==", basename(f), "--", nrow(d), "rows\n")
    cat("   removed:", if (length(drop)) paste(ifelse(trimws(nm[drop]) == "", "(unnamed)", nm[drop]), collapse = " | ") else "-", "\n")
    cat("   kept:   ", length(keep), "columns\n")
  } else {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    d <- d[, keep]
    if (!is.na(COORD_DIGITS)) {
      for (j in which(names(d) %in% COORD_COLUMNS)) d[[j]] <- round_text(d[[j]], COORD_DIGITS)
      for (j in which(names(d) %in% UTM_COLUMNS))   d[[j]] <- round_text(d[[j]], -3)
    }
    write_csv(d, file.path(out_dir, basename(f)), na = "")
    cat("written:", basename(f), "-", nrow(d), "rows,", ncol(d), "columns (removed", length(drop), ")\n")
  }
}
if (INSPECT) cat("\nIf the lists are correct, set INSPECT <- FALSE and run again.\n")