# ------------------------------------------------------------------------------
# 02b_gbif_visual_check.R
# Module: pipeline/01_data_prep/05_species_distribution/01_gbif_occurrences
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: EXPLORATORY / not part of the numbered pipeline output chain --
# visual sanity check of what 02_gbif_clean_occurrences.R did, for a random
# sample of species: raw vs. spatially-thinned occurrence points, and the
# convex hull used for the range-size trait. Uses base R plotting only (no
# ggplot2/sf) to keep this simple and robust to run interactively.
#
# Scope: North America only (see 02_gbif_clean_occurrences.R and the
# accompanying methods note -- the global scope was dropped).
#
# Change SEED or N_SPECIES below to look at a different random sample.
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)

# 1. Settings ---------------------------------------------------------------------
SEED       <- 123  # change this to see a different random sample of species
N_SPECIES  <- 6

# 2. Load data ----------------------------------------------------------------------
occ_raw <- read_csv(
  here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_raw",
       "gbif_occurrences_raw_north_america.csv"),
  show_col_types = FALSE
)
occ_clean <- read_csv(
  here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "data_clean",
       "gbif_occurrences_clean_north_america.csv"),
  show_col_types = FALSE
)
range_size <- read_csv(
  here("pipeline", "01_data_prep", "05_species_distribution", "01_gbif_occurrences", "outputs",
       "gbif_species_range_size_north_america.csv"),
  show_col_types = FALSE
)

# 3. Pick a random sample of species with a computable convex hull (>=3 pts) --------
set.seed(SEED)
species_sample <- range_size %>%
  filter(!is.na(range_size_km2)) %>%
  slice_sample(n = N_SPECIES) %>%
  pull(species)

cat("Selected species (SEED =", SEED, "):\n")
print(species_sample)

# 4. Plot: one panel per species, base R graphics only -------------------------------
n_row <- ceiling(sqrt(N_SPECIES))
n_col <- ceiling(N_SPECIES / n_row)

op <- par(mfrow = c(n_row, n_col), mar = c(2, 2, 2, 1))

for (sp in species_sample) {
  raw_pts   <- occ_raw   %>% filter(species == sp) %>% select(decimalLongitude, decimalLatitude) %>% drop_na()
  clean_pts <- occ_clean %>% filter(species == sp) %>% select(decimalLongitude, decimalLatitude) %>% drop_na()
  
  plot(raw_pts$decimalLongitude, raw_pts$decimalLatitude,
       col = adjustcolor("grey50", alpha.f = 0.3), pch = 16, cex = 0.5,
       xlab = "", ylab = "", main = paste0(sp, "\n(raw n=", nrow(raw_pts), ", clean n=", nrow(clean_pts), ")"),
       cex.main = 0.8)
  
  points(clean_pts$decimalLongitude, clean_pts$decimalLatitude,
         col = adjustcolor("firebrick", alpha.f = 0.5), pch = 16, cex = 0.5)
  
  if (nrow(clean_pts) >= 3) {
    hull_idx <- chull(clean_pts$decimalLongitude, clean_pts$decimalLatitude)
    hull_idx <- c(hull_idx, hull_idx[1])  # close the polygon
    lines(clean_pts$decimalLongitude[hull_idx], clean_pts$decimalLatitude[hull_idx],
          col = "steelblue", lwd = 1.5)
  }
}

par(op)

cat("\nGrey = raw (before cleaning), red = clean (thinned), blue line = convex hull.\n")
cat("The maps are drawn in the RStudio Plots pane.\n")
