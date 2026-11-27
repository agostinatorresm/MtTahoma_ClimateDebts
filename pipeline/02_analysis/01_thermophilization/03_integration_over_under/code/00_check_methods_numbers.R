# ------------------------------------------------------------------------------
# 00_check_methods_numbers.R
# Module: pipeline (diagnostic)
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: print every number quoted in the manuscript Methods that is not
# already reported by another script's console output, so each value can be
# checked against the files the pipeline actually produced. Diagnostic only;
# nothing is written to disk.
#
# Sections:
#   1. Species pools (broad / strict / broad with full climatic coverage)
#   2. GBIF records per pool species (n_occurrences after cleaning)
#   3. Plots and survey interval entering CTI/CMI, per layer
#   4. Trait coverage per layer: direct / literature / genus-imputed /
#      unresolved, by trait
#   5. Leave-one-out validation files, if present
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

prep_dir     <- function(...) here("pipeline", "01_data_prep", ...)
thermo_dir   <- function(...) here("pipeline", "02_analysis", "01_thermophilization", ...)

read_if_exists <- function(path) {
  if (!file.exists(path)) {
    cat("  [missing]", path, "\n")
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE)
}

section <- function(title) cat("\n==============================\n", title, "\n==============================\n", sep = "")

# 1. Species pools ---------------------------------------------------------------
section("1. SPECIES POOLS")

pool_broad  <- read_if_exists(prep_dir("02_species_pool", "data_clean", "species_pool_broad.csv"))
pool_strict <- read_if_exists(prep_dir("02_species_pool", "data_clean", "species_pool_strict.csv"))

if (!is.null(pool_broad))  { cat("BROAD  (species x stratum):\n");  print(count(pool_broad,  stratum)) }
if (!is.null(pool_strict)) { cat("STRICT (species x stratum):\n"); print(count(pool_strict, stratum)) }

full_cov <- list(
  overstory  = read_if_exists(prep_dir("02_species_pool", "outputs", "overstory_species_pool_broad_full_clim_cov.csv")),
  understory = read_if_exists(prep_dir("02_species_pool", "outputs", "understory_species_pool_broad_full_clim_cov.csv"))
)
cat("\nBROAD with full climatic coverage (used for CTI/CMI):\n")
iwalk(full_cov, ~ if (!is.null(.x)) cat("  ", .y, ":", nrow(.x), "species\n"))

# Broad-pool species lost to the full-coverage requirement
if (!is.null(pool_broad)) {
  for (st in c("overstory", "understory")) {
    if (is.null(full_cov[[st]])) next
    lost <- setdiff(pool_broad %>% filter(stratum == st) %>% pull(species), full_cov[[st]]$species)
    cat("\n", st, ": broad-pool species without full climatic coverage (", length(lost), "):\n", sep = "")
    if (length(lost) > 0) print(lost)
  }
}

# 2. GBIF records per pool species ------------------------------------------------
section("2. GBIF RECORDS PER POOL SPECIES (after cleaning and thinning)")

affinity <- read_if_exists(prep_dir("05_species_distribution", "02_integrate_gbif_chelsa", "data_clean", "species_climatic_affinity.csv"))

if (!is.null(affinity)) {
  for (st in names(full_cov)) {
    if (is.null(full_cov[[st]])) next
    sub <- affinity %>% filter(species %in% full_cov[[st]]$species)
    cat("\n", st, " (n species = ", nrow(sub), ")\n", sep = "")
    cat("  n_occurrences: min", min(sub$n_occurrences), "| median", median(sub$n_occurrences),
        "| max", max(sub$n_occurrences), "\n")
    cat("  species with fewer than 20 records:", sum(sub$n_occurrences < 20), "\n")
    cat("  species with fewer than 50 records:", sum(sub$n_occurrences < 50), "\n")
    cat("  n_truncated_MAT: min", min(sub$n_truncated_MAT), "| median", median(sub$n_truncated_MAT), "\n")
  }
  cat("\nSpecies in the affinity file with no MAT_mean (no usable records):",
      sum(is.na(affinity$MAT_mean)), "\n")
}

# 3. Plots and survey interval -------------------------------------------------------
section("3. PLOTS AND SURVEY INTERVAL ENTERING CTI / CMI")

rates_over  <- read_if_exists(thermo_dir("01_overstory",  "outputs", "overstory_thermo_xero_rates.csv"))
rates_under <- read_if_exists(thermo_dir("02_understory", "outputs", "understory_thermo_xero_rates.csv"))

summarise_rates <- function(df, primary_weight, label) {
  if (is.null(df)) return(invisible(NULL))
  d <- df %>% filter(weight_var == primary_weight)
  cat("\n", label, " (weight_var = ", primary_weight, ")\n", sep = "")
  cat("  plots with rates in both periods:", n_distinct(d$plot_id), "\n")
  cat("  year_gap_decades: min", round(min(d$year_gap_decades, na.rm = TRUE), 2),
      "| mean", round(mean(d$year_gap_decades, na.rm = TRUE), 3),
      "| max", round(max(d$year_gap_decades, na.rm = TRUE), 2), "\n")
  cat("  distribution of year_gap_decades (rounded to 0.1):\n")
  print(count(d, gap = round(year_gap_decades, 1)))
  cat("  mean rate_CTI (degC/decade):", round(mean(d$rate_CTI, na.rm = TRUE), 4),
      "| mean rate_CMI (mm/decade):", round(mean(d$rate_CMI, na.rm = TRUE), 2), "\n")
}

summarise_rates(rates_over,  "basal_area_spp_perc", "OVERSTORY")
summarise_rates(rates_under, "cover",               "UNDERSTORY")

# LOESS settings used by the lag scripts (see each layer's 02_thermo_xero_elevation_lag.R)
cat("\nLOESS span used for the elevation curves and the space-for-time relationship: 0.75\n")

# 4. Trait coverage ------------------------------------------------------------------
section("4. TRAIT COVERAGE (final trait tables)")

classify_src <- function(src) {
  case_when(
    is.na(src)                                   ~ "unresolved",
    src == "unresolved"                          ~ "unresolved",
    str_detect(str_to_lower(src), "literature|lit:|mitchell|melton") ~ "literature",
    str_detect(str_to_lower(src), "genus")       ~ "genus-imputed",
    TRUE                                         ~ "direct (TRY and/or GIFT)"
  )
}

trait_tables <- list(
  overstory  = read_if_exists(prep_dir("04_traits", "outputs", "overstory_traits.csv")),
  understory = read_if_exists(prep_dir("04_traits", "outputs", "understory_traits.csv"))
)

for (st in names(trait_tables)) {
  tt <- trait_tables[[st]]
  if (is.null(tt)) next
  # Restrict to the broad pool with full climatic coverage, as in the analyses
  if (!is.null(full_cov[[st]])) tt <- tt %>% filter(species %in% full_cov[[st]]$species)
  
  cat("\n", st, " (n species = ", nrow(tt), ")\n", sep = "")
  cov <- tt %>%
    select(species, ends_with("_src")) %>%
    pivot_longer(-species, names_to = "trait", values_to = "src") %>%
    mutate(trait = str_remove(trait, "_src$"), class = classify_src(src)) %>%
    count(trait, class) %>%
    pivot_wider(names_from = class, values_from = n, values_fill = 0)
  print(cov, width = Inf)
  
  cat("\nSpecies still unresolved, by trait:\n")
  unresolved <- tt %>%
    select(species, ends_with("_src")) %>%
    pivot_longer(-species, names_to = "trait", values_to = "src") %>%
    filter(classify_src(src) == "unresolved") %>%
    mutate(trait = str_remove(trait, "_src$")) %>%
    group_by(trait) %>% summarise(n = n(), species = paste(species, collapse = "; "), .groups = "drop")
  print(unresolved, width = Inf)
  
  cat("\nSpecies filled from literature:\n")
  lit <- tt %>%
    select(species, ends_with("_src")) %>%
    pivot_longer(-species, names_to = "trait", values_to = "src") %>%
    filter(classify_src(src) == "literature")
  print(lit, width = Inf)
}

# 5. Leave-one-out validation outputs -------------------------------------------------
section("5. LEAVE-ONE-OUT VALIDATION FILES")

loo_files <- list.files(prep_dir("04_traits"), pattern = "loo|validation", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
loo_files <- loo_files[str_detect(loo_files, "\\.csv$")]

if (length(loo_files) == 0) {
  cat("No LOO csv found under 04_traits; re-run 09_trait_imputation.R and read the\n",
      "printed R2 / accuracy / n values from its console output.\n")
} else {
  for (f in loo_files) {
    cat("\n", basename(f), "\n", sep = "")
    print(read_csv(f, show_col_types = FALSE), n = 30, width = Inf)
  }
}

cat("\n========================================\n")
cat("METHODS NUMBERS CHECK -- COMPLETE\n")
cat("========================================\n")
