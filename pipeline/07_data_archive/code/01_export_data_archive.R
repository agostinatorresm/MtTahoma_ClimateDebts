# ------------------------------------------------------------------------------
# 01_export_data_archive.R
# Module: pipeline/07_data_archive
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: export the data needed to reproduce the analyses into one folder
# for the public data archive (Dryad), with clear file and column
# names. Only final tables are exported; intermediate and diagnostic files
# stay in the pipeline.
#
# Files written to pipeline/07_data_archive/archive/:
#   vegetation_overstory.csv     plot x survey x species: occurrence, stem
#                                count, basal area
#   vegetation_understory.csv    plot x survey x species: occurrence, cover
#   plots.csv                    one row per plot: location, topography,
#                                survey years, climate in both periods
#   species_climatic_affinity.csv            one row per species
#   species_climatic_affinity_resamples.csv  100 resamples per species
#   species_traits.csv           one row per layer x species, with source
#   species_pools.csv            one row per layer x species, pool membership
#   phylogeny.tre                Newick tree of all species
#
# Coordinates are exported at full precision unless COORD_DIGITS is set
# (e.g. 2 decimal places, about 1 km), in case exact plot locations cannot
# be published.
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(ape)

COORD_DIGITS <- NA   # NA = full precision; e.g. 2 to round to ~1 km

p   <- function(...) here("pipeline", ...)
out <- function(...) p("07_data_archive", "archive", ...)
dir.create(out(), recursive = TRUE, showWarnings = FALSE)

# 1. Survey years per plot -----------------------------------------------------
plot_info <- read_csv(p("01_data_prep", "01_cleaning", "outputs", "plot_info.csv"), show_col_types = FALSE)
past <- plot_info %>% filter(period == "past") %>% distinct(plot_id, .keep_all = TRUE)
present_year <- read_csv(p("01_data_prep", "01_cleaning", "outputs", "present_census_year.csv"),
                         show_col_types = FALSE)
survey_year <- bind_rows(
  past %>% transmute(plot_id, period = "past", year),
  present_year %>% transmute(plot_id, period = "present", layer = stratum, year)
)
year_for <- function(layer_name) {
  survey_year %>% filter(is.na(layer) | layer == layer_name) %>% select(plot_id, period, survey_year = year)
}

# 2. Vegetation tables -----------------------------------------------------------
overstory <- read_csv(p("01_data_prep", "01_cleaning", "outputs", "overstory_long.csv"), show_col_types = FALSE) %>%
  select(-year) %>%
  left_join(year_for("overstory"), by = c("plot_id", "period")) %>%
  transmute(plot_id, period, survey_year, species, occurrence, stem_count = count,
            basal_area_m2 = basal_area_spp_m2, plot_area_m2 = plot_area) %>%
  arrange(plot_id, period, species)
write_csv(overstory, out("vegetation_overstory.csv"))

understory <- read_csv(p("01_data_prep", "01_cleaning", "outputs", "understory_long.csv"), show_col_types = FALSE) %>%
  select(-year) %>%
  left_join(year_for("understory"), by = c("plot_id", "period")) %>%
  transmute(plot_id, period, survey_year, species, occurrence, cover_percent = cover) %>%
  arrange(plot_id, period, species)
write_csv(understory, out("vegetation_understory.csv"))

# 3. Plots -----------------------------------------------------------------------
climate <- read_csv(p("01_data_prep", "03_climate", "02_explore_climate", "data_clean",
                      "climate_past_present_wide.csv"), show_col_types = FALSE) %>%
  select(plot_id, MAT_past, MAT_present, MAT_min_past, MAT_min_present,
         MAP_past, MAP_present, PAS_past, PAS_present,
         n_years_temp_past, n_years_temp_present, n_years_precip_past, n_years_precip_present)

present_wide <- present_year %>%
  transmute(plot_id, layer = stratum, year, selection) %>%
  pivot_wider(names_from = layer, values_from = c(year, selection), names_glue = "{.value}_present_{layer}")

round_coord <- function(x) if (is.na(COORD_DIGITS)) x else round(x, COORD_DIGITS)
plots <- past %>%
  transmute(plot_id, latitude = round_coord(lat), longitude = round_coord(long),
            elevation_m = elevation, aspect_deg = aspect, slope_deg = slope,
            plot_area_past_m2 = plotarea, year_past = year) %>%
  left_join(present_wide, by = "plot_id") %>%
  left_join(climate, by = "plot_id") %>%
  filter(plot_id %in% c(overstory$plot_id, understory$plot_id)) %>%
  arrange(plot_id)
write_csv(plots, out("plots.csv"))

# 4. Species tables ------------------------------------------------------------------
aff_dir <- function(...) p("01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", "data_clean", ...)
read_csv(aff_dir("species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  write_csv(out("species_climatic_affinity.csv"))
read_csv(aff_dir("species_climatic_affinity_resamples.csv"), show_col_types = FALSE) %>%
  write_csv(out("species_climatic_affinity_resamples.csv"))

traits <- bind_rows(
  read_csv(p("01_data_prep", "04_traits", "outputs", "overstory_traits.csv"), show_col_types = FALSE) %>%
    mutate(layer = "overstory"),
  read_csv(p("01_data_prep", "04_traits", "outputs", "understory_traits.csv"), show_col_types = FALSE) %>%
    mutate(layer = "understory")
) %>%
  select(layer, species, any_of(c("seed_mass", "seed_mass_src", "seed_mass_n_records",
                                  "SLA", "SLA_src", "SLA_n_records",
                                  "wood_density", "wood_density_src", "wood_density_n_records",
                                  "woodiness", "woodiness_src", "woodiness_n_records",
                                  "growth_form", "growth_form_src", "growth_form_n_records",
                                  "longevity", "longevity_src", "longevity_n_records"))) %>%
  arrange(layer, species)
write_csv(traits, out("species_traits.csv"))

pool_dir <- function(...) p("01_data_prep", "02_species_pool", "data_clean", ...)
pools <- bind_rows(
  read_csv(pool_dir("species_pool_candidates_overstory.csv"), show_col_types = FALSE) %>% mutate(layer = "overstory"),
  read_csv(pool_dir("species_pool_candidates_understory.csv"), show_col_types = FALSE) %>% mutate(layer = "understory")
) %>%
  transmute(layer, species, n_plots_past, n_plots_present,
            shared_pool = included_broad, frequent_pool = included_strict) %>%
  arrange(layer, species)
write_csv(pools, out("species_pools.csv"))

# 5. Phylogeny -------------------------------------------------------------------------
load(p("01_data_prep", "04_traits", "outputs", "MORA_combined_phylo_tree.rda"))   # object: phy
write.tree(phy, out("phylogeny.tre"))

# 6. Checks ------------------------------------------------------------------------------
stopifnot(!any(is.na(overstory$survey_year)), !any(is.na(understory$survey_year)))
cat("Files written to", out(), "\n")
print(tibble(file = list.files(out()),
             size_kb = round(file.size(list.files(out(), full.names = TRUE)) / 1024)))
cat("Plots:", nrow(plots), "| overstory species:", n_distinct(overstory$species),
    "| understory species:", n_distinct(understory$species),
    "| tips in phylogeny:", length(phy$tip.label), "\n")