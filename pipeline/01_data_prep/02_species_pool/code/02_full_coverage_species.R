# ------------------------------------------------------------------------------
# 02_full_coverage_species.R
# Module: pipeline/01_data_prep/02_species_pool
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Adds climatic-affinity and functional-trait coverage to the species pools
# defined in 01_define_species_pool.R and writes the species lists used by
# each analysis:
#   thermophilization (CTI/CMI): species with climatic affinity (MAT and MAP);
#       functional traits are not needed.
#   HMSC: species with climatic affinity AND every functional trait in the
#       trait matrix of the layer (HMSC needs a complete species x trait
#       table).
# Coverage is assessed for the shared pool; frequent-pool lists are the
# subset of frequent species.
#
# Functional traits required per layer: the traits entering TrFormula in
# 02_analysis/02_hmsc/code/01_prepare_hmsc_inputs.R (traits_model, without
# MAT/MAP). The two settings must match, so that the HMSC pool is not
# restricted by a trait the models do not use.
#   overstory:  seed_mass, wood_density
#   understory: SLA, seed_mass, woodiness
# growth_form is reported for the understory but not required.
#
# Spore-bearing plants (ferns, horsetails, clubmosses) have no seeds: seed
# mass is undefined for them, not missing, so any value in the trait table
# (spore mass or a genus imputation) is set to NA before coverage is
# assessed. Where seed mass is a required trait, these species are therefore
# outside the HMSC set; they remain in the thermophilization set.
#
# Trait values are read from the final trait tables (04_traits/outputs/),
# where a value is either a species-level record, a genus-level imputation
# or a literature value, and NA only when unresolved.
#
# Names: the trait and climatic-affinity tables use current accepted names;
# they are translated to the names used in the vegetation tables and pools
# (name_to_pool) before matching.
#
# Inputs : data_clean/<stratum>_final_list_broad.csv / _strict.csv
#          04_traits/outputs/<stratum>_traits.csv
#          05_species_distribution/02_integrate_gbif_chelsa/data_clean/
#            species_climatic_affinity.csv
# Outputs (data_clean/):
#   <stratum>_trait_climatic_coverage.csv  one row per shared-pool species:
#       resolved_<trait>, MAT_mean, MAP_mean, functional_ok, climatic_ok,
#       full_coverage_hmsc, in_broad_pool, in_strict_pool, missing_traits
#   <stratum>_species_full_coverage_hmsc.csv  shared pool, functional + climatic
#   <stratum>_species_full_climatic.csv       shared pool, climatic
# Outputs (outputs/):
#   <stratum>_species_pool_<broad|strict>_full_clim_cov.csv
#   <stratum>_species_pool_<broad|strict>_full_func_clim_cov.csv
#   species_pool_coverage_summary.csv
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

traits_dir    <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)
climaffinity_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution",
                                       "02_integrate_gbif_chelsa", "data_clean", ...)
this_module_dir  <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)

dir.create(this_module_dir("data_clean"), recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("outputs"),    recursive = TRUE, showWarnings = FALSE)

strata <- c("overstory", "understory")

# Functional traits required for the HMSC pool (must match traits_model in
# 01_prepare_hmsc_inputs.R)
functional_traits <- list(
  overstory  = c("seed_mass", "wood_density"),
  understory = c("SLA", "seed_mass", "woodiness")
)
# Reported but not required
optional_traits <- list(
  overstory  = character(0),
  understory = c("growth_form")
)

# Accepted names in the trait / affinity tables -> names in the pools
name_to_pool <- c("Mahonia nervosa"            = "Berberis nervosa",
                  "Mahonia aquifolium"         = "Berberis aquifolium",
                  "Trautvetteria carolinensis" = "Trautvetteria caroliniensis",
                  "Hypopitys monotropa"        = "Monotropa hypopitys",
                  "Mycelis muralis"            = "Lactuca muralis",
                  "Ozomelis trifida"           = "Mitella trifida",
                  "Xanthocyparis nootkatensis" = "Callitropsis nootkatensis",
                  "Brewerimitella breweri"     = "Mitella breweri",
                  "Mitellastra caulescens"     = "Mitella caulescens",
                  "Pectiantia pentandra"       = "Mitella pentandra",
                  "Dodecatheon jeffreyanum"    = "Primula jeffreyi")
to_pool_name <- function(df) df %>% mutate(species = recode(species, !!!name_to_pool))

# Spore-bearing genera (seed mass not defined)
spore_genera <- c("Adiantum", "Asplenium", "Athyrium", "Blechnum", "Cystopteris",
                  "Dryopteris", "Equisetum", "Gymnocarpium", "Lycopodium",
                  "Polypodium", "Polystichum", "Pteridium")

read_species_col <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  if ("species" %in% names(df)) df$species else df[[1]]
}

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------
# Climatic affinity: one row per pool name; where an accepted name and its
# pool name both have a record, the one with MAT and MAP is kept
affinity <- read_csv(climaffinity_dir("species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean) %>%
  to_pool_name() %>%
  arrange(species, is.na(MAT_mean) | is.na(MAP_mean)) %>%
  distinct(species, .keep_all = TRUE)

traits <- map(set_names(strata), ~
                read_csv(traits_dir("outputs", paste0(.x, "_traits.csv")), show_col_types = FALSE) %>%
                to_pool_name() %>%
                select(-any_of(c("MAT_mean", "MAP_mean"))) %>%
                distinct(species, .keep_all = TRUE) %>%
                mutate(across(any_of("seed_mass"), ~ if_else(word(species, 1) %in% spore_genera, NA_real_, .x))))

pools <- map(set_names(strata), ~ list(
  broad  = read_species_col(this_module_dir("data_clean", paste0(.x, "_final_list_broad.csv"))),
  strict = read_species_col(this_module_dir("data_clean", paste0(.x, "_final_list_strict.csv")))
))

# ------------------------------------------------------------------------------
# 2. Coverage table (shared pool)
# ------------------------------------------------------------------------------
build_coverage_table <- function(stratum) {
  tr     <- traits[[stratum]]
  broad  <- pools[[stratum]]$broad
  strict <- pools[[stratum]]$strict
  req    <- functional_traits[[stratum]]
  opt    <- intersect(optional_traits[[stratum]], names(tr))
  
  missing_cols <- setdiff(req, names(tr))
  if (length(missing_cols) > 0) {
    stop(stratum, " trait table lacks required trait columns: ", paste(missing_cols, collapse = ", "))
  }
  no_row <- setdiff(broad, tr$species)
  if (length(no_row) > 0) {
    message(str_to_title(stratum), ": shared-pool species without a trait-table row: ",
            paste(no_row, collapse = ", "))
  }
  not_in_pool <- setdiff(tr$species, broad)
  if (length(not_in_pool) > 0) {
    message(str_to_title(stratum), ": trait-table species not in the shared pool (ignored): ",
            paste(not_in_pool, collapse = ", "))
  }
  
  tibble(species = broad) %>%
    left_join(tr %>% select(species, all_of(c(req, opt))), by = "species") %>%
    left_join(affinity, by = "species") %>%
    mutate(across(all_of(c(req, opt)), ~ !is.na(.x), .names = "resolved_{.col}")) %>%
    mutate(
      functional_ok      = if_all(all_of(paste0("resolved_", req))),
      climatic_ok        = !is.na(MAT_mean) & !is.na(MAP_mean),
      full_coverage_hmsc = functional_ok & climatic_ok,
      in_broad_pool      = TRUE,
      in_strict_pool     = species %in% strict
    ) %>%
    rowwise() %>%
    mutate(missing_traits = paste(c(req[!c_across(all_of(paste0("resolved_", req)))],
                                    if (!climatic_ok) "climatic_affinity"),
                                  collapse = ", ")) %>%
    ungroup() %>%
    mutate(missing_traits = na_if(missing_traits, "")) %>%
    select(species, functional_ok, climatic_ok, full_coverage_hmsc, in_broad_pool,
           in_strict_pool, missing_traits, MAT_mean, MAP_mean, starts_with("resolved_"))
}

coverage <- map(set_names(strata), build_coverage_table)

# ------------------------------------------------------------------------------
# 3. Species lists
# ------------------------------------------------------------------------------
save_species_list <- function(cov, keep, path) {
  cov %>% filter(keep) %>% select(species) %>% write_csv(path)
}

iwalk(coverage, function(cov, s) {
  write_csv(cov, this_module_dir("data_clean", paste0(s, "_trait_climatic_coverage.csv")))
  
  save_species_list(cov, cov$full_coverage_hmsc,
                    this_module_dir("data_clean", paste0(s, "_species_full_coverage_hmsc.csv")))
  save_species_list(cov, cov$climatic_ok,
                    this_module_dir("data_clean", paste0(s, "_species_full_climatic.csv")))
  
  out <- function(pool, cov_type) this_module_dir("outputs",
                                                  paste0(s, "_species_pool_", pool, "_", cov_type, ".csv"))
  save_species_list(cov, cov$climatic_ok,                           out("broad",  "full_clim_cov"))
  save_species_list(cov, cov$in_strict_pool & cov$climatic_ok,      out("strict", "full_clim_cov"))
  save_species_list(cov, cov$full_coverage_hmsc,                    out("broad",  "full_func_clim_cov"))
  save_species_list(cov, cov$in_strict_pool & cov$full_coverage_hmsc, out("strict", "full_func_clim_cov"))
})

# ------------------------------------------------------------------------------
# 4. Summary
# ------------------------------------------------------------------------------
coverage_summary <- imap_dfr(coverage, ~ tibble(
  stratum                       = .y,
  functional_traits             = paste(functional_traits[[.y]], collapse = ", "),
  shared                        = nrow(.x),
  shared_climatic               = sum(.x$climatic_ok),        # thermophilization
  shared_functional             = sum(.x$functional_ok),
  shared_climatic_functional    = sum(.x$full_coverage_hmsc),
  frequent                      = sum(.x$in_strict_pool),
  frequent_climatic             = sum(.x$in_strict_pool & .x$climatic_ok),
  frequent_climatic_functional  = sum(.x$in_strict_pool & .x$full_coverage_hmsc)  # HMSC
))
write_csv(coverage_summary, this_module_dir("outputs", "species_pool_coverage_summary.csv"))
print(coverage_summary %>% select(-functional_traits))

iwalk(coverage, function(cov, s) {
  cat("\n--", s, ": shared-pool species without full coverage --\n")
  print(cov %>% filter(!full_coverage_hmsc) %>% select(species, in_strict_pool, missing_traits), n = Inf)
})
