# ------------------------------------------------------------------------------
# 01_define_species_pool.R
# Module: pipeline/01_data_prep/02_species_pool
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Defines the species pools used by all downstream modules (trait table,
# climatic affinity, phylogeny, thermophilization, HMSC), so that every module
# works on the same species lists.
#
# Starting set: all species in overstory_long.csv / understory_long.csv
# (species after cleaning, in plots surveyed in both periods).
#
# Pools (nested):
#   shared pool   (code: broad)  -- species recorded in both surveys. A species
#                                   recorded in one survey only cannot be
#                                   compared between surveys.
#   frequent pool (code: strict) -- shared species also recorded in more than
#                                   PLOT_THRESHOLD plots in each survey. Species
#                                   in fewer plots give too little information
#                                   to estimate their responses and residual
#                                   associations in HMSC.
#
# Coverage of climatic affinity and functional traits is added in
# 02_full_coverage_species.R; the analysis sets are
#   thermophilization: shared pool with climatic affinity
#   HMSC:              frequent pool with climatic affinity and functional traits
#
# Inputs (01_cleaning/outputs/): overstory_long.csv, understory_long.csv
# Outputs (data_clean/):
#   species_pool_candidates_<stratum>.csv  every species after cleaning, with
#       n_plots_past, n_plots_present, present_both_periods,
#       meets_plot_threshold, included_broad, included_strict,
#       exclusion_reason_strict
#   <stratum>_final_list_broad.csv / _strict.csv  species lists per pool
#   species_pool_broad.csv / species_pool_strict.csv  species + layer,
#       both layers combined
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)
cleaning_dir    <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)

dir.create(this_module_dir("data_clean"), recursive = TRUE, showWarnings = FALSE)

strata         <- c("overstory", "understory")
PLOT_THRESHOLD <- 5   # frequent pool: > 5 plots in each survey (both strata)

# ------------------------------------------------------------------------------
# 1. Candidate table: plots occupied per survey and pool membership
# ------------------------------------------------------------------------------
build_species_pool_candidates <- function(long_df, plot_threshold) {
  per_period <- long_df %>%
    filter(occurrence == 1) %>%
    distinct(species, period, plot_id) %>%
    count(species, period, name = "n_plots") %>%
    pivot_wider(names_from = period, values_from = n_plots, values_fill = 0,
                names_glue = "n_plots_{period}")
  for (col in c("n_plots_past", "n_plots_present")) {
    if (!col %in% names(per_period)) per_period[[col]] <- 0L
  }
  
  per_period %>%
    mutate(
      present_both_periods = n_plots_past > 0 & n_plots_present > 0,
      meets_plot_threshold = n_plots_past > plot_threshold & n_plots_present > plot_threshold,
      included_broad       = present_both_periods,
      included_strict      = present_both_periods & meets_plot_threshold,
      exclusion_reason_strict = case_when(
        included_strict       ~ NA_character_,
        !present_both_periods ~ "not present in both periods",
        !meets_plot_threshold ~ paste0("<= ", plot_threshold, " plots in at least one period")
      )
    ) %>%
    select(species, n_plots_past, n_plots_present, present_both_periods,
           meets_plot_threshold, included_broad, included_strict, exclusion_reason_strict) %>%
    arrange(desc(included_strict), desc(included_broad), desc(n_plots_past + n_plots_present))
}

candidates <- map(set_names(strata), function(s) {
  long_df <- read_csv(cleaning_dir("outputs", paste0(s, "_long.csv")), show_col_types = FALSE)
  cand <- build_species_pool_candidates(long_df, PLOT_THRESHOLD)
  
  write_csv(cand, this_module_dir("data_clean", paste0("species_pool_candidates_", s, ".csv")))
  write_csv(cand %>% filter(included_broad)  %>% select(species),
            this_module_dir("data_clean", paste0(s, "_final_list_broad.csv")))
  write_csv(cand %>% filter(included_strict) %>% select(species),
            this_module_dir("data_clean", paste0(s, "_final_list_strict.csv")))
  cand
})

# ------------------------------------------------------------------------------
# 2. Combined pools (species + stratum)
# ------------------------------------------------------------------------------
pool_table <- function(flag) {
  imap_dfr(candidates, ~ .x %>% filter(.data[[flag]]) %>% transmute(species, stratum = .y))
}
write_csv(pool_table("included_broad"),  this_module_dir("data_clean", "species_pool_broad.csv"))
write_csv(pool_table("included_strict"), this_module_dir("data_clean", "species_pool_strict.csv"))

# ------------------------------------------------------------------------------
# 3. Summary
# ------------------------------------------------------------------------------
pool_summary <- imap_dfr(candidates, ~ tibble(
  stratum                 = .y,
  after_cleaning          = nrow(.x),
  shared                  = sum(.x$included_broad),
  not_in_both_surveys     = sum(!.x$present_both_periods),
  frequent                = sum(.x$included_strict),
  shared_below_threshold  = sum(.x$included_broad & !.x$included_strict)
))
print(pool_summary)
iwalk(candidates, ~ { cat("\n--", .y, "shared pool --\n"); print(filter(.x, included_broad), n = Inf) })
