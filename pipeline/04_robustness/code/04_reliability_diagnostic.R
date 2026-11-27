# =============================================================================
# 04_reliability_diagnostic.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Combines the sampling-reliability diagnostics into one table per
# understory species of the evaluated set (00: recorded in more than five
# plots in each survey, as the frequent pool of the HMSC models) and writes the species lists for the HMSC robustness
# models. The flags are used to refit the models without flagged species and
# compare conclusions; they are not filters of the main analysis.
#
# Flags:
#   Abundance (cover when present)
#     patchy             : mean present-day subplot occupancy < 0.4 (01b)
#     cover_unreliable   : relative-cover change not distinguishable from
#                          present-day subplot sampling variation (fewer than
#                          half of the plot tests valid, or relative shrinkage
#                          >= 0.75; 02)
#   Occurrence
#     area_sensitive     : past occurrence increases with plot area
#                          (positive, FDR-significant area coefficient; 03)
#   Frequency (both strata)
#     rare               : recorded in fewer than min_records plots on average
#                          over the two surveys, or in five plots or fewer in
#                          either survey. The mean is used rather than the
#                          minimum so that common species that declined
#                          strongly are not flagged as rare. Rare species are
#                          kept in the models; the flag marks species whose
#                          changes are not interpreted individually (after
#                          Wilson et al. 2024, Ecology 105: e4242, who
#                          restricted interpretation to species with at least
#                          10 records in each survey)
#
# NA means that the species was not evaluated by that diagnostic (too few
# records), not that it passed it.
#
# Inputs : pipeline/01_data_prep/01_cleaning/outputs/<stratum>_long.csv
#          pipeline/01_data_prep/02_species_pool/outputs/<stratum>_species_pool_strict_full_func_clim_cov.csv
#            (species of the HMSC models)
#          pipeline/04_robustness/outputs/accumulation_curve/accumulation_curve_summary.csv (01b)
#          pipeline/04_robustness/outputs/bootstrapping/bootstrap_species.csv (02)
#          pipeline/04_robustness/outputs/occurrence/area_sensitivity_species.csv (03)
# Outputs: pipeline/04_robustness/outputs/
#            composite_reliability_diagnostic.csv
#            species_to_exclude_robustness_abundance.csv
#            species_to_exclude_robustness_occurrence.csv
#            rare_species_flag.csv (overstory and understory)
# =============================================================================

library(tidyverse)
library(here)

out_dir   <- here("pipeline", "04_robustness", "outputs")
input_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
min_records <- 10

# -----------------------------------------------------------------------------
# Rare species: number of plots with the species in each survey
# -----------------------------------------------------------------------------
rare_flag <- map_dfr(c("overstory", "understory"), function(st) {
  read_csv(file.path(input_dir, paste0(st, "_long.csv")), show_col_types = FALSE) %>%
    filter(occurrence > 0) %>%
    distinct(period, plot_id, species) %>%
    count(species, period) %>%
    pivot_wider(names_from = period, values_from = n, names_prefix = "n_plots_", values_fill = 0) %>%
    mutate(stratum = st, .before = 1)
}) %>%
  mutate(n_plots_mean = (n_plots_past + n_plots_present) / 2,
         rare = n_plots_mean < min_records | pmin(n_plots_past, n_plots_present) <= 5) %>%
  arrange(stratum, species)
write_csv(rare_flag, file.path(out_dir, "rare_species_flag.csv"))
cat("Rare species (mean of the two surveys below", min_records, "plots, or five plots or fewer in a survey):\n")
print(rare_flag %>% count(stratum, rare))

patchy <- read_csv(file.path(out_dir, "accumulation_curve", "accumulation_curve_summary.csv"), show_col_types = FALSE)
boot   <- read_csv(file.path(out_dir, "bootstrapping", "bootstrap_species.csv"), show_col_types = FALSE) %>%
  filter(response == "relative_cover") %>%
  transmute(species, boot_n_tests = n_tests, boot_pct_valid = pct_valid, boot_rel_shrinkage = rel_shrinkage,
            cover_unreliable = !bootstrap_reliable)
area   <- read_csv(file.path(out_dir, "occurrence", "area_sensitivity_species.csv"), show_col_types = FALSE) %>%
  transmute(species, area_coef, area_p_adj, area_sensitive)

# Species of the HMSC models (frequent pool with complete traits and affinity)
pool_file <- function(st) here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                               paste0(st, "_species_pool_strict_full_func_clim_cov.csv"))
hmsc_species <- map_dfr(c("overstory", "understory"), function(st) {
  if (!file.exists(pool_file(st))) return(NULL)
  read_csv(pool_file(st), show_col_types = FALSE) %>% distinct(species) %>% mutate(stratum = st)
}) %>%
  bind_rows(tibble(species = character(), stratum = character()))
rare_flag <- rare_flag %>%
  left_join(hmsc_species %>% mutate(in_hmsc = TRUE), by = c("stratum", "species")) %>%
  mutate(in_hmsc = coalesce(in_hmsc, FALSE))
write_csv(rare_flag, file.path(out_dir, "rare_species_flag.csv"))
cat("\nRare species among the HMSC species:\n")
print(rare_flag %>% filter(in_hmsc) %>% count(stratum, rare))
cat("\nRare HMSC species:\n")
print(rare_flag %>% filter(in_hmsc, rare) %>% select(stratum, species, n_plots_past, n_plots_present), n = Inf)

# Evaluated set: understory species in more than five plots in each survey
composite <- rare_flag %>%
  filter(stratum == "understory", n_plots_past > 5, n_plots_present > 5) %>%
  select(-stratum) %>%
  left_join(patchy, by = "species") %>%
  left_join(boot, by = "species") %>%
  left_join(area, by = "species") %>%
  mutate(concerns_abundance  = rowSums(cbind(patchy, cover_unreliable), na.rm = TRUE),
         concerns_occurrence = as.integer(coalesce(area_sensitive, FALSE))) %>%
  arrange(desc(concerns_abundance + concerns_occurrence), species)

write_csv(composite, file.path(out_dir, "composite_reliability_diagnostic.csv"))

cat("\nFlags among the understory HMSC species:\n")
print(composite %>% filter(in_hmsc) %>%
        summarise(n = n(), rare = sum(rare, na.rm = TRUE), patchy = sum(patchy, na.rm = TRUE),
                  cover_unreliable = sum(cover_unreliable, na.rm = TRUE),
                  area_sensitive = sum(area_sensitive, na.rm = TRUE)))

cat("Species evaluated:", nrow(composite), "\n")
cat("\nAbundance flags:\n"); print(count(composite, patchy, cover_unreliable))
cat("\nOccurrence flags:\n"); print(count(composite, area_sensitive))

write_csv(composite %>% filter(concerns_abundance >= 1) %>% select(species),
          file.path(out_dir, "species_to_exclude_robustness_abundance.csv"))
write_csv(composite %>% filter(concerns_occurrence >= 1) %>% select(species),
          file.path(out_dir, "species_to_exclude_robustness_occurrence.csv"))
cat("\nSpecies lists for the robustness models: abundance", sum(composite$concerns_abundance >= 1),
    "| occurrence", sum(composite$concerns_occurrence >= 1), "\n")
