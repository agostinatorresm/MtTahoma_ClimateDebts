# =============================================================================
# 04_regional_affinity_index_rates.R
#
# Pipeline: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/
# Author: Agostina Torres
#           (after 03_regional_affinity_sensitivity.R)
#
# Effect of the geographic extent of the occurrence records on the community
# climatic-affinity shifts. Plot-level CTI and CMI and their rates of change are
# computed with range-wide (North America) and regional (500 and 1000 km)
# point affinities, using the same weights and species pool as the main
# analysis (overstory: relative basal area; understory: cover; shared pool).
# Point affinities are used for all three sets, so the comparison isolates the
# effect of the affinities; the main analysis uses the median over 100 affinity
# resamples, whose mean rates are very close to those of point affinities.
#
# Outputs: 02_integrate_gbif_chelsa/outputs/regional_affinity_index_rates.csv
#          02_integrate_gbif_chelsa/outputs/regional_affinity_index_rates_plot_cor.csv
# =============================================================================

library(tidyverse)
library(here)

module_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
STRATA <- tibble(stratum = c("overstory", "understory"), weight = c("basal_area_spp_perc", "cover"))

# Affinity sets
aff_sets <- list(
  rangewide = read_csv(module_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE),
  regional_500km = read_csv(module_dir("outputs", "species_climatic_affinity_regional_500km.csv"), show_col_types = FALSE),
  regional_1000km = read_csv(module_dir("outputs", "species_climatic_affinity_regional_1000km.csv"), show_col_types = FALSE)) %>%
  imap_dfr(~ transmute(.x, species = str_replace_all(species, "_", " "), MAT_mean, MAP_mean, affinity_set = .y))

# Plot-level indices (shared pool, abundance-weighted means of point affinities)
plot_rates <- pmap_dfr(STRATA, function(stratum, weight) {
  d <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", paste0(stratum, "_long.csv")), show_col_types = FALSE) %>%
    mutate(species = str_replace_all(species, "_", " "))
  tcol <- intersect(c("period", "time", "survey", "census"), names(d))[1]
  shared <- d %>% filter(occurrence > 0) %>% group_by(species) %>% filter(n_distinct(.data[[tcol]]) == 2) %>% pull(species) %>% unique()
  d <- d %>% filter(species %in% shared, occurrence > 0, !is.na(.data[[weight]]), .data[[weight]] > 0) %>%
    mutate(period = .data[[tcol]], w = .data[[weight]])
  d %>% inner_join(aff_sets, by = "species", relationship = "many-to-many") %>%
    group_by(affinity_set, plot_id, period) %>%
    summarise(CTI = weighted.mean(MAT_mean, w), CMI = weighted.mean(MAP_mean, w), year = mean(year), .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = c(CTI, CMI, year)) %>%
    { nm <- names(.); past <- grep("past", nm, value = TRUE); pres <- grep("present", nm, value = TRUE)
    mutate(., gap = (.data[[pres[grepl("year", pres)]]] - .data[[past[grepl("year", past)]]]) / 10,
           rate_CTI = (.data[[pres[grepl("CTI", pres)]]] - .data[[past[grepl("CTI", past)]]]) / gap,
           rate_CMI = (.data[[pres[grepl("CMI", pres)]]] - .data[[past[grepl("CMI", past)]]]) / gap) } %>%
    filter(is.finite(rate_CTI), is.finite(rate_CMI)) %>%
    mutate(stratum = stratum) %>%
    select(stratum, affinity_set, plot_id, rate_CTI, rate_CMI)
})

# Mean rate across plots with 95% confidence interval
mean_ci <- function(x) { t <- t.test(x); tibble(n = length(x), mean = mean(x), lower_95 = t$conf.int[1], upper_95 = t$conf.int[2]) }
summ <- plot_rates %>%
  pivot_longer(c(rate_CTI, rate_CMI), names_to = "index", values_to = "rate") %>%
  group_by(stratum, index, affinity_set) %>%
  group_modify(~ mean_ci(.x$rate)) %>% ungroup() %>%
  mutate(across(c(mean, lower_95, upper_95), ~ signif(.x, 3)))
write_csv(summ, module_dir("outputs", "regional_affinity_index_rates.csv"))
cat("\nMean rates of CTI (°C/decade) and CMI (mm/decade) with range-wide and regional affinities:\n")
print(summ, n = Inf)

# Agreement of plot-level rates between affinity sets
pc <- plot_rates %>%
  pivot_longer(c(rate_CTI, rate_CMI), names_to = "index", values_to = "rate") %>%
  pivot_wider(names_from = affinity_set, values_from = rate) %>%
  group_by(stratum, index) %>%
  summarise(r_500km = cor(rangewide, regional_500km, use = "complete.obs"),
            r_1000km = cor(rangewide, regional_1000km, use = "complete.obs"), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
write_csv(pc, module_dir("outputs", "regional_affinity_index_rates_plot_cor.csv"))
cat("\nCorrelation of plot-level rates, range-wide vs regional affinities:\n")
print(pc)



# -----------------------------------------------------------------------------
# 5. Species contributions to the overstory CMI rate with each affinity set
# -----------------------------------------------------------------------------
# With relative weights p (summing to 1 within a plot), ΔCMI = Σ Δp_i (a_i − c)
# for any constant c; c = past plot CMI, so a species contributes positively when
# it declines and is drier than the past community, or increases and is wetter.
aff_pt <- list(
  rangewide       = read_csv(module_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE),
  regional_1000km = read_csv(module_dir("outputs", "species_climatic_affinity_regional_1000km.csv"), show_col_types = FALSE),
  regional_500km  = read_csv(module_dir("outputs", "species_climatic_affinity_regional_500km.csv"), show_col_types = FALSE)) %>%
  imap_dfr(~ transmute(.x, species = str_replace_all(species, "_", " "), a = MAP_mean, affinity_set = .y))

ov <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", "overstory_long.csv"), show_col_types = FALSE) %>%
  mutate(species = str_replace_all(species, "_", " "))
tcol <- intersect(c("period", "time", "survey", "census"), names(ov))[1]
shared_ov <- ov %>% filter(occurrence > 0) %>% group_by(species) %>% filter(n_distinct(.data[[tcol]]) == 2) %>% pull(species) %>% unique()
ov <- ov %>% filter(species %in% shared_ov, occurrence > 0, basal_area_spp_perc > 0) %>%
  transmute(plot_id, period = .data[[tcol]], year, species, w = basal_area_spp_perc)
gap <- ov %>% group_by(plot_id, period) %>% summarise(year = mean(year), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = year) %>%
  transmute(plot_id, gap = (present - past) / 10)

contrib <- ov %>% select(-year) %>%
  inner_join(aff_pt, by = "species", relationship = "many-to-many") %>%
  group_by(affinity_set, plot_id, period) %>% mutate(p = w / sum(w)) %>% ungroup() %>%
  select(affinity_set, plot_id, species, a, period, p) %>%
  pivot_wider(names_from = period, values_from = p, values_fill = 0) %>%
  group_by(affinity_set, plot_id) %>% mutate(c_past = sum(past * a)) %>% ungroup() %>%
  left_join(gap, by = "plot_id") %>%
  mutate(contrib = (present - past) * (a - c_past) / gap) %>%
  group_by(affinity_set, species) %>%
  summarise(MAP_affinity = first(a), mean_contrib = mean(contrib), .groups = "drop") %>%
  pivot_wider(names_from = affinity_set, values_from = c(MAP_affinity, mean_contrib)) %>%
  arrange(desc(abs(mean_contrib_rangewide)))
write_csv(contrib, module_dir("outputs", "overstory_CMI_species_contributions_by_affinity_set.csv"))
cat("\nSpecies contributions to the mean overstory CMI rate (mm/decade), point affinities:\n")
print(mutate(contrib, across(where(is.double), ~ round(.x, 1))), n = 12, width = Inf)
