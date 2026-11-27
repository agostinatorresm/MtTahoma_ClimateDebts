# =============================================================================
# 09_rarefied_thermophilization.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; after 05)
# Author: Agostina Torres
#
# Understory thermophilization recomputed with the past data reduced to the
# present sampling effort (Wilson et al. 2024; procedure as in 05, applied to
# the species set of the thermophilization analyses).
#
#   1. Species: shared pool with climatic affinity (as in 01_compute_cti_cmi
#      of the thermophilization module). Richness per plot over these species
#      in each survey; the difference in mean richness (d, rounded) is removed
#      from every past plot (plots with d species or fewer unchanged).
#   2. Retained species drawn without replacement with probability
#      proportional to past cover (cover-weighted removal; missing cover set to
#      0.1%), and with equal probability (uniform removal) as a comparison;
#      n_rare repetitions each. Present data unchanged.
#   3. For each reduced dataset, the past community indices are recomputed
#      with the retained species, for both weights of the thermophilization
#      analyses: cover-weighted (species with cover) and occurrence-based
#      (unweighted mean over the species present, including those recorded
#      only in the 50-m2 census). Rates = (present - past) / survey interval
#      in decades, per plot.
# Indices use point affinities (mean of each species' truncated
# distribution); the observed indices are recomputed in the same way, so that
# observed and reduced rates are directly comparable. The observed point-
# affinity rates are checked against the median rates of 01.
#
# Outputs per index and weighting: mean rate across plots for each reduction
# (distribution over reductions) and the per-plot mean rate over reductions
# (for the elevational pattern).
#
# Inputs : pipeline/01_data_prep/01_cleaning/outputs/understory_long.csv, plot_info.csv,
#          present_census_year.csv
#          pipeline/01_data_prep/02_species_pool/outputs/understory_species_pool_broad_full_clim_cov.csv
#          pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/data_clean/species_climatic_affinity.csv
#          pipeline/02_analysis/01_thermophilization/02_understory/outputs/understory_thermo_xero_rates.csv (check)
# Outputs: pipeline/04_robustness/outputs/rarefaction_thermo/
#            rarefaction_thermo_richness.csv     mean richness per survey, d
#            rarefaction_thermo_mean_rates.csv   mean rate per reduction (index x weight x weighting x rep)
#            rarefaction_thermo_summary.csv      observed vs reduced mean rates (mean, 2.5-97.5%)
#            rarefaction_thermo_plot_rates.csv   per-plot observed rate and mean reduced rate
#            rarefaction_thermo_plot_rates_by_rep.csv  per-plot rate of every reduction and of the
#                                                observed data (input of 12_rarefied_rate_models_bayes.R)
#          pipeline/04_robustness/figures/figSR5_rarefaction_thermophilization.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

n_rare      <- 100
proxy_cover <- 0.1
weightings  <- c("cover", "uniform")
set.seed(20261012)

cleaning_dir <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
out_dir <- here("pipeline", "04_robustness", "outputs", "rarefaction_thermo")
sup_dir <- here("pipeline", "04_robustness", "figures")
walk(c(out_dir, sup_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(),
          axis.line  = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 11),
          axis.title.y = element_text(colour = "black", size = 11),
          legend.text  = element_text(colour = "black"),
          strip.text   = element_text(colour = "black", face = "bold", size = 10),
          plot.title   = element_text(colour = "black", size = 11, face = "bold"),
          plot.margin  = margin(12, 10, 10, 10))
}
col_data <- c(observed = "grey40", cover = "#E07B39", uniform = "#3A7CA5")
lab_data <- c(observed = "Observed", cover = "Reduced, cover-weighted", uniform = "Reduced, uniform")
lab_weight <- c(cover = "Cover-weighted index", occurrence = "Occurrence-based index")

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
pool <- read_csv(here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                      "understory_species_pool_broad_full_clim_cov.csv"), show_col_types = FALSE) %>%
  pull(species) %>% unique()
affinity <- read_csv(here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa",
                          "data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)

veg <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE) %>%
  filter(occurrence == 1, species %in% pool, species %in% affinity$species) %>%
  distinct(period, plot_id, species, .keep_all = TRUE) %>%
  select(period, plot_id, species, cover) %>%
  left_join(affinity, by = "species")

plot_info <- read_csv(cleaning_dir("outputs", "plot_info.csv"), show_col_types = FALSE)
plot_fixed <- plot_info %>% filter(period == "past") %>% distinct(plot_id, .keep_all = TRUE) %>%
  select(plot_id, elevation, year_past = year)
present_year <- read_csv(cleaning_dir("outputs", "present_census_year.csv"), show_col_types = FALSE) %>%
  filter(stratum == "understory") %>% select(plot_id, year_present = year)
year_gap <- plot_fixed %>% inner_join(present_year, by = "plot_id") %>%
  mutate(year_gap_decades = (year_present - year_past) / 10)

plots <- intersect(unique(veg$plot_id[veg$period == "past"]), unique(veg$plot_id[veg$period == "present"]))
plots <- intersect(plots, year_gap$plot_id)
veg   <- veg %>% filter(plot_id %in% plots)
cat("Species in the shared pool with affinity:", n_distinct(veg$species), "| plots:", length(plots), "\n")

# -----------------------------------------------------------------------------
# 1. Richness and number of species to remove
# -----------------------------------------------------------------------------
richness <- veg %>% count(period, plot_id, name = "richness")
rich_summary <- richness %>% group_by(period) %>%
  summarise(mean_richness = mean(richness), .groups = "drop")
d <- round(rich_summary$mean_richness[rich_summary$period == "past"] -
             rich_summary$mean_richness[rich_summary$period == "present"])
rich_summary <- rich_summary %>% mutate(n_removed_per_plot = d)
print(rich_summary)
cat("Species removed per past plot:", d, "\n")
write_csv(rich_summary, file.path(out_dir, "rarefaction_thermo_richness.csv"))

# -----------------------------------------------------------------------------
# 2. Indices
# -----------------------------------------------------------------------------
# Per plot: cover-weighted (species with cover) and occurrence-based indices
index_by_plot <- function(df) {
  df %>%
    group_by(plot_id) %>%
    summarise(CTI_cover      = if (any(cover > 0, na.rm = TRUE)) weighted.mean(MAT_mean[cover > 0 & !is.na(cover)], cover[cover > 0 & !is.na(cover)]) else NA_real_,
              CMI_cover      = if (any(cover > 0, na.rm = TRUE)) weighted.mean(MAP_mean[cover > 0 & !is.na(cover)], cover[cover > 0 & !is.na(cover)]) else NA_real_,
              CTI_occurrence = mean(MAT_mean),
              CMI_occurrence = mean(MAP_mean),
              .groups = "drop")
}
to_long <- function(idx) {
  idx %>% pivot_longer(-plot_id, names_to = c("index", "weight"), names_sep = "_", values_to = "value")
}
present_idx <- to_long(index_by_plot(filter(veg, period == "present"))) %>% rename(present = value)
past_obs    <- to_long(index_by_plot(filter(veg, period == "past")))    %>% rename(past = value)

rates_from_past <- function(past_long) {
  past_long %>%
    inner_join(present_idx, by = c("plot_id", "index", "weight")) %>%
    inner_join(year_gap %>% select(plot_id, year_gap_decades), by = "plot_id") %>%
    mutate(rate = (present - past) / year_gap_decades) %>%
    filter(!is.na(rate))
}
obs_rates <- rates_from_past(past_obs)

# Check against the median rates of 01 (affinity resamples)
rates_01_file <- here("pipeline", "02_analysis", "01_thermophilization", "02_understory", "outputs",
                      "understory_thermo_xero_rates.csv")
if (file.exists(rates_01_file)) {
  r01 <- read_csv(rates_01_file, show_col_types = FALSE) %>%
    select(plot_id, weight = weight_var, CTI = rate_CTI, CMI = rate_CMI) %>%
    pivot_longer(c(CTI, CMI), names_to = "index", values_to = "rate_01")
  chk <- obs_rates %>% inner_join(r01, by = c("plot_id", "index", "weight")) %>%
    group_by(index, weight) %>%
    summarise(n = n(), mean_rate_point = mean(rate), mean_rate_01 = mean(rate_01),
              r = cor(rate, rate_01), .groups = "drop")
  cat("\nObserved rates, point affinities vs medians of 01:\n"); print(chk)
}

# -----------------------------------------------------------------------------
# 3. Reductions
# -----------------------------------------------------------------------------
past_split <- veg %>% filter(period == "past") %>%
  mutate(w = if_else(is.na(cover) | cover <= 0, proxy_cover, cover)) %>%
  split(.$plot_id)

reduce_once <- function(weighting) {
  map_dfr(past_split, function(x) {
    n <- nrow(x)
    if (n <= d) return(x)
    keep <- sample(seq_len(n), n - d, prob = if (weighting == "cover") x$w else NULL)
    x[keep, ]
  })
}

plot_rates_red <- list()
mean_rates <- map_dfr(weightings, function(wt) {
  map_dfr(seq_len(n_rare), function(k) {
    r <- rates_from_past(to_long(index_by_plot(reduce_once(wt))) %>% rename(past = value))
    plot_rates_red[[paste(wt, k)]] <<- r %>% mutate(weighting = wt, rep = k)
    r %>% group_by(index, weight) %>%
      summarise(mean_rate = mean(rate), n_plots = n(), .groups = "drop") %>%
      mutate(weighting = wt, rep = k)
  })
})
plot_rates_red <- bind_rows(plot_rates_red)
write_csv(mean_rates, file.path(out_dir, "rarefaction_thermo_mean_rates.csv"))
# Per-plot rates of every reduction and of the observed data (point
# affinities), for the Bayesian rate models pooled across reductions (10)
write_csv(bind_rows(plot_rates_red %>% select(plot_id, index, weight, weighting, rep, rate),
                    obs_rates %>% transmute(plot_id, index, weight, weighting = "observed", rep = 0L, rate)),
          file.path(out_dir, "rarefaction_thermo_plot_rates_by_rep.csv"))

obs_mean <- obs_rates %>% group_by(index, weight) %>%
  summarise(mean_rate = mean(rate), n_plots = n(), .groups = "drop")
summary_tbl <- bind_rows(
  obs_mean %>% mutate(data = "observed", q025 = NA_real_, q975 = NA_real_),
  mean_rates %>% group_by(index, weight, data = weighting) %>%
    summarise(q025 = quantile(mean_rate, 0.025), q975 = quantile(mean_rate, 0.975),
              mean_rate = mean(mean_rate), n_plots = round(mean(n_plots)), .groups = "drop")
) %>%
  arrange(index, weight, data)
cat("\nMean rate across plots (CTI: degC per decade; CMI: mm per decade):\n")
print(summary_tbl %>% mutate(across(c(mean_rate, q025, q975), ~ signif(.x, 3))), n = Inf)
write_csv(summary_tbl, file.path(out_dir, "rarefaction_thermo_summary.csv"))

plot_tbl <- plot_rates_red %>%
  group_by(plot_id, index, weight, weighting) %>%
  summarise(rate_reduced = mean(rate), .groups = "drop") %>%
  pivot_wider(names_from = weighting, values_from = rate_reduced, names_prefix = "rate_reduced_") %>%
  left_join(obs_rates %>% select(plot_id, index, weight, rate_observed = rate), by = c("plot_id", "index", "weight")) %>%
  left_join(plot_fixed %>% select(plot_id, elevation), by = "plot_id")
write_csv(plot_tbl, file.path(out_dir, "rarefaction_thermo_plot_rates.csv"))
cat("\nCorrelation across plots, observed vs mean reduced (cover-weighted removal):\n")
print(plot_tbl %>% group_by(index, weight) %>%
        summarise(r = round(cor(rate_observed, rate_reduced_cover, use = "complete.obs"), 3), .groups = "drop"))

# -----------------------------------------------------------------------------
# Diagnostic figure (final version: Figure S4.1, 12_rarefied_rate_models_bayes.R): (A) mean rate observed vs reduced; (B) rate along elevation
# -----------------------------------------------------------------------------
y_lab <- c(CTI = "Thermophilization rate (°C per decade)", CMI = "Moisture-affinity shift (mm per decade)")
panel_mean <- function(ix) {
  summary_tbl %>% filter(index == ix) %>%
    mutate(data = factor(data, levels = names(lab_data)), weight = lab_weight[weight]) %>%
    ggplot(aes(data, mean_rate, colour = data)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 2.6) +
    facet_wrap(~ weight) +
    scale_colour_manual(values = col_data, labels = lab_data, guide = "none") +
    scale_x_discrete(labels = c(observed = "Observed", cover = "Cover-\nweighted", uniform = "Uniform")) +
    labs(x = NULL, y = y_lab[[ix]], title = paste0(ix, ": mean rate across plots")) +
    theme_mora()
}
panel_elev <- function(ix) {
  plot_tbl %>% filter(index == ix) %>%
    select(plot_id, weight, elevation, observed = rate_observed, cover = rate_reduced_cover) %>%
    pivot_longer(c(observed, cover), names_to = "data", values_to = "rate") %>%
    mutate(data = factor(data, levels = c("observed", "cover")), weight = lab_weight[weight]) %>%
    ggplot(aes(elevation, rate, colour = data, fill = data)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "loess", span = 0.75, formula = y ~ x, alpha = 0.15, linewidth = 0.9) +
    facet_wrap(~ weight) +
    scale_colour_manual(values = col_data, labels = lab_data, name = NULL) +
    scale_fill_manual(values = col_data, labels = lab_data, name = NULL) +
    labs(x = "Elevation (m)", y = y_lab[[ix]], title = paste0(ix, ": rate along elevation")) +
    theme_mora() + theme(legend.position = "bottom")
}
fig <- (panel_mean("CTI") | panel_elev("CTI")) / (panel_mean("CMI") | panel_elev("CMI")) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(sup_dir, "figSR5_rarefaction_thermophilization.png"), fig,
       width = 14, height = 10, dpi = 300, bg = "white")
cat("\nFigure written to", sup_dir, "\n")
