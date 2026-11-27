# ------------------------------------------------------------------------------
# 06_community_affinity_edges.R
# Module: pipeline/02_analysis/01_thermophilization/01_overstory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Cold and warm tails of the overstory community, to test whether the change
# in CTI comes from the loss of cold-affinity species (cold tail rises) or
# from the gain of warm-affinity species (warm tail rises); the same for the
# dry and wet tails of the moisture affinity (CMI).
#
# For each plot and survey, the abundance-weighted 10th, 50th and 90th
# percentiles of the affinities of the species present are computed (point
# affinities: mean of each species' truncated distribution). Weights:
# basal area as % of plot area. Weighted percentiles are interpolated between the
# cumulative-weight midpoints of the species ordered by affinity; with a
# single species all three equal its affinity.
#   q10: cold (dry) tail;  q50: weighted median;  q90: warm (wet) tail.
# Rates are (present - past) / survey interval in decades, with the plots and
# intervals of 01 (plots with the index in both surveys).
#
# Interpretation:
#   cold tail rises more than the warm tail  -> loss / decline of cold-affinity
#                                               species (retraction);
#   warm tail rises more than the cold tail  -> gain / increase of
#                                               warm-affinity species (expansion).
# The difference between tail rates is tested with a paired t-test across
# plots.
#
# Inputs : 01_cleaning/outputs/overstory_long.csv
#          02_species_pool/outputs/overstory_species_pool_broad_full_clim_cov.csv
#          05_species_distribution/02_integrate_gbif_chelsa/data_clean/species_climatic_affinity.csv
#          outputs/overstory_thermo_xero_rates.csv (01)
# Outputs: data_clean/overstory_affinity_edges.csv      edges per plot x survey
#          outputs/overstory_affinity_edge_rates.csv     rates per plot
#          outputs/overstory_affinity_edge_summary.csv   mean rates, CI, paired tests
#          figures/overstory_affinity_edges.png
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(hrbrthemes)
library(ggpubr)

this_module_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)
cleaning_dir     <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
affinity_dir     <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
species_pool_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)
for (d in c("data_clean", "outputs", "figures")) dir.create(this_module_dir(d), recursive = TRUE, showWarnings = FALSE)

stratum        <- "overstory"
PRIMARY_WEIGHT <- "basal_area_spp_perc"
EDGE_PROBS     <- c(q10 = 0.10, q50 = 0.50, q90 = 0.90)
LOESS_SPAN     <- 0.75
col_temp <- "coral"
col_prec <- "steelblue"

# 1. Inputs ---------------------------------------------------------------------
veg_long <- read_csv(cleaning_dir("outputs", paste0(stratum, "_long.csv")), show_col_types = FALSE)
species_pool <- read_csv(species_pool_dir("outputs", paste0(stratum, "_species_pool_broad_full_clim_cov.csv")),
                         show_col_types = FALSE) %>% pull(species)
affinity <- read_csv(affinity_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)
rates_01 <- read_csv(this_module_dir("outputs", paste0(stratum, "_thermo_xero_rates.csv")), show_col_types = FALSE) %>%
  filter(weight_var == PRIMARY_WEIGHT) %>%
  select(plot_id, elevation, year_gap_decades, rate_CTI, rate_CMI)

# 2. Weighted percentiles per plot x survey ----------------------------------------
weighted_quantile <- function(x, w, p) {
  o  <- order(x); x <- x[o]; w <- w[o] / sum(w)
  if (length(x) == 1) return(rep(x, length(p)))
  cw <- cumsum(w) - w / 2
  approx(cw, x, xout = p, rule = 2, ties = mean)$y
}

edges <- veg_long %>%
  filter(occurrence == 1, species %in% species_pool) %>%
  select(plot_id, period, species, weight = all_of(PRIMARY_WEIGHT)) %>%
  filter(!is.na(weight), weight > 0) %>%
  inner_join(affinity, by = "species") %>%
  group_by(plot_id, period) %>%
  reframe(edge      = names(EDGE_PROBS),
          MAT       = weighted_quantile(MAT_mean, weight, EDGE_PROBS),
          MAP       = weighted_quantile(MAP_mean, weight, EDGE_PROBS),
          n_species = n_distinct(species))
write_csv(edges, this_module_dir("data_clean", paste0(stratum, "_affinity_edges.csv")))

# 3. Rates per plot ---------------------------------------------------------------
edge_rates <- edges %>%
  pivot_longer(c(MAT, MAP), names_to = "affinity", values_to = "value") %>%
  select(-n_species) %>%
  pivot_wider(names_from = period, values_from = value) %>%
  inner_join(rates_01 %>% select(plot_id, elevation, year_gap_decades), by = "plot_id") %>%
  filter(!is.na(past), !is.na(present)) %>%
  mutate(rate  = (present - past) / year_gap_decades,
         index = if_else(affinity == "MAT", "CTI", "CMI"),
         edge  = factor(edge, levels = names(EDGE_PROBS)))
write_csv(edge_rates, this_module_dir("outputs", paste0(stratum, "_affinity_edge_rates.csv")))

cat("Plots with edges in both surveys:", n_distinct(edge_rates$plot_id), "\n")
cat("Species per plot and survey (median, range):",
    median(edges$n_species), paste(range(edges$n_species), collapse = "-"), "\n")
cat("Plots with a single species in a survey (edges equal):",
    edges %>% filter(n_species == 1) %>% distinct(plot_id) %>% nrow(), "\n")

# 4. Mean rates and paired comparison of edges -----------------------------------------
mean_ci <- function(x) {
  x <- x[!is.na(x)]; tt <- t.test(x)
  tibble(n = length(x), mean = mean(x), ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2], p_value = tt$p.value)
}
edge_summary <- edge_rates %>%
  group_by(index, edge) %>%
  summarise(mean_ci(rate), .groups = "drop")

paired_tests <- edge_rates %>%
  select(plot_id, index, edge, rate) %>%
  pivot_wider(names_from = edge, values_from = rate) %>%
  group_by(index) %>%
  summarise(mean_q10_minus_q90 = mean(q10 - q90, na.rm = TRUE),
            p_q10_vs_q90       = t.test(q10, q90, paired = TRUE)$p.value,
            .groups = "drop")

write_csv(left_join(edge_summary, paired_tests, by = "index"),
          this_module_dir("outputs", paste0(stratum, "_affinity_edge_summary.csv")))
cat("\nMean rate per edge (CTI: degC per decade; CMI: mm per decade):\n")
print(edge_summary, width = Inf)
cat("\nCold/dry edge minus warm/wet edge (paired across plots):\n")
print(paired_tests, width = Inf)

# 5. Figure: mean rates per edge and rates along elevation ------------------------------
if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = 0.3),
          axis.line.y = element_line(colour = "black", linewidth = 0.3),
          plot.margin = margin(10, 10, 10, 10))
}
edge_labels <- list(CTI = c(q10 = "Cold edge (q10)", q50 = "Median (q50)", q90 = "Warm edge (q90)"),
                    CMI = c(q10 = "Dry edge (q10)",  q50 = "Median (q50)", q90 = "Wet edge (q90)"))
y_labels <- c(CTI = "Rate (°C per decade)", CMI = "Rate (mm per decade)")
cols     <- c(CTI = col_temp, CMI = col_prec)

mean_panel <- function(ix) {
  edge_summary %>%
    filter(index == ix) %>%
    mutate(edge_lab = factor(edge_labels[[ix]][as.character(edge)], levels = edge_labels[[ix]])) %>%
    ggplot(aes(x = edge_lab, y = mean, ymin = ci_lo, ymax = ci_hi)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(colour = cols[[ix]], size = 0.5) +
    labs(x = NULL, y = y_labels[[ix]], title = ix) +
    theme_mora()
}
elev_panel <- function(ix) {
  edge_rates %>%
    filter(index == ix) %>%
    mutate(edge_lab = factor(edge_labels[[ix]][as.character(edge)], levels = edge_labels[[ix]])) %>%
    ggplot(aes(x = elevation, y = rate, linetype = edge_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "loess", span = LOESS_SPAN, se = TRUE, colour = cols[[ix]],
                fill = cols[[ix]], alpha = 0.12, linewidth = 0.8, formula = y ~ x) +
    scale_linetype_manual(values = c("dotted", "solid", "longdash"), name = NULL) +
    labs(x = "Elevation (m)", y = y_labels[[ix]]) +
    theme_mora() + theme(legend.position = "top")
}
fig <- ggarrange(mean_panel("CTI"), elev_panel("CTI"), mean_panel("CMI"), elev_panel("CMI"),
                 ncol = 2, nrow = 2, widths = c(1, 1.6), labels = c("A", "B", "C", "D"))
ggsave(this_module_dir("figures", paste0(stratum, "_affinity_edges.png")), fig,
       width = 9, height = 7, dpi = 300, bg = "white")
cat("\nSaved: figures/", stratum, "_affinity_edges.png\n", sep = "")
