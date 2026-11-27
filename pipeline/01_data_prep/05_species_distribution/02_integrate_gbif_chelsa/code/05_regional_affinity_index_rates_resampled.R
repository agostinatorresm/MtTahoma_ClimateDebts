# =============================================================================
# 05_regional_affinity_index_rates_resampled.R
#
# Pipeline: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/
# Author: Agostina Torres
#           (after 03_regional_affinity_sensitivity.R)
#
# Community climatic-affinity shifts computed with range-wide (North America)
# and regional (within 500 and 1000 km of Mount Rainier) affinities, following
# the procedure of the main analysis: for each species, 100 values of MAT and of
# MAP are sampled with replacement from its truncated (5th-95th percentile)
# distribution of occurrence climates; plot-level CTI and CMI are computed for
# each resample (abundance-weighted means; overstory: relative basal area;
# understory: cover; shared pool) and summarised as the median over the 100
# resamples; rates are the change between surveys divided by each plot's survey
# interval (per decade). Species with fewer than min_records regional records
# keep their range-wide distribution.
#
# Inputs : 02_integrate_gbif_chelsa/data_clean/gbif_occurrences_with_bioclim.csv
#          01_data_prep/01_cleaning/outputs/<stratum>_long.csv
# Outputs: 02_integrate_gbif_chelsa/outputs/regional_affinity_index_rates_resampled.csv
#          02_integrate_gbif_chelsa/outputs/regional_affinity_index_rates_resampled_plot_cor.csv
#          02_integrate_gbif_chelsa/figures/regional_affinity_index_rates.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)

module_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
dir.create(module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

park_lat <- 46.853; park_lon <- -121.760
radii_km <- c(1000, 500)
min_records <- 30
n_res <- 100
trunc_q <- c(0.05, 0.95)
STRATA <- tibble(stratum = c("overstory", "understory"), weight = c("basal_area_spp_perc", "cover"))
set.seed(2026)

set_levels <- c("rangewide" = "Range-wide", "regional_1000km" = "1000 km", "regional_500km" = "500 km")
col_index <- c(CTI = "#F26B5B", CMI = "#1A759F")
if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.text = element_text(colour = "black", size = 11),
          axis.title.x = element_text(colour = "black", size = 12, hjust = 0.5),
          axis.title.y = element_text(colour = "black", size = 12, hjust = 0.5),
          strip.text = element_text(colour = "black", size = 12, face = "bold"))
}

# -----------------------------------------------------------------------------
# 1. Occurrence records (same filters and thinning as the main affinities)
# -----------------------------------------------------------------------------
haversine <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  a <- sin((lat2 - lat1) * rad / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin((lon2 - lon1) * rad / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}
occ <- read_csv(module_dir("data_clean", "gbif_occurrences_with_bioclim.csv"),
                col_select = c(species, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, MAT, MAP),
                show_col_types = FALSE) %>%
  rename(lat = decimalLatitude, lon = decimalLongitude) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(MAT), !is.na(MAP),
         is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 1e5) %>%
  mutate(species = str_replace_all(species, "_", " "),
         cell_x = floor(lon * 120), cell_y = floor(lat * 120)) %>%
  distinct(species, cell_x, cell_y, .keep_all = TRUE) %>%
  mutate(dist_km = haversine(park_lat, park_lon, lat, lon))

# -----------------------------------------------------------------------------
# 2. Resampled affinities (100 per species and variable) for each set
# -----------------------------------------------------------------------------
# sample by index (sample() on a length-one vector would sample from 1:x)
draw <- function(x) { q <- quantile(x, trunc_q); xt <- x[x >= q[1] & x <= q[2]]; xt[sample.int(length(xt), n_res, replace = TRUE)] }
resample_set <- function(d) d %>% group_by(species) %>%
  summarise(n = n(), rep = list(seq_len(n_res)), MAT = list(draw(MAT)), MAP = list(draw(MAP)), .groups = "drop") %>%
  unnest(c(rep, MAT, MAP))

rw <- resample_set(occ) %>% mutate(affinity_set = "rangewide")
sets <- list(rw)
for (r in radii_km) {
  occ_r <- filter(occ, dist_km <= r) %>% group_by(species) %>% filter(n() >= min_records) %>% ungroup()
  reg <- resample_set(occ_r)
  fallback <- filter(rw, !species %in% reg$species)              # too few regional records
  sets[[length(sets) + 1]] <- bind_rows(reg, fallback) %>% mutate(affinity_set = paste0("regional_", r, "km"))
}
aff_res <- bind_rows(sets) %>% select(affinity_set, species, rep, MAT, MAP)

# -----------------------------------------------------------------------------
# 3. Plot-level indices (median over resamples) and rates
# -----------------------------------------------------------------------------
plot_rates <- pmap_dfr(STRATA, function(stratum, weight) {
  d <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", paste0(stratum, "_long.csv")), show_col_types = FALSE) %>%
    mutate(species = str_replace_all(species, "_", " "))
  tcol <- intersect(c("period", "time", "survey", "census"), names(d))[1]
  shared <- d %>% filter(occurrence > 0) %>% group_by(species) %>% filter(n_distinct(.data[[tcol]]) == 2) %>% pull(species) %>% unique()
  d <- d %>% filter(species %in% shared, occurrence > 0, !is.na(.data[[weight]]), .data[[weight]] > 0) %>%
    transmute(plot_id, period = .data[[tcol]], year, species, w = .data[[weight]])
  yrs <- d %>% group_by(plot_id, period) %>% summarise(year = mean(year), .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = year, names_prefix = "year_")
  idx <- d %>% select(-year) %>%
    inner_join(aff_res, by = "species", relationship = "many-to-many") %>%
    group_by(affinity_set, plot_id, period, rep) %>%
    summarise(CTI = weighted.mean(MAT, w), CMI = weighted.mean(MAP, w), .groups = "drop") %>%
    group_by(affinity_set, plot_id, period) %>%
    summarise(CTI = median(CTI), CMI = median(CMI), .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = c(CTI, CMI))
  past <- grep("past", names(yrs), value = TRUE); pres <- grep("present", names(yrs), value = TRUE)
  idx %>% left_join(yrs, by = "plot_id") %>%
    mutate(gap = (.data[[pres]] - .data[[past]]) / 10,
           rate_CTI = (.data[[sub("year", "CTI", pres)]] - .data[[sub("year", "CTI", past)]]) / gap,
           rate_CMI = (.data[[sub("year", "CMI", pres)]] - .data[[sub("year", "CMI", past)]]) / gap,
           stratum = stratum) %>%
    filter(is.finite(rate_CTI), is.finite(rate_CMI)) %>%
    select(stratum, affinity_set, plot_id, rate_CTI, rate_CMI)
})

mean_ci <- function(x) { t <- t.test(x); tibble(n = length(x), mean = mean(x), lower_95 = t$conf.int[1], upper_95 = t$conf.int[2]) }
long <- plot_rates %>% pivot_longer(c(rate_CTI, rate_CMI), names_to = "index", values_to = "rate") %>%
  mutate(index = str_remove(index, "rate_"))
summ <- long %>% group_by(stratum, index, affinity_set) %>%
  group_modify(~ mean_ci(.x$rate)) %>% ungroup() %>%
  mutate(across(c(mean, lower_95, upper_95), ~ signif(.x, 3)))
write_csv(summ, module_dir("outputs", "regional_affinity_index_rates_resampled.csv"))
cat("\nMean rates (CTI, °C/decade; CMI, mm/decade), median over", n_res, "affinity resamples:\n")
print(arrange(summ, stratum, index, match(affinity_set, names(set_levels))), n = Inf)

pc <- long %>% pivot_wider(names_from = affinity_set, values_from = rate) %>%
  group_by(stratum, index) %>%
  summarise(r_1000km = cor(rangewide, regional_1000km, use = "complete.obs"),
            r_500km = cor(rangewide, regional_500km, use = "complete.obs"), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
write_csv(pc, module_dir("outputs", "regional_affinity_index_rates_resampled_plot_cor.csv"))
cat("\nCorrelation of plot-level rates with those based on range-wide affinities:\n"); print(pc)

# -----------------------------------------------------------------------------
# 4. Figure: mean rate (95% CI) with each set of affinities
# -----------------------------------------------------------------------------
fig_d <- summ %>%
  mutate(set = factor(set_levels[affinity_set], levels = set_levels),
         index = factor(index, levels = c("CTI", "CMI")),
         stratum = str_to_title(stratum),
         panel = factor(paste(stratum, if_else(index == "CTI", "ΔCTI (°C decade⁻¹)", "ΔCMI (mm decade⁻¹)"), sep = " – "),
                        levels = c(paste("Overstory", c("ΔCTI (°C decade⁻¹)", "ΔCMI (mm decade⁻¹)"), sep = " – "),
                                   paste("Understory", c("ΔCTI (°C decade⁻¹)", "ΔCMI (mm decade⁻¹)"), sep = " – "))))
p <- ggplot(fig_d, aes(set, mean, colour = index)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), linewidth = 1.6, lineend = "round") +
  geom_point(size = 3) +
  scale_colour_manual(values = col_index, guide = "none") +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  labs(x = "Species affinities", y = "Mean rate across plots (95% CI)") +
  theme_mora()
fig_file <- module_dir("figures", "regional_affinity_index_rates.png")
ggsave(fig_file, p, width = 9, height = 7, dpi = 300, bg = "white")
cat("\nFigure written to", fig_file, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
