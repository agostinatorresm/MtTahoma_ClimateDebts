# =============================================================================
# 03_regional_affinity_sensitivity.R
#
# Pipeline: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/
# Author: Agostina Torres
#           (after the script that writes gbif_occurrences_with_bioclim.csv and
#            species_climatic_affinity.csv)
#
# Sensitivity of species' climatic affinities to the geographic extent of the
# occurrence records. The main affinities use occurrences across North America,
# which describe each species' broad climatic tolerance but may not represent
# the populations of the study region (e.g. species with coastal and interior
# populations). Here affinities are recomputed from occurrences within a radius
# of Mount Rainier and compared with the range-wide affinities. A broad radius
# is used rather than the Cascades alone: a small region truncates niches to the
# regional climate and makes affinities reflect species' positions along the
# same elevational gradient used to quantify community shifts.
#
# Procedure (as for the main affinities):
#   - records with coordinate uncertainty >= 100 km removed
#   - one record per species per 30-arc-second (~1 km) grid cell
#   - point affinity = mean of MAT and of MAP within the 5th-95th percentiles
#     of each species' values
# The range-wide affinities are first recomputed with this procedure and checked
# against species_climatic_affinity.csv, so that range-wide and regional values
# are computed identically. Regional affinities require at least min_records
# thinned records within the radius; species below the threshold keep their
# range-wide affinity (flagged).
#
# Inputs : 02_integrate_gbif_chelsa/data_clean/gbif_occurrences_with_bioclim.csv
#          02_integrate_gbif_chelsa/data_clean/species_climatic_affinity.csv
#          01_data_prep/01_cleaning/outputs/<stratum>_long.csv (species of the shared pools)
# Outputs: 02_integrate_gbif_chelsa/outputs/species_climatic_affinity_regional_<radius>km.csv
#          02_integrate_gbif_chelsa/outputs/regional_affinity_comparison.csv
#          02_integrate_gbif_chelsa/figures/regional_vs_rangewide_affinity.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)

module_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
dir.create(module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

park_lat <- 46.853; park_lon <- -121.760    # Mount Rainier
radii_km <- c(500, 1000)
min_records <- 30
trunc_q <- c(0.05, 0.95)

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
# 1. Occurrence records with climate
# -----------------------------------------------------------------------------
occ <- read_csv(module_dir("data_clean", "gbif_occurrences_with_bioclim.csv"),
                col_select = c(species, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, MAT, MAP),
                show_col_types = FALSE) %>%
  rename(lat = decimalLatitude, lon = decimalLongitude) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(MAT), !is.na(MAP),
         is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 1e5) %>%
  mutate(species = str_replace_all(species, "_", " "))

# one record per species per 30-arc-second cell
occ <- occ %>%
  mutate(cell_x = floor(lon * 120), cell_y = floor(lat * 120)) %>%
  distinct(species, cell_x, cell_y, .keep_all = TRUE)

# great-circle distance to the park (km)
haversine <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad; dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}
occ <- occ %>% mutate(dist_km = haversine(park_lat, park_lon, lat, lon))
cat("Thinned records:", nrow(occ), "; species:", n_distinct(occ$species), "\n")

# -----------------------------------------------------------------------------
# 2. Affinities: truncated means, range-wide and regional
# -----------------------------------------------------------------------------
trunc_mean <- function(x) { q <- quantile(x, trunc_q, na.rm = TRUE); mean(x[x >= q[1] & x <= q[2]], na.rm = TRUE) }
affinity <- function(d) d %>% group_by(species) %>%
  summarise(n_records = n(), MAT_mean = trunc_mean(MAT), MAP_mean = trunc_mean(MAP), .groups = "drop")

range_wide <- affinity(occ)

# Check against the main affinity table (same species names expected)
main <- read_csv(module_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  mutate(species = str_replace_all(species, "_", " ")) %>% distinct(species, .keep_all = TRUE)
chk <- inner_join(select(main, species, MAT_main = MAT_mean, MAP_main = MAP_mean), range_wide, by = "species")
cat("\nRecomputed vs main range-wide affinities (", nrow(chk), "of", nrow(main), "species matched):\n",
    "  MAT r =", round(cor(chk$MAT_main, chk$MAT_mean), 4), "; max |diff| =", round(max(abs(chk$MAT_main - chk$MAT_mean)), 3), "°C\n",
    "  MAP r =", round(cor(chk$MAP_main, chk$MAP_mean), 4), "; max |diff| =", round(max(abs(chk$MAP_main - chk$MAP_mean)), 1), "mm\n")
unmatched <- setdiff(main$species, range_wide$species)
if (length(unmatched)) cat("Species of the main table without records here (names may differ):\n ", paste(unmatched, collapse = ", "), "\n")
if (cor(chk$MAT_main, chk$MAT_mean) < 0.99) warning("Recomputed range-wide affinities differ from the main table; check the filters before interpreting")

# Species of the shared pools (species recorded at both times in each stratum)
shared <- map_dfr(c("overstory", "understory"), function(st) {
  d <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", paste0(st, "_long.csv")), show_col_types = FALSE)
  tcol <- intersect(c("period", "time", "survey", "census"), names(d))[1]
  d %>% filter(occurrence > 0) %>% group_by(species) %>% filter(n_distinct(.data[[tcol]]) == 2) %>%
    ungroup() %>% distinct(species) %>% mutate(species = str_replace_all(species, "_", " "), stratum = st)
})

# Regional affinities, with range-wide values for species with too few regional records
regional <- map_dfr(radii_km, function(r) {
  reg <- affinity(filter(occ, dist_km <= r)) %>% rename(n_regional = n_records, MAT_reg = MAT_mean, MAP_reg = MAP_mean)
  out <- main %>% select(species, MAT_rw = MAT_mean, MAP_rw = MAP_mean) %>%
    left_join(reg, by = "species") %>%
    mutate(n_regional = replace_na(n_regional, 0L),
           regional_used = n_regional >= min_records,
           MAT_mean = if_else(regional_used, MAT_reg, MAT_rw),
           MAP_mean = if_else(regional_used, MAP_reg, MAP_rw),
           radius_km = r)
  write_csv(select(out, species, MAT_mean, MAP_mean, n_regional, regional_used),
            module_dir("outputs", paste0("species_climatic_affinity_regional_", r, "km.csv")))
  out
})

# -----------------------------------------------------------------------------
# 3. Comparison with range-wide affinities (shared-pool species)
# -----------------------------------------------------------------------------
cmp <- regional %>% inner_join(shared, by = "species", relationship = "many-to-many")
summ <- cmp %>% group_by(radius_km, stratum) %>%
  summarise(n_species = n_distinct(species),
            n_regional = sum(regional_used),
            MAT_r = cor(MAT_rw[regional_used], MAT_reg[regional_used]),
            MAT_rho = cor(MAT_rw[regional_used], MAT_reg[regional_used], method = "spearman"),
            MAT_mean_diff = mean(MAT_reg[regional_used] - MAT_rw[regional_used]),
            MAP_r = cor(MAP_rw[regional_used], MAP_reg[regional_used]),
            MAP_rho = cor(MAP_rw[regional_used], MAP_reg[regional_used], method = "spearman"),
            MAP_mean_diff = mean(MAP_reg[regional_used] - MAP_rw[regional_used]),
            # thermal-moisture affinity correlation with regional affinities
            MAT_MAP_r_regional = cor(MAT_mean, MAP_mean),
            MAT_MAP_r_rangewide = cor(MAT_rw, MAP_rw),
            .groups = "drop") %>%
  mutate(across(where(is.double), ~ signif(.x, 3)))
write_csv(summ, module_dir("outputs", "regional_affinity_comparison.csv"))
cat("\nRegional vs range-wide affinities (shared pools; correlations over species with regional affinities):\n")
print(summ, width = Inf)

cat("\nLargest changes in moisture affinity (", max(radii_km), "km radius, shared pools):\n")
cmp %>% filter(radius_km == max(radii_km), regional_used) %>% distinct(species, .keep_all = TRUE) %>%
  mutate(dMAP = MAP_reg - MAP_rw, dMAT = MAT_reg - MAT_rw) %>% arrange(desc(abs(dMAP))) %>%
  select(species, n_regional, MAT_rw, MAT_reg, dMAT, MAP_rw, MAP_reg, dMAP) %>%
  mutate(across(where(is.double), ~ round(.x, 1))) %>% print(n = 15)

# -----------------------------------------------------------------------------
# 4. Figure: regional vs range-wide affinities
# -----------------------------------------------------------------------------
fig_d <- cmp %>% filter(regional_used) %>% distinct(radius_km, species, .keep_all = TRUE) %>%
  select(radius_km, species, MAT_rw, MAT_reg, MAP_rw, MAP_reg) %>%
  pivot_longer(-c(radius_km, species), names_to = c("var", ".value"), names_pattern = "(MAT|MAP)_(rw|reg)") %>%
  mutate(var = recode(var, MAT = "Thermal affinity (MAT, °C)", MAP = "Moisture affinity (MAP, mm)"),
         radius = paste0(radius_km, " km"))
lab_d <- fig_d %>% group_by(radius, var) %>%
  summarise(r = cor(rw, reg), n = n(), .groups = "drop") %>%
  mutate(lab = sprintf("r = %.2f (n = %d)", r, n))
p <- ggplot(fig_d, aes(rw, reg)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey50") +
  geom_point(alpha = 0.6, size = 1.6, colour = "grey25") +
  geom_text(data = lab_d, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.3, size = 3.6, inherit.aes = FALSE) +
  facet_wrap(radius ~ var, scales = "free", ncol = 2) +
  labs(x = "Range-wide affinity (North America)", y = "Regional affinity") +
  theme_mora()
fig_file <- module_dir("figures", "regional_vs_rangewide_affinity.png")
ggsave(fig_file, p, width = 9, height = 8.5, dpi = 300, bg = "white")
cat("\nFigure written to", fig_file, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
