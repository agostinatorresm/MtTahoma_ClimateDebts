# =============================================================================
# 07_ecoregion_affinity_sensitivity.R
#
# Pipeline: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/
# Author: Agostina Torres
#           (after the script that writes gbif_occurrences_with_bioclim.csv)
#
# Sensitivity of the community climatic-affinity shifts to the geographic extent
# of the occurrence records used to estimate species' affinities. Regional
# extents are defined with the CEC North American terrestrial ecoregions
# (Level III), so that their limits follow climate and vegetation rather than a
# fixed distance from the study area:
#   - Pacific Northwest (PNW): Cascade ranges and Pacific slope forests
#     (6.2.5, 6.2.7, 6.2.8, 6.2.11, 7.1.5-7.1.9); Rocky Mountain and interior
#     ecoregions are excluded.
#   - Western forests: all Level III ecoregions of Level II 6.2 (Western
#     Cordillera) and 7.1 (Marine West Coast Forest).
# Continent-wide affinities use all records (North America).
#
# Records are filtered and thinned as for the main affinities (coordinate
# uncertainty < 100 km; one record per species per 30-arc-second cell). Regional
# affinities require at least min_records records within the region; species
# below the threshold keep their continent-wide values.
#
# Steps
#   1. Regions and map (extents and records of a focal species).
#   2. Point affinities (truncated 5th-95th percentile means) per extent.
#   3. Community rates: for each species, 100 values of MAT and MAP are sampled
#      from its truncated distribution; plot CTI and CMI (abundance-weighted
#      means; overstory relative basal area, understory cover; shared pool) are
#      the median over resamples; rates = change / survey interval (per decade).
#   4. Attribution of the difference in mean rate between extents to species:
#      mean plot rate = sum_i mean(dp_i / gap) * a_i, so replacing species i's
#      continent-wide affinity by its regional affinity changes the mean rate by
#      mean(dp_i / gap) * (a_i,regional - a_i,continent-wide) (point affinities;
#      exact and additive).
#
# Inputs : 02_integrate_gbif_chelsa/data_clean/gbif_occurrences_with_bioclim.csv
#          02_integrate_gbif_chelsa/data_raw/ecoregions/**/NA_CEC_Eco_Level3.shp
#          01_data_prep/01_cleaning/outputs/<stratum>_long.csv
# Outputs: 02_integrate_gbif_chelsa/outputs/species_climatic_affinity_<extent>.csv
#          02_integrate_gbif_chelsa/outputs/ecoregion_affinity_index_rates.csv
#          02_integrate_gbif_chelsa/outputs/ecoregion_affinity_index_rates_plot_cor.csv
#          02_integrate_gbif_chelsa/outputs/ecoregion_affinity_rate_attribution.csv
#          02_integrate_gbif_chelsa/figures/ecoregion_extents_map.png
#          02_integrate_gbif_chelsa/figures/ecoregion_affinity_index_rates.png
# =============================================================================

pkgs <- c("tidyverse", "here", "sf", "rnaturalearth", "rnaturalearthdata", "patchwork", "hrbrthemes", "ggtext")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install first: install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))")

library(tidyverse)
library(here)
library(sf)
library(rnaturalearth)
library(patchwork)
library(hrbrthemes)
sf_use_s2(FALSE)

module_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
dir.create(module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

park_lat <- 46.853; park_lon <- -121.760
crs_aeqd <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m", park_lat, park_lon)
crs_map  <- "+proj=laea +lat_0=50 +lon_0=-112 +datum=WGS84 +units=m"
min_records <- 30
n_res <- 100
trunc_q <- c(0.05, 0.95)
focal_sp <- "Pseudotsuga menziesii"
STRATA <- tibble(stratum = c("overstory", "understory"), weight = c("basal_area_spp_perc", "cover"))
set.seed(2026)

pnw_codes <- c("6.2.5", "6.2.7", "6.2.8", "6.2.11", "7.1.5", "7.1.6", "7.1.7", "7.1.8", "7.1.9")
set_levels <- c(rangewide = "Continent-wide", western = "Western forests", pnw = "Pacific Northwest")
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
theme_map <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title = element_blank(), axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(size = 13, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9))
}
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))

# -----------------------------------------------------------------------------
# 1. Regions, records and map
# -----------------------------------------------------------------------------
eco_shp <- list.files(module_dir("data_raw", "ecoregions"), pattern = "NA_CEC_Eco_Level3\\.shp$", recursive = TRUE, full.names = TRUE)[1]
if (is.na(eco_shp)) stop("NA_CEC_Eco_Level3.shp not found in data_raw/ecoregions")
eco <- st_read(eco_shp, quiet = TRUE) %>%
  filter(NA_L2CODE %in% c("6.2", "7.1")) %>%
  st_make_valid() %>% st_transform(crs_aeqd)

regions <- list(western = unique(eco$NA_L3CODE), pnw = pnw_codes) %>%
  imap(~ eco %>% filter(NA_L3CODE %in% .x) %>% st_union() %>% st_make_valid() %>% st_simplify(dTolerance = 1000))
cat("Level III ecoregions in the PNW extent:\n")
eco %>% st_drop_geometry() %>% filter(NA_L3CODE %in% pnw_codes) %>% distinct(NA_L3CODE, NA_L3NAME) %>% arrange(NA_L3CODE) %>% print()

occ <- read_csv(module_dir("data_clean", "gbif_occurrences_with_bioclim.csv"),
                col_select = c(species, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, MAT, MAP),
                show_col_types = FALSE) %>%
  rename(lat = decimalLatitude, lon = decimalLongitude) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(MAT), !is.na(MAP),
         is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 1e5) %>%
  mutate(species = str_replace_all(species, "_", " "),
         cell_x = floor(lon * 120), cell_y = floor(lat * 120)) %>%
  distinct(species, cell_x, cell_y, .keep_all = TRUE)
occ_pts <- st_as_sf(occ, coords = c("lon", "lat"), crs = 4326) %>% st_transform(crs_aeqd)
occ <- occ %>% mutate(in_western = lengths(st_intersects(occ_pts, regions$western)) > 0,
                      in_pnw = lengths(st_intersects(occ_pts, regions$pnw)) > 0)
cat("\nThinned records: all", nrow(occ), "; Western forests", sum(occ$in_western), "; PNW", sum(occ$in_pnw), "\n")

# Map: extents (A) and records of the focal species (B)
land <- ne_countries(scale = 50, continent = "North America", returnclass = "sf") %>% st_transform(crs_map)
reg_map <- imap_dfr(regions, ~ st_sf(extent = set_levels[.y], geometry = st_transform(.x, crs_map))) %>%
  mutate(extent = factor(extent, levels = set_levels[c("western", "pnw")]))
park_m <- st_sfc(st_point(c(park_lon, park_lat)), crs = 4326) %>% st_transform(crs_map)
foc <- occ %>% filter(species == focal_sp) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% st_transform(crs_map)
bb <- st_bbox(st_union(st_geometry(foc), st_geometry(filter(reg_map, extent == set_levels["pnw"]))))
xl <- c(bb["xmin"], bb["xmax"]) + c(-2e5, 2e5); yl <- c(bb["ymin"], bb["ymax"]) + c(-2e5, 2e5)

trunc_mean <- function(x) { q <- quantile(x, trunc_q); mean(x[x >= q[1] & x <= q[2]]) }
foc_tab <- occ %>% filter(species == focal_sp) %>%
  summarise(across(c(MAT, MAP), list(rangewide = ~ trunc_mean(.x),
                                     western = ~ trunc_mean(.x[in_western]),
                                     pnw = ~ trunc_mean(.x[in_pnw]))),
            n_rangewide = n(), n_western = sum(in_western), n_pnw = sum(in_pnw)) %>%
  pivot_longer(everything(), names_to = c(".value", "extent"), names_sep = "_")
cat("\n", focal_sp, ": records and affinities per extent\n", sep = ""); print(mutate(foc_tab, MAT = round(MAT, 1), MAP = round(MAP)))

pA <- ggplot() +
  geom_sf(data = land, fill = "grey93", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = reg_map, aes(fill = extent), colour = NA, alpha = 0.6) +
  geom_sf(data = park_m, shape = 24, size = 3, fill = "black") +
  coord_sf(xlim = xl, ylim = yl, crs = crs_map, expand = FALSE) +
  scale_fill_manual(values = setNames(c("#A9CCE3", "#1A759F"), set_levels[c("western", "pnw")]), name = "Regional extent") +
  labs(title = "A. Regional extents (CEC ecoregions)") +
  theme_map() + theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
                      legend.background = element_rect(fill = "white", colour = NA))
lab_B <- with(foc_tab, paste0(set_levels[extent], ": ", round(MAP), " mm (n = ", format(n, big.mark = ","), ")", collapse = "\n"))
pB <- ggplot() +
  geom_sf(data = land, fill = "grey93", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = reg_map, aes(linetype = extent), fill = NA, colour = "black", linewidth = 0.4, show.legend = FALSE) +
  geom_sf(data = arrange(foc, MAP), aes(colour = MAP), size = 0.5, alpha = 0.7) +
  geom_sf(data = park_m, shape = 24, size = 3, fill = "white") +
  annotate("label", x = xl[2] - 5e4, y = yl[2] - 5e4, label = lab_B, hjust = 1, vjust = 1, size = 3.2, fill = alpha("white", 0.85)) +
  coord_sf(xlim = xl, ylim = yl, crs = crs_map, expand = FALSE) +
  scale_colour_gradientn(colours = c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D"), trans = "log10", name = "MAP (mm)") +
  scale_linetype_manual(values = setNames(c("22", "solid"), set_levels[c("western", "pnw")])) +
  labs(title = paste0("B. Records of *", focal_sp, "*")) +
  theme_map() + theme(plot.title = ggtext::element_markdown(size = 13, face = "bold"))
map_file <- module_dir("figures", "ecoregion_extents_map.png")
ggsave(map_file, pA + pB, width = 13, height = 7, dpi = 300, bg = "white")
open_fig(map_file)

# -----------------------------------------------------------------------------
# 2. Point affinities per extent
# -----------------------------------------------------------------------------
point_aff <- function(d) d %>% group_by(species) %>%
  summarise(n = n(), MAT = trunc_mean(MAT), MAP = trunc_mean(MAP), .groups = "drop")
rw_pt <- point_aff(occ)
aff_pt <- bind_rows(
  mutate(rw_pt, affinity_set = "rangewide", regional_used = FALSE),
  map_dfr(c("western", "pnw"), function(s) {
    reg <- point_aff(filter(occ, .data[[paste0("in_", s)]])) %>% filter(n >= min_records)
    out <- rw_pt %>% select(species, MAT_rw = MAT, MAP_rw = MAP) %>%
      left_join(reg, by = "species") %>%
      mutate(regional_used = !is.na(n), n = replace_na(n, 0L),
             MAT = coalesce(MAT, MAT_rw), MAP = coalesce(MAP, MAP_rw), affinity_set = s) %>%
      select(species, n, MAT, MAP, regional_used, affinity_set)
    write_csv(rename(out, n_regional = n, MAT_mean = MAT, MAP_mean = MAP) %>% select(-affinity_set),
              module_dir("outputs", paste0("species_climatic_affinity_ecoregion_", s, ".csv")))
    out
  }))

# -----------------------------------------------------------------------------
# 3. Community rates with resampled affinities
# -----------------------------------------------------------------------------
draw <- function(x) { q <- quantile(x, trunc_q); xt <- x[x >= q[1] & x <= q[2]]; xt[sample.int(length(xt), n_res, replace = TRUE)] }
resample_set <- function(d) d %>% group_by(species) %>%
  summarise(rep = list(seq_len(n_res)), MAT = list(draw(MAT)), MAP = list(draw(MAP)), .groups = "drop") %>%
  unnest(c(rep, MAT, MAP))
rw_res <- resample_set(occ)
aff_res <- bind_rows(
  mutate(rw_res, affinity_set = "rangewide"),
  map_dfr(c("western", "pnw"), function(s) {
    reg <- occ %>% filter(.data[[paste0("in_", s)]]) %>% group_by(species) %>% filter(n() >= min_records) %>% ungroup() %>% resample_set()
    bind_rows(reg, filter(rw_res, !species %in% reg$species)) %>% mutate(affinity_set = s)
  }))

read_stratum <- function(stratum, weight) {
  d <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", paste0(stratum, "_long.csv")), show_col_types = FALSE) %>%
    mutate(species = str_replace_all(species, "_", " "))
  tcol <- intersect(c("period", "time", "survey", "census"), names(d))[1]
  shared <- d %>% filter(occurrence > 0) %>% group_by(species) %>% filter(n_distinct(.data[[tcol]]) == 2) %>% pull(species) %>% unique()
  d <- d %>% filter(species %in% shared, occurrence > 0, !is.na(.data[[weight]]), .data[[weight]] > 0) %>%
    transmute(plot_id, period = .data[[tcol]], year, species, w = .data[[weight]])
  gap <- d %>% group_by(plot_id, period) %>% summarise(year = mean(year), .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = year) %>%
    transmute(plot_id, gap = (present - past) / 10) %>% filter(is.finite(gap), gap > 0)
  list(d = semi_join(d, gap, by = "plot_id"), gap = gap, shared = shared)
}
strata <- pmap(STRATA, read_stratum) %>% set_names(STRATA$stratum)

# species of the shared pools with regional affinities
cat("\nShared-pool species with regional affinities (>=", min_records, "records in the region):\n")
imap_dfr(strata, function(s, st) aff_pt %>% filter(species %in% s$shared, affinity_set != "rangewide") %>%
           group_by(affinity_set) %>% summarise(stratum = st, n_species = n(), n_regional = sum(regional_used), .groups = "drop")) %>% print()

plot_rates <- imap_dfr(strata, function(s, st) {
  s$d %>% select(-year) %>%
    inner_join(aff_res, by = "species", relationship = "many-to-many") %>%
    group_by(affinity_set, plot_id, period, rep) %>%
    summarise(CTI = weighted.mean(MAT, w), CMI = weighted.mean(MAP, w), .groups = "drop") %>%
    group_by(affinity_set, plot_id, period) %>%
    summarise(CTI = median(CTI), CMI = median(CMI), .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = c(CTI, CMI)) %>%
    inner_join(s$gap, by = "plot_id") %>%
    transmute(stratum = st, affinity_set, plot_id,
              CTI = (CTI_present - CTI_past) / gap, CMI = (CMI_present - CMI_past) / gap) %>%
    filter(is.finite(CTI), is.finite(CMI))
}) %>% pivot_longer(c(CTI, CMI), names_to = "index", values_to = "rate")

mean_ci <- function(x) { t <- t.test(x); tibble(n = length(x), mean = mean(x), lower_95 = t$conf.int[1], upper_95 = t$conf.int[2]) }
summ <- plot_rates %>% group_by(stratum, index, affinity_set) %>%
  group_modify(~ mean_ci(.x$rate)) %>% ungroup() %>%
  mutate(across(c(mean, lower_95, upper_95), ~ signif(.x, 3)))
write_csv(summ, module_dir("outputs", "ecoregion_affinity_index_rates.csv"))
cat("\nMean rates (CTI, °C/decade; CMI, mm/decade), median over", n_res, "affinity resamples:\n")
print(arrange(summ, stratum, desc(index), match(affinity_set, names(set_levels))), n = Inf)

pc <- plot_rates %>% pivot_wider(names_from = affinity_set, values_from = rate) %>%
  group_by(stratum, index) %>%
  summarise(r_western = cor(rangewide, western, use = "complete.obs"),
            r_pnw = cor(rangewide, pnw, use = "complete.obs"), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
write_csv(pc, module_dir("outputs", "ecoregion_affinity_index_rates_plot_cor.csv"))
cat("\nCorrelation of plot-level rates with those based on continent-wide affinities:\n"); print(pc)

fig_d <- summ %>%
  mutate(set = factor(set_levels[affinity_set], levels = set_levels),
         index = factor(index, levels = c("CTI", "CMI")),
         panel = paste(str_to_title(stratum), if_else(index == "CTI", "ΔCTI (°C decade⁻¹)", "ΔCMI (mm decade⁻¹)"), sep = " – "),
         panel = factor(panel, levels = c(paste("Overstory", c("ΔCTI (°C decade⁻¹)", "ΔCMI (mm decade⁻¹)"), sep = " – "),
                                          paste("Understory", c("ΔCTI (°C decade⁻¹)", "ΔCMI (mm decade⁻¹)"), sep = " – "))))
p_rates <- ggplot(fig_d, aes(set, mean, colour = index)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), linewidth = 1.6, lineend = "round") +
  geom_point(size = 3) +
  scale_colour_manual(values = col_index, guide = "none") +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  labs(x = "Extent of occurrence records", y = "Mean rate across plots (95% CI)") +
  theme_mora()
rates_file <- module_dir("figures", "ecoregion_affinity_index_rates.png")
ggsave(rates_file, p_rates, width = 9, height = 7, dpi = 300, bg = "white")
open_fig(rates_file)

# -----------------------------------------------------------------------------
# 4. Attribution of the difference in mean rate to species (point affinities)
# -----------------------------------------------------------------------------
dp <- imap_dfr(strata, function(s, st) {
  n_plots <- n_distinct(s$gap$plot_id)
  s$d %>% select(-year) %>%
    group_by(plot_id, period) %>% mutate(p = w / sum(w)) %>% ungroup() %>%
    select(plot_id, species, period, p) %>%
    pivot_wider(names_from = period, values_from = p, values_fill = 0) %>%
    inner_join(s$gap, by = "plot_id") %>%
    group_by(species) %>%
    summarise(mean_dp_rate = sum((present - past) / gap) / n_plots, .groups = "drop") %>%   # absent plots contribute 0
    mutate(stratum = st)
})
attr_tab <- aff_pt %>% select(species, affinity_set, MAT, MAP) %>%
  pivot_longer(c(MAT, MAP), names_to = "var", values_to = "a") %>%
  pivot_wider(names_from = affinity_set, values_from = a) %>%
  inner_join(dp, by = "species", relationship = "many-to-many") %>%
  mutate(index = if_else(var == "MAT", "CTI", "CMI"),
         attr_western = mean_dp_rate * (western - rangewide),
         attr_pnw = mean_dp_rate * (pnw - rangewide))
write_csv(attr_tab, module_dir("outputs", "ecoregion_affinity_rate_attribution.csv"))

cat("\nMean rates with point affinities (sum over species of mean dp rate x affinity):\n")
attr_tab %>% group_by(stratum, index) %>%
  summarise(rangewide = sum(mean_dp_rate * rangewide), western = sum(mean_dp_rate * western),
            pnw = sum(mean_dp_rate * pnw), .groups = "drop") %>%
  mutate(across(where(is.double), ~ signif(.x, 3))) %>% print()

cat("\nSpecies contributing most to the difference PNW - continent-wide (mean_dp_rate in % relative abundance per decade):\n")
attr_tab %>% group_by(stratum, index) %>% slice_max(abs(attr_pnw), n = 6) %>% ungroup() %>%
  transmute(stratum, index, species, mean_dp_rate = round(mean_dp_rate * 100, 2),
            a_rangewide = round(rangewide, 1), a_pnw = round(pnw, 1),
            attr_pnw = signif(attr_pnw, 2), attr_western = signif(attr_western, 2)) %>%
  print(n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 5. Regional vs continent-wide affinities of shared-pool species (S2)
# -----------------------------------------------------------------------------
shared_sp <- imap_dfr(strata, ~ tibble(species = .x$shared, stratum = .y))
aff_cmp <- aff_pt %>% filter(affinity_set != "rangewide", regional_used) %>%
  select(species, affinity_set, MAT, MAP) %>%
  inner_join(select(rw_pt, species, MAT_rw = MAT, MAP_rw = MAP), by = "species") %>%
  inner_join(shared_sp, by = "species", relationship = "many-to-many") %>%
  group_by(stratum, affinity_set) %>%
  summarise(n_species = n(),
            MAT_r = cor(MAT, MAT_rw), MAT_mean_diff = mean(MAT - MAT_rw),
            MAP_r = cor(MAP, MAP_rw), MAP_mean_diff = mean(MAP - MAP_rw), .groups = "drop") %>%
  mutate(across(where(is.double), ~ signif(.x, 3)))
write_csv(aff_cmp, module_dir("outputs", "ecoregion_vs_rangewide_affinity_comparison.csv"))
cat("\nRegional vs continent-wide affinities (shared-pool species with regional affinities):\n"); print(aff_cmp)
