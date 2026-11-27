# =============================================================================
# 06_map_affinity_extents.R
#
# Pipeline: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/
# Author: Agostina Torres
#           (after 03_regional_affinity_sensitivity.R)
#
# Map of the geographic extents used to estimate species' climatic affinities:
# range-wide (all thinned GBIF records) and regional (records within 500 and
# 1000 km of Mount Rainier). Distances are great-circle distances to the park
# centre; circles are drawn as buffers in an azimuthal equidistant projection
# centred on the park, so they represent true distances.
#   A. Extents on the map of North America.
#   B. Thinned records of a focal species coloured by MAP, showing which part
#      of its climatic distribution enters each affinity set.
#
# Inputs : 02_integrate_gbif_chelsa/data_clean/gbif_occurrences_with_bioclim.csv
# Outputs: 02_integrate_gbif_chelsa/figures/affinity_extents_map.png
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

module_dir <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
dir.create(module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

park_lat <- 46.853; park_lon <- -121.760
radii_km <- c(500, 1000)
focal_sp <- "Pseudotsuga menziesii"
crs_map  <- "+proj=laea +lat_0=50 +lon_0=-112 +datum=WGS84 +units=m"
crs_aeqd <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m", park_lat, park_lon)

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_map <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_blank(), axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(size = 13, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9))
}

# -----------------------------------------------------------------------------
# 1. Base map, park and extents
# -----------------------------------------------------------------------------
land <- ne_countries(scale = 50, continent = "North America", returnclass = "sf") %>% st_transform(crs_map)
park <- st_sfc(st_point(c(park_lon, park_lat)), crs = 4326)
circles <- map_dfr(radii_km, function(r) {
  st_sf(radius = paste(r, "km"), geometry = park %>% st_transform(crs_aeqd) %>% st_buffer(r * 1000, nQuadSegs = 90))
}) %>% st_transform(crs_map) %>% mutate(radius = factor(radius, levels = paste(radii_km, "km")))
park_m <- st_transform(park, crs_map)

# -----------------------------------------------------------------------------
# 2. Records (same filters and thinning as the affinities)
# -----------------------------------------------------------------------------
occ <- read_csv(module_dir("data_clean", "gbif_occurrences_with_bioclim.csv"),
                col_select = c(species, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, MAT, MAP),
                show_col_types = FALSE) %>%
  rename(lat = decimalLatitude, lon = decimalLongitude) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(MAT), !is.na(MAP),
         is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 1e5) %>%
  mutate(species = str_replace_all(species, "_", " "),
         cell_x = floor(lon * 120), cell_y = floor(lat * 120)) %>%
  distinct(species, cell_x, cell_y, .keep_all = TRUE)

foc <- occ %>% filter(species == focal_sp) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  mutate(dist_km = as.numeric(st_distance(geometry, park)) / 1000) %>%
  st_transform(crs_map)

# Affinity (truncated mean MAP) and number of records of the focal species per extent
trunc_mean <- function(x) { q <- quantile(x, c(0.05, 0.95)); mean(x[x >= q[1] & x <= q[2]]) }
foc_summ <- bind_rows(
  tibble(extent = "Range-wide", n = nrow(foc), MAP = trunc_mean(foc$MAP)),
  map_dfr(radii_km, ~ tibble(extent = paste(.x, "km"), n = sum(foc$dist_km <= .x), MAP = trunc_mean(foc$MAP[foc$dist_km <= .x]))))
cat("\n", focal_sp, ": thinned records and moisture affinity (truncated mean MAP, mm) per extent\n", sep = "")
print(mutate(foc_summ, MAP = round(MAP)))

# Map window: all records of the focal species plus the largest circle
bb <- st_bbox(st_union(st_geometry(foc), st_geometry(circles)))
pad <- 2e5
xl <- c(bb["xmin"] - pad, bb["xmax"] + pad); yl <- c(bb["ymin"] - pad, bb["ymax"] + pad)

# -----------------------------------------------------------------------------
# 3. Figure
# -----------------------------------------------------------------------------
col_circ <- c("500 km" = "#1A759F", "1000 km" = "#F26B5B")
base <- function() list(
  geom_sf(data = land, fill = "grey93", colour = "grey60", linewidth = 0.2),
  coord_sf(xlim = xl, ylim = yl, crs = crs_map, expand = FALSE))

pA <- ggplot() + base() +
  geom_sf(data = circles, aes(colour = radius), fill = NA, linewidth = 0.9) +
  geom_sf(data = park_m, shape = 24, size = 3, fill = "black") +
  coord_sf(xlim = xl, ylim = yl, crs = crs_map, expand = FALSE) +
  scale_colour_manual(values = col_circ, name = "Regional extent") +
  labs(title = "A. Extents of the occurrence records") +
  theme_map() + theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
                      legend.background = element_rect(fill = "white", colour = NA))

lab_B <- with(foc_summ, paste0(extent, ": ", round(MAP), " mm (n = ", format(n, big.mark = ","), ")", collapse = "\n"))
pB <- ggplot() + base() +
  geom_sf(data = arrange(foc, MAP), aes(colour = MAP), size = 0.5, alpha = 0.7) +
  geom_sf(data = circles, aes(linetype = radius), colour = "black", fill = NA, linewidth = 0.6, show.legend = FALSE) +
  geom_sf(data = park_m, shape = 24, size = 3, fill = "white") +
  annotate("label", x = xl[2] - 5e4, y = yl[2] - 5e4, label = lab_B, hjust = 1, vjust = 1,
           size = 3.2, label.size = 0, fill = alpha("white", 0.85)) +
  coord_sf(xlim = xl, ylim = yl, crs = crs_map, expand = FALSE) +
  scale_colour_gradientn(colours = c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D"),
                         trans = "log10", name = "MAP (mm)") +
  scale_linetype_manual(values = c("500 km" = "solid", "1000 km" = "22")) +
  labs(title = paste0("B. Records of *", focal_sp, "*")) +
  theme_map() + theme(plot.title = ggtext::element_markdown(size = 13, face = "bold"))

p <- pA + pB + plot_layout(ncol = 2)
fig_file <- module_dir("figures", "affinity_extents_map.png")
ggsave(fig_file, p, width = 13, height = 7, dpi = 300, bg = "white")
cat("\nFigure written to", fig_file, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
