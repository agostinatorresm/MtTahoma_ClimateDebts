# =============================================================================
# 02_climate_change_spatial_structure.R
#
# Pipeline: pipeline/01_data_prep/03_climate/02_explore_climate/code/
# Author: Agostina Torres
#
# Describes how plot-level climate change (present - past) varies across the
# landscape: along elevation, from west to east and from south to north, and
# between the west and east sides of the mountain (split at the longitude of
# the summit). Used to describe the spatial structure of climate change in the
# Results.
#
# Analyses, for each change variable (Delta = present - past):
#   - Spearman correlation with elevation, easting and northing;
#   - linear model Delta ~ elevation + easting + northing with predictors per
#     standard deviation (standardised effects, comparable among predictors)
#     and per natural unit (per 100 m of elevation, per 10 km);
#   - mean change on the west and east sides, and the east - west difference
#     at equal elevation (Delta ~ side + elevation).
#   - linear trend of each change along elevation, across all plots and on
#     each side (slope per 100 m and R^2; reported in the figure caption).
# Only the plots used in the analyses are included (the 324 plots of the
# overstory climatic-affinity analyses, which contain the 316 understory plots).
# Plots are spatially autocorrelated, so p-values are optimistic; effects and
# R^2 are used descriptively.
#
# Inputs : 02_explore_climate/data_clean/climate_past_present_wide.csv
#          01_data_prep/01_cleaning/outputs/plot_info.csv (coordinates, elevation)
#          02_analysis/01_thermophilization/01_overstory/outputs/overstory_thermo_xero_rates.csv
#            (plots used in the analyses)
# Outputs: 02_explore_climate/outputs/climate_change_spatial_correlations.csv
#          02_explore_climate/outputs/climate_change_spatial_models.csv
#          02_explore_climate/outputs/climate_change_by_side.csv
#          02_explore_climate/outputs/climate_change_elevation_trend.csv
#          02_explore_climate/figures/climate_change_spatial_structure.png
# =============================================================================

library(tidyverse)
library(here)
library(sf)
library(hrbrthemes)
library(patchwork)

clim_dir <- function(...) here("pipeline", "01_data_prep", "03_climate", "02_explore_climate", ...)
out_dir <- clim_dir("outputs"); fig_dir <- clim_dir("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

summit <- c(lon = -121.7603, lat = 46.8529)          # Mt. Tahoma summit
crs_utm <- "EPSG:32610"                               # UTM zone 10N

# Change variables: name in climate_past_present_wide -> label (those present are used)
vars <- c(delta_MAT = "ΔMAT (°C)", delta_MAP = "ΔMAP (mm)", delta_MAP_pct = "ΔMAP (%)",
          delta_temp_JJA = "ΔT summer (°C)", delta_temp_DJF = "ΔT winter (°C)",
          delta_PAS = "ΔPAS (mm)")

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.text = element_text(colour = "black", size = 11),
          axis.title.x = element_text(colour = "black", size = 12, hjust = 0.5),
          axis.title.y = element_text(colour = "black", size = 12, hjust = 0.5),
          plot.title = element_text(colour = "black", size = 13, face = "bold"))
}
pal_moist <- c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D")
col_side  <- c(West = "#1A759F", East = "#C98B2B")

# -----------------------------------------------------------------------------
# 1. Data: climate change, elevation and coordinates per plot
# -----------------------------------------------------------------------------
clim <- read_csv(clim_dir("data_clean", "climate_past_present_wide.csv"), show_col_types = FALSE) %>%
  mutate(plot_id = as.character(plot_id))
if (all(c("delta_MAP", "MAP_past") %in% names(clim))) clim <- mutate(clim, delta_MAP_pct = 100 * delta_MAP / MAP_past)
vars <- vars[names(vars) %in% names(clim)]

# Plots used in the analyses
analysis_plots <- read_csv(here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", "outputs",
                                "overstory_thermo_xero_rates.csv"), show_col_types = FALSE) %>%
  distinct(plot_id = as.character(plot_id)) %>% pull(plot_id)

info <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", "plot_info.csv"), show_col_types = FALSE)
pick <- function(cand) { x <- intersect(cand, names(info)); if (!length(x)) stop("None of ", paste(cand, collapse = "/"), " in plot_info"); x[1] }
geo <- info %>%
  transmute(plot_id = as.character(.data[[pick(c("plot_id", "plot", "Plot", "PlotID"))]]),
            lat = .data[[pick(c("lat", "latitude", "Latitude", "LAT", "lat_dd"))]],
            lon = .data[[pick(c("lon", "long", "longitude", "Longitude", "LON", "lon_dd"))]],
            elev = .data[[pick(c("elevation", "elev", "Elevation", "elev_m", "ELEV"))]]) %>%
  distinct(plot_id, .keep_all = TRUE) %>% filter(!is.na(lat), !is.na(lon)) %>%
  filter(plot_id %in% analysis_plots)   # before projecting, so coordinates and plots stay aligned
xy <- st_coordinates(st_transform(st_as_sf(geo, coords = c("lon", "lat"), crs = 4326), crs_utm))
summit_xy <- st_coordinates(st_transform(st_sfc(st_point(summit), crs = 4326), crs_utm))

plots <- geo %>%
  mutate(easting_km = (xy[, 1] - summit_xy[1]) / 1000,               # km east of the summit
         northing_km = (xy[, 2] - summit_xy[2]) / 1000,              # km north of the summit
         side = factor(if_else(lon < summit["lon"], "West", "East"), levels = c("West", "East"))) %>%
  inner_join(select(clim, plot_id, all_of(names(vars))), by = "plot_id")
missing_clim <- setdiff(analysis_plots, plots$plot_id)
if (length(missing_clim)) cat("Analysis plots without climate or coordinates:", paste(missing_clim, collapse = ", "), "\n")
cat("Plots:", nrow(plots), "(west", sum(plots$side == "West"), ", east", sum(plots$side == "East"), ")\n")

long <- plots %>% pivot_longer(all_of(names(vars)), names_to = "variable", values_to = "delta") %>%
  mutate(label = vars[variable])

# -----------------------------------------------------------------------------
# 2. Correlations and models with elevation, easting and northing
# -----------------------------------------------------------------------------
cors <- long %>% group_by(variable, label) %>%
  summarise(mean = mean(delta), sd = sd(delta), min = min(delta), max = max(delta),
            rho_elevation = cor(delta, elev, method = "spearman"),
            rho_easting = cor(delta, easting_km, method = "spearman"),
            rho_northing = cor(delta, northing_km, method = "spearman"), .groups = "drop")

models <- long %>% group_by(variable, label) %>%
  group_modify(function(d, k) {
    f_sd <- lm(scale(delta) ~ scale(elev) + scale(easting_km) + scale(northing_km), d)       # standardised
    f_un <- lm(delta ~ I(elev / 100) + I(easting_km / 10) + I(northing_km / 10), d)          # natural units
    tibble(predictor = c("elevation", "easting", "northing"),
           unit = c("per 100 m", "per 10 km east", "per 10 km north"),
           effect_per_unit = coef(f_un)[-1], effect_std = coef(f_sd)[-1],
           p = summary(f_un)$coefficients[-1, 4], r2 = summary(f_un)$r.squared)
  }) %>% ungroup()

side_tab <- long %>% group_by(variable, label) %>%
  group_modify(function(d, k) {
    f <- lm(delta ~ side + elev, d)
    tibble(mean_west = mean(d$delta[d$side == "West"]), mean_east = mean(d$delta[d$side == "East"]),
           diff_east_minus_west_same_elev = coef(f)["sideEast"],
           p = summary(f)$coefficients["sideEast", 4])
  }) %>% ungroup()

# Linear trend along elevation, across all plots and on each side
elev_trend <- long %>%
  filter(variable %in% c("delta_MAT", "delta_MAP")) %>%
  mutate(group = as.character(side)) %>%
  bind_rows(mutate(filter(long, variable %in% c("delta_MAT", "delta_MAP")), group = "All")) %>%
  group_by(variable, label, group) %>%
  group_modify(function(d, k) {
    f <- lm(delta ~ I(elev / 100), d)
    tibble(n = nrow(d), slope_per_100m = coef(f)[2],
           ci_lower = confint(f)[2, 1], ci_upper = confint(f)[2, 2],
           r2 = summary(f)$r.squared, p = summary(f)$coefficients[2, 4],
           elev_range_m = diff(range(d$elev)),
           change_over_range = coef(f)[2] * diff(range(d$elev)) / 100)
  }) %>% ungroup() %>%
  mutate(group = factor(group, levels = c("All", "West", "East"))) %>% arrange(variable, group)

write_csv(mutate(cors, across(where(is.double), ~ signif(.x, 3))), file.path(out_dir, "climate_change_spatial_correlations.csv"))
write_csv(mutate(models, across(where(is.double), ~ signif(.x, 3))), file.path(out_dir, "climate_change_spatial_models.csv"))
write_csv(mutate(side_tab, across(where(is.double), ~ signif(.x, 3))), file.path(out_dir, "climate_change_by_side.csv"))
write_csv(mutate(elev_trend, across(where(is.double), ~ signif(.x, 3))), file.path(out_dir, "climate_change_elevation_trend.csv"))

cat("\nClimate change (present - past): summary and Spearman correlations\n")
print(mutate(cors, across(where(is.double), ~ signif(.x, 3))), width = Inf)
cat("\nDelta ~ elevation + easting + northing (effects per unit and per SD)\n")
print(mutate(models, across(where(is.double), ~ signif(.x, 3))), n = Inf, width = Inf)
cat("\nWest vs east side of the summit (difference at equal elevation)\n")
print(mutate(side_tab, across(where(is.double), ~ signif(.x, 3))), width = Inf)
cat("\nChange along elevation (Delta ~ elevation), all plots and by side\n")
print(mutate(elev_trend, across(where(is.double), ~ signif(.x, 3))), width = Inf)

# -----------------------------------------------------------------------------
# 3. Figure: maps of ΔMAT and ΔMAP, and change along elevation by side
# -----------------------------------------------------------------------------
map_panel <- function(v, title, scale_fill) {
  ggplot(plots, aes(easting_km, northing_km, fill = .data[[v]])) +
    geom_vline(xintercept = 0, linetype = "22", colour = "grey60") +
    geom_point(shape = 21, size = 2.2, colour = "grey30", stroke = 0.2) +
    annotate("point", x = 0, y = 0, shape = 24, size = 3, fill = "black") +
    scale_fill + coord_equal() +
    labs(title = title, x = "Easting from summit (km)", y = "Northing from summit (km)", fill = vars[v]) +
    theme_mora() + theme(legend.position = "right")
}
elev_panel <- function(v, title) {
  ggplot(plots, aes(elev, .data[[v]], colour = side)) +
    geom_point(alpha = 0.6, size = 1.4) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9, alpha = 0.15, aes(fill = side)) +
    scale_colour_manual(values = col_side, name = NULL) + scale_fill_manual(values = col_side, name = NULL) +
    labs(title = title, x = "Elevation (m)", y = vars[v]) +
    theme_mora() + theme(legend.position = "top", legend.justification = "left",          # legend above the panel, clear of the points
                         legend.text = element_text(size = 11), legend.key.width = unit(1.2, "cm"))
}
fig <- (map_panel("delta_MAT", "A. Change in temperature", scale_fill_viridis_c(option = "plasma")) |
          map_panel("delta_MAP", "B. Change in precipitation", scale_fill_gradientn(colours = pal_moist))) /
  (elev_panel("delta_MAT", "C. Temperature change along elevation") |
     elev_panel("delta_MAP", "D. Precipitation change along elevation"))
fig_file <- file.path(fig_dir, "climate_change_spatial_structure.png")
ggsave(fig_file, fig, width = 12, height = 10, dpi = 300, bg = "white")
cat("\nFigure written to", fig_file, "\n")

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)