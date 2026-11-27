# =============================================================================
# figureS1_study_area.R
#
# Pipeline: pipeline/05_manuscript_figures/figureS1/
# Author: Agostina Torres
#
# Figure S1.1: study area and survey design.
#   A. Study plots over hillshade and CHELSA V2.1 mean annual temperature
#      (1981-2010, 30 arc-sec; the climatology used for species' affinities),
#      with an inset locating the study area in North America.
#   B. Same extent with CHELSA V2.1 mean annual precipitation (1981-2010).
#   C. Survey design: the sampled units of the past (~1975) and present
#      (2014-2020) surveys. Past: overstory stems tallied in 10-cm DBH classes
#      and understory presence and cover recorded over the whole plot
#      (250-1000 m2). Present: DBH measured per stem in a 500-m2 plot,
#      understory presence in a 50-m2 census and understory cover in four
#      1-m2 quadrats within it. The example plot is marked in A and B
#      and linked to C by zoom lines.
# Plots are those with data in both surveys in at least one layer.
#
# Inputs : 01_data_prep/01_cleaning/outputs/plot_info.csv
#          01_data_prep/01_cleaning/outputs/<stratum>_long.csv
#          figureS1/data_raw/: optional park boundary (*boundary*.shp/.gpkg/.geojson)
# Downloads (first run, cached in figureS1/data_raw/):
#          CHELSA bio1 and bio12 windows around the park (only the window is read)
#          AWS terrain tiles (elevatr) for the hillshade
# Outputs: figureS1/figureS1_study_area.png, figureS1/figureS1_study_area.svg
# =============================================================================

pkgs <- c("tidyverse", "here", "sf", "terra", "tidyterra", "ggnewscale", "ggspatial", "elevatr",
          "rnaturalearth", "rnaturalearthdata", "patchwork", "hrbrthemes", "ragg")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install first: install.packages(c(", paste0('"', missing_pkgs, '"', collapse = ", "), "))")

library(tidyverse)
library(here)
library(sf)
library(terra)
library(tidyterra)
library(ggnewscale)
library(ggspatial)
library(rnaturalearth)
library(patchwork)
library(hrbrthemes)
sf_use_s2(FALSE)

fig_dir <- here("pipeline", "05_manuscript_figures", "figureS1")
raw_dir <- file.path(fig_dir, "data_raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

crs_utm <- "EPSG:32610"                                    # UTM zone 10N
crs_map <- "+proj=laea +lat_0=50 +lon_0=-105 +datum=WGS84 +units=m"
chelsa_url <- function(v) sprintf("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_%s_1981-2010_V.2.1.tif", v)
buffer_m <- 3000                                           # margin around the plots
pal_moist <- c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D")
col_tree <- "#2D5A27"; col_under <- "#8DBF5A"; col_plot <- "#F3EEE2"
emoji_tree <- "\U0001F332"; emoji_leaf <- "\U0001F33F"    # as in Figure 1 (colour emoji rendered by ragg)

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_map <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title = element_blank(), axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          plot.title = element_text(size = 13, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          plot.margin = margin(4, 4, 4, 4))
}
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))

# Zoom connectors between panels: invisible markers record their device
# position when the figure is drawn; lines are then drawn on top of the figure.
zoom_env <- new.env()
zoom_mark <- function(id, x, y) annotation_custom(grid::grob(id = id, cl = "zoommark"), xmin = x, xmax = x, ymin = y, ymax = y)
drawDetails.zoommark <- function(x, recording) {
  zoom_env[[x$id]] <- unlist(grid::deviceLoc(grid::unit(0.5, "npc"), grid::unit(0.5, "npc"), valueOnly = TRUE))
}

# -----------------------------------------------------------------------------
# 1. Plots with data in both surveys
# -----------------------------------------------------------------------------
clean_dir <- here("pipeline", "01_data_prep", "01_cleaning")
info_file <- file.path(clean_dir, "outputs", "plot_info.csv")
if (!file.exists(info_file)) info_file <- list.files(clean_dir, "plot_info.*\\.csv$", recursive = TRUE, full.names = TRUE)[1]
info <- read_csv(info_file, show_col_types = FALSE)
pick <- function(cand) { x <- intersect(cand, names(info)); if (!length(x)) stop("None of ", paste(cand, collapse = "/"), " in plot_info; columns: ", paste(names(info), collapse = ", ")); x[1] }
id_col <- pick(c("plot_id", "plot", "Plot", "PlotID"))
lat_col <- pick(c("lat", "latitude", "Latitude", "LAT", "lat_dd"))
lon_col <- pick(c("lon", "long", "longitude", "Longitude", "LON", "lon_dd"))
elev_col <- pick(c("elevation", "elev", "Elevation", "elev_m", "ELEV"))

strata_plots <- map_dfr(c("overstory", "understory"), function(st) {
  d <- read_csv(file.path(clean_dir, "outputs", paste0(st, "_long.csv")), show_col_types = FALSE)
  tcol <- intersect(c("period", "time", "survey", "census"), names(d))[1]
  d %>% group_by(plot_id) %>% filter(n_distinct(.data[[tcol]]) == 2) %>% ungroup() %>%
    distinct(plot_id) %>% mutate(stratum = st)
})
plots <- strata_plots %>% group_by(plot_id) %>% summarise(strata = paste(sort(stratum), collapse = " + "), .groups = "drop") %>%
  left_join(info %>% transmute(plot_id = .data[[id_col]], lat = .data[[lat_col]], lon = .data[[lon_col]],
                               elevation = .data[[elev_col]]) %>% distinct(plot_id, .keep_all = TRUE),
            by = "plot_id")
cat("Plots with data in both surveys:\n"); print(count(strata_plots, stratum))
cat("Any stratum:", nrow(plots), "; without coordinates:", sum(is.na(plots$lat) | is.na(plots$lon)), "\n")
cat("Elevation range (m):", paste(round(range(plots$elevation, na.rm = TRUE)), collapse = "-"), "\n")
if (median(plots$lon, na.rm = TRUE) > 0) stop("Longitudes are positive; check the sign in plot_info")

plots_utm <- plots %>% filter(!is.na(lat), !is.na(lon)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% st_transform(crs_utm)
# example plot for panel C: surveyed in both layers, closest to the median elevation
example <- plots_utm %>% filter(strata == "overstory + understory", !is.na(elevation)) %>%
  slice_min(abs(elevation - median(elevation)), n = 1, with_ties = FALSE)
cat("Example plot:", example$plot_id, "(", round(example$elevation), "m )\n")

bb <- st_bbox(plots_utm) + c(-buffer_m, -buffer_m, buffer_m, buffer_m)
ext_utm <- st_as_sfc(bb)
bb_ll <- st_bbox(st_transform(ext_utm, 4326)) + c(-0.05, -0.05, 0.05, 0.05)

boundary_file <- list.files(raw_dir, "boundary.*\\.(shp|gpkg|geojson)$", ignore.case = TRUE, full.names = TRUE)[1]
boundary <- if (!is.na(boundary_file)) st_read(boundary_file, quiet = TRUE) %>% st_transform(crs_utm) else NULL

# -----------------------------------------------------------------------------
# 2. Rasters: CHELSA MAT and MAP, hillshade
# -----------------------------------------------------------------------------
get_chelsa <- function(v) {
  f <- file.path(raw_dir, sprintf("chelsa_%s_1981-2010_window.tif", v))
  if (!file.exists(f)) {
    r <- crop(rast(paste0("/vsicurl/", chelsa_url(v))), ext(bb_ll["xmin"], bb_ll["xmax"], bb_ll["ymin"], bb_ll["ymax"]))
    writeRaster(r, f, overwrite = TRUE)
  }
  rast(f)
}
bio1 <- get_chelsa("bio1")
if (global(bio1, "max", na.rm = TRUE)[[1]] > 100) bio1 <- bio1 * 0.1 - 273.15      # stored as K x 10
bio12 <- get_chelsa("bio12")
if (global(bio12, "max", na.rm = TRUE)[[1]] > 10000) bio12 <- bio12 * 0.1          # stored as mm x 10
to_utm <- function(r, nm) { r <- project(r, crs_utm, method = "bilinear", res = 100) %>% crop(vect(ext_utm)); names(r) <- nm; r }
mat_utm <- to_utm(bio1, "MAT"); map_utm <- to_utm(bio12, "MAP")

dem_file <- file.path(raw_dir, "dem_window.tif")
if (!file.exists(dem_file)) {
  dem <- elevatr::get_elev_raster(locations = st_sf(geometry = st_transform(ext_utm, 4326)), z = 10, clip = "bbox")
  writeRaster(rast(dem), dem_file, overwrite = TRUE)
}
dem <- rast(dem_file) %>% project(crs_utm, method = "bilinear", res = 50) %>% crop(vect(ext_utm))
hill <- shade(terrain(dem, "slope", unit = "radians"), terrain(dem, "aspect", unit = "radians"), angle = 40, direction = 315)
names(hill) <- "hill"

# -----------------------------------------------------------------------------
# 3. Panels A and B: climate maps (A with locator inset)
# -----------------------------------------------------------------------------
climate_map <- function(r, scale_fill, title, legend_title) {
  ggplot() +
    geom_spatraster(data = hill, show.legend = FALSE) +
    scale_fill_gradient(low = "grey15", high = "white", na.value = NA, guide = "none") +
    new_scale_fill() +
    geom_spatraster(data = r, alpha = 0.6) +
    scale_fill +
    { if (!is.null(boundary)) geom_sf(data = boundary, fill = NA, colour = "black", linewidth = 0.5) } +
    geom_sf(data = plots_utm, shape = 21, fill = "black", colour = "white", size = 1.7, stroke = 0.3) +   # black with a white rim, visible on both palettes
    geom_sf(data = example, shape = 22, fill = "black", colour = "white", size = 3.8, stroke = 1) +
    annotation_scale(location = "bl", width_hint = 0.25, height = unit(0.15, "cm"), text_cex = 0.8) +
    annotation_north_arrow(location = "tr", style = north_arrow_minimal(), height = unit(0.9, "cm"), width = unit(0.6, "cm")) +
    coord_sf(xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")], crs = crs_utm, expand = FALSE, datum = NA) +
    labs(title = title, fill = legend_title) +
    theme_map() + theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
                        legend.key.height = unit(0.45, "cm"), legend.key.width = unit(0.35, "cm"))
}
ex_xy <- st_coordinates(example)
pA_main <- climate_map(mat_utm, scale_fill_viridis_c(option = "plasma", end = 0.9, na.value = NA),
                       paste0("A. Mean annual temperature (n = ", nrow(plots_utm), " plots)"), "MAT (°C)") +
  zoom_mark("src", ex_xy[1, "X"], ex_xy[1, "Y"])
pB <- climate_map(map_utm, scale_fill_gradientn(colours = pal_moist, na.value = NA),
                  "B. Mean annual precipitation", "MAP (mm)")

land <- ne_countries(scale = 50, continent = "North America", returnclass = "sf") %>% st_transform(crs_map)
lim_loc <- st_bbox(st_transform(st_as_sfc(st_bbox(c(xmin = -135, xmax = -65, ymin = 22, ymax = 60), crs = 4326)), crs_map))
park_pt <- st_centroid(st_transform(ext_utm, crs_map))
p_loc <- ggplot() +
  geom_sf(data = land, fill = "grey88", colour = "grey55", linewidth = 0.15) +
  geom_sf(data = park_pt, shape = 22, fill = "black", colour = "white", size = 2.6, stroke = 0.6) +
  coord_sf(xlim = lim_loc[c("xmin", "xmax")], ylim = lim_loc[c("ymin", "ymax")], crs = crs_map, expand = FALSE, datum = NA) +
  theme_void() + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
                       panel.background = element_rect(fill = "white", colour = NA))
pA <- pA_main + inset_element(p_loc, left = 0.01, bottom = 0.70, right = 0.36, top = 0.99, align_to = "panel")

# -----------------------------------------------------------------------------
# 4. Panel C: survey design (schematic, arbitrary units)
# -----------------------------------------------------------------------------
# The two surveys differ in how the understory was sampled and how stems were
# recorded, which is why the past understory data are reduced to the present
# sampling effort (Supporting Information S3) and present DBH values are
# assigned to the past DBH classes. The schematic shows the sampled units of
# each survey: in the past, the whole plot for both layers; in the present,
# the whole plot for the overstory, a 50-m2 census at the plot centre for
# understory presence and four 1-m2 quadrats within it, in the cardinal
# directions, for understory cover (quadrats drawn larger than to scale).
set.seed(11)
circle <- function(cx, cy, r, n = 200) tibble(x = cx + r * cos(seq(0, 2 * pi, length.out = n)), y = cy + r * sin(seq(0, 2 * pi, length.out = n)))
blob <- function(cx, cy, r, id) {
  a <- seq(0, 2 * pi, length.out = 60)
  rr <- r * (1 + 0.25 * sin(3 * a + runif(1, 0, 2 * pi)) + 0.12 * sin(5 * a + runif(1, 0, 2 * pi)))
  tibble(x = cx + rr * cos(a), y = cy + rr * sin(a), id = id)
}
square <- function(cx, cy, h, id) tibble(x = cx + c(-h, h, h, -h), y = cy + c(-h, -h, h, h), id = id)
R <- 1.35; cx <- c(past = 2.2, present = 6.4); cy <- 0
r_census <- R * sqrt(50 / 500)                 # 50 m2 within a 500 m2 plot (to scale)
d_quad <- 0.62 * r_census; h_quad <- 0.065      # quadrats inside the census, cardinal directions (not to scale)
in_circle <- function(n, rmax) { t <- runif(n, 0, 2 * pi); u <- rmax * sqrt(runif(n)); tibble(dx = u * cos(t), dy = u * sin(t)) }
trees <- in_circle(14, R - 0.2) %>% mutate(size = runif(14, 1.5, 6.5), tree = row_number())
trees_past <- trees %>% filter(!tree %in% c(3, 9)) %>% mutate(x = cx["past"] + dx, y = cy + dy)
trees_pres <- bind_rows(trees %>% filter(!tree %in% c(5, 12)) %>% mutate(size = pmin(size * 1.15, 7.5)),
                        in_circle(3, R - 0.2) %>% mutate(size = 1.4)) %>% mutate(x = cx["present"] + dx, y = cy + dy)
under <- in_circle(9, R - 0.3)
under_past <- map_dfr(1:9, ~ blob(cx["past"] + under$dx[.x], cy + under$dy[.x], runif(1, 0.15, 0.30), paste0("p", .x)))
under_pres <- map_dfr(1:9, ~ blob(cx["present"] + under$dx[.x], cy + under$dy[.x], runif(1, 0.10, 0.24), paste0("q", .x)))
plot_outline <- bind_rows(circle(cx["past"], cy, R) %>% mutate(id = "past"), circle(cx["present"], cy, R) %>% mutate(id = "present"))
census <- circle(cx["present"], cy, r_census)
quads <- bind_rows(square(cx["present"], cy + d_quad, h_quad, "N"), square(cx["present"] + d_quad, cy, h_quad, "E"),
                   square(cx["present"], cy - d_quad, h_quad, "S"), square(cx["present"] - d_quad, cy, h_quad, "W"))

txt <- function(x, y, label, ...) annotate("text", x = x, y = y, label = label, family = "Arial Narrow", ...)
y_tl <- 2.05
y_area <- -R - 0.25                             # plot area
y_desc <- -R - 0.55                             # what was recorded in each survey
box_C <- c(xmin = 0.85, xmax = 8.05, ymin = -R - 1.75, ymax = y_tl + 0.65)   # extent reached by the zoom lines
pC <- ggplot() +
  zoom_mark("dst1", box_C["xmin"], box_C["ymax"]) + zoom_mark("dst2", box_C["xmax"], box_C["ymax"]) +
  # timeline
  annotate("segment", x = 1.0, xend = 7.8, y = y_tl, yend = y_tl, linewidth = 0.7,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  annotate("segment", x = cx, xend = cx, y = y_tl - 0.12, yend = y_tl + 0.12, linewidth = 0.7) +
  txt(cx, y_tl + 0.35, c("Past survey (~1975)", "Resurvey (2014–2020)"), fontface = "bold", size = 4.3) +
  annotate("segment", x = cx, xend = cx, y = y_tl - 0.15, yend = R + 0.12, linewidth = 0.3, linetype = "22", colour = "grey40") +
  # plots: the past plot is shaded because the understory was recorded over the whole plot
  geom_polygon(data = plot_outline, aes(x, y, group = id, fill = id), colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = c(past = alpha(col_under, 0.22), present = col_plot), guide = "none") +
  geom_polygon(data = under_past, aes(x, y, group = id), fill = col_under, colour = NA, alpha = 0.8) +
  geom_polygon(data = under_pres, aes(x, y, group = id), fill = col_under, colour = NA, alpha = 0.45) +
  # present understory sampling units
  geom_polygon(data = census, aes(x, y), fill = alpha(col_under, 0.35), colour = "black", linewidth = 0.45, linetype = "22") +
  geom_point(data = bind_rows(trees_past, trees_pres), aes(x, y, size = size), shape = 21, fill = col_tree, colour = "white", stroke = 0.3) +
  geom_polygon(data = quads, aes(x, y, group = id), fill = alpha("white", 0.85), colour = "black", linewidth = 0.45) +
  scale_size_identity() +
  txt(cx, y_area, c("250–1000 m² (most 500 m²)", "500 m²"), size = 3.8) +
  txt(cx, y_desc, vjust = 1, size = 3.6, lineheight = 0.95,
      label = c("Overstory: stems tallied in\n10-cm DBH classes\nUnderstory: presence and cover\nover the whole plot",
                "Overstory: DBH measured per stem\nUnderstory: presence in 50-m² census;\ncover in four 1-m² quadrats\nwithin the census")) +
  # legend
  annotate("point", x = 8.55, y = 1.15, shape = 21, size = 5, fill = col_tree, colour = "white") +
  txt(8.85, 1.15, hjust = 0, size = 4, lineheight = 0.95, label = paste0("Overstory stem ", emoji_tree)) +
  annotate("polygon", x = blob(8.55, 0.4, 0.2, "k")$x, y = blob(8.55, 0.4, 0.2, "k")$y, fill = col_under, alpha = 0.8) +
  txt(8.85, 0.4, hjust = 0, size = 4, lineheight = 0.95, label = paste0("Understory species ", emoji_leaf)) +
  geom_polygon(data = circle(8.55, -0.35, 0.2), aes(x, y), fill = alpha(col_under, 0.35), colour = "black", linewidth = 0.45, linetype = "22") +
  txt(8.85, -0.35, hjust = 0, size = 4, label = "50-m² census") +
  geom_polygon(data = square(8.55, -1.1, 0.1, "leg"), aes(x, y), fill = "white", colour = "black", linewidth = 0.45) +
  txt(8.85, -1.1, hjust = 0, size = 4, lineheight = 0.95, label = "1-m² quadrat\n(not to scale)") +
  coord_equal(xlim = c(0.7, 12.0), ylim = c(box_C["ymin"] - 0.05, box_C["ymax"] + 0.05), expand = FALSE, clip = "off") +
  labs(title = "C. Survey design") +
  theme_void() + theme(plot.title = element_text(family = "Arial Narrow", size = 13, face = "bold"),
                       plot.margin = margin(6, 6, 6, 6))

# -----------------------------------------------------------------------------
# 5. Assemble and save
# -----------------------------------------------------------------------------
p <- (pA | pB) / pC + plot_layout(heights = c(1.5, 1.25))

# draw the figure, then the zoom lines from the example plot in A to panel C
draw_fig <- function() {
  print(p)
  s <- zoom_env$src; d1 <- zoom_env$dst1; d2 <- zoom_env$dst2
  if (is.null(s) || is.null(d1) || is.null(d2)) { warning("Zoom markers not found; figure drawn without connectors"); return(invisible()) }
  grid::pushViewport(grid::viewport())
  u <- function(v) grid::unit(v, "in")
  for (d in list(d1, d2)) grid::grid.lines(u(c(s[1], d[1])), u(c(s[2], d[2])), gp = grid::gpar(col = "grey35", lwd = 0.9))
  grid::grid.points(u(s[1]), u(s[2]), pch = 22, size = grid::unit(3.2, "mm"),
                    gp = grid::gpar(fill = "black", col = "white", lwd = 1.2))
  grid::popViewport()
}
out_png <- file.path(fig_dir, "figureS1_study_area.png")
ragg::agg_png(out_png, width = 13, height = 12, units = "in", res = 600, background = "white"); draw_fig(); dev.off()
if (requireNamespace("svglite", quietly = TRUE)) {        # vector version (editable in PowerPoint/Illustrator)
  svglite::svglite(file.path(fig_dir, "figureS1_study_area.svg"), width = 13, height = 12); draw_fig(); dev.off()
}
cat("\nFigure written to", out_png, "\n")
open_fig(out_png)