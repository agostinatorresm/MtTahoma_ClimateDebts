# =============================================================================
# 05d_species_change_elevation_heatmap.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the local computer, after
#           04_hmsc_results.R for the final overstory and understory models)
# Author: Agostina Torres
#
# Change between past and present of each species along the elevational
# gradient, as heat maps (overstory top, understory bottom): species in rows
# (ordered by the mean elevation of the plots where they occur, names coloured
# by thermal affinity, mountain profile as in 05b), elevation in columns,
# colour = predicted change at that elevation (fixed effects, ~ time * elevation;
# 04_hmsc_results.R, delta_elevation_species.csv). Cells are shown only within
# the elevational range of the plots where each species occurs. Changes with
# f >= 0.90 in full colour, 0.80 <= f < 0.90 faded, f < 0.80 in light grey. Elevation axis and
# colour limits are shared between layers (one legend per response).
#
# Occurrence: change in occurrence probability (understory model fitted to
# reduced past data). Abundance: change in log10 relative abundance; for the
# understory it is shown relative to the mean change across the species present
# at each elevation, centered within each posterior sample
# (25_species_change_elevation_centered.R), to remove the common increase in
# relative cover (as in 05b).
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/<stratum>_..._rel_res_delta_elevation_species.csv (04)
#          pipeline/02_analysis/02_hmsc/outputs/centered/understory_abundance_change_elevation_centered.csv (25)
#          pipeline/01_data_prep/01_cleaning/outputs/plot_info.csv, <stratum>_long.csv
# Outputs: pipeline/02_analysis/02_hmsc/figures/m1_final/ms1__fig2b_species_change_elevation_heatmap_m1_final.png
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)
library(patchwork)

tag        <- "m1_final"
fig_prefix <- "ms1__"
model_files <- c(
  overstory  = "overstory_strict_medium_m1_time_elev_nf2_nsp_nos_rel_res_delta_elevation_species.csv",
  understory = "understory_strict_long_m1_time_elev_nsp_nos_redpool_rel_res_delta_elevation_species.csv")
centered_strata <- "understory"     # abundance centred across species only here (as in 05b)
strata_titles <- c(overstory = "\U0001F332 Overstory", understory = "\U0001F33F Understory")
sig_level  <- 0.90   # full colour
weak_level <- 0.80   # faded colour (0.80 <= f < 0.90); below, light grey
col_affinity <- viridisLite::plasma(256, end = 0.80, direction = 1)   # as in 05b
mat_limits   <- c(0, 12.5)
col_change   <- c("#2B5C8A", "white", "#C4561D")                      # decrease, none, increase

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))

theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text = element_text(colour = "black", size = 13),
          axis.title.x = element_text(colour = "black", size = 16, hjust = 1),
          legend.text = element_text(colour = "black", size = 12),
          legend.title = element_text(colour = "black", size = 13),
          plot.title = element_text(colour = "black", size = 17, face = "bold"),
          plot.margin = margin(12, 10, 10, 10))
}
sp_label <- function(x) str_replace_all(x, "_", " ")

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
fig_dir    <- file.path(module_dir, "figures", tag)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cleaning_dir <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)

# -----------------------------------------------------------------------------
# 1. Inputs
# -----------------------------------------------------------------------------
missing <- model_files[!file.exists(file.path(module_dir, "outputs", model_files))]
if (length(missing) > 0) stop("Files not found (edit model_files):\n", paste(missing, collapse = "\n"))

plot_info <- read_csv(cleaning_dir("outputs", "plot_info.csv"), show_col_types = FALSE)
plot_key  <- intersect(c("plot_id", "plot", "Plot", "PlotID"), names(plot_info))[1]
elev_col  <- intersect(c("elevation", "elevation_m", "elev", "elev_m"), names(plot_info))[1]
plot_elev <- plot_info %>%
  transmute(plot_id = as.character(.data[[plot_key]]), elevation = as.numeric(.data[[elev_col]])) %>%
  filter(!is.na(elevation)) %>% distinct(plot_id, .keep_all = TRUE)

read_stratum <- function(st) {
  d <- read_csv(file.path(module_dir, "outputs", model_files[[st]]), show_col_types = FALSE) %>%
    filter(elev_sc >= elev_min, elev_sc <= elev_max) %>%      # within each species' range
    mutate(response = recode(response, abundance = "log10_abundance"), stratum = st)
  stopifnot("Expected occurrence and log10_abundance" = all(c("occurrence", "log10_abundance") %in% d$response))
  
  # Understory abundance: change relative to the mean change across the species
  # present at each elevation, centered within each posterior sample
  # (25_species_change_elevation_centered.R); other panels on the original scale
  if (st %in% centered_strata) {
    cen <- read_csv(file.path(module_dir, "outputs", "centered", "understory_abundance_change_elevation_centered.csv"),
                    show_col_types = FALSE) %>%
      filter(version == "centered") %>% transmute(species, elev_key = round(elev_sc, 6), mean_c = mean, support_c = support)
    d <- d %>%
      mutate(elev_key = round(elev_sc, 6)) %>%
      left_join(cen, by = c("species", "elev_key")) %>%
      mutate(mean    = if_else(response == "log10_abundance", mean_c, mean),
             support = if_else(response == "log10_abundance", support_c, support)) %>%
      select(-mean_c, -support_c, -elev_key)
    stopifnot("Centered change missing for some cells" = !any(is.na(d$mean)))
  }
  d <- d %>%
    mutate(value = mean,
           f = pmax(support, 1 - support),
           support_class = case_when(f >= sig_level ~ "strong", f >= weak_level ~ "weak", TRUE ~ "none"))
  
  # Species order: mean elevation of the plots where each species occurs (as in 05)
  species_elev <- read_csv(cleaning_dir("outputs", paste0(st, "_long.csv")), show_col_types = FALSE) %>%
    filter(occurrence > 0) %>% mutate(plot_id = as.character(plot_id)) %>%
    left_join(plot_elev, by = "plot_id") %>%
    group_by(species) %>% summarise(mean_elevation = mean(elevation, na.rm = TRUE), .groups = "drop") %>%
    mutate(species = str_replace_all(species, " ", "_"))
  sp_info <- d %>% distinct(species, MAT_mean) %>%
    left_join(species_elev, by = "species") %>%
    arrange(mean_elevation) %>%
    mutate(y = row_number(), label = sp_label(species))
  list(d = left_join(d, select(sp_info, species, y), by = "species"), sp = sp_info)
}
strata <- set_names(names(model_files)) %>% map(read_stratum)
all_d  <- map_dfr(strata, "d")
# Each layer has its own elevation grid, so the tile width is computed per layer
steps    <- map_dbl(strata, ~ median(diff(sort(unique(.x$d$elevation)))))
elev_lim <- range(map2(strata, steps, ~ range(.x$d$elevation) + c(-1, 1) * .y / 2))
fill_lim <- all_d %>% group_by(response) %>% summarise(m = max(abs(value))) %>% deframe()

# -----------------------------------------------------------------------------
# 2. Mountain profile and species names (as in 05b)
# -----------------------------------------------------------------------------
displace <- function(x, y, iters, sd, decay = 0.55) {
  for (i in seq_len(iters)) {
    k  <- length(x)
    mx <- (x[-k] + x[-1]) / 2
    my <- (y[-k] + y[-1]) / 2 + rnorm(k - 1, 0, sd)
    x  <- c(rbind(x[-k], mx), x[k]); y <- c(rbind(y[-k], my), y[k])
    sd <- sd * decay
  }
  tibble(x = x, y = y)
}
set.seed(7)
ridge <- displace(c(0, 0.16, 0.30, 0.42, 0.52, 0.60, 0.70, 0.80, 0.90, 1.00),
                  c(0.58, 0.66, 0.80, 0.90, 1.00, 0.88, 0.70, 0.48, 0.24, 0.00), iters = 6, sd = 0.045)
rocks <- imap_dfr(list(list(x = c(0.52, 0.50, 0.48), y = c(0.99, 0.93, 0.88)),
                       list(x = c(0.56, 0.58, 0.57), y = c(0.94, 0.89, 0.84))),
                  ~ displace(.x$x, .x$y, iters = 3, sd = 0.015) %>% mutate(g = .y))

names_panel <- function(st) {
  sp_info <- strata[[st]]$sp; n <- nrow(sp_info)
  ggplot() +
    geom_path(data = mutate(filter(ridge, x >= 0.30), y = 0.5 + y * n), aes(x, y), colour = "grey50",
              alpha = 0.6, linewidth = 0.45, lineend = "round", linejoin = "round") +
    geom_path(data = mutate(rocks, y = 0.5 + y * n), aes(x, y, group = g), colour = "grey50", alpha = 0.6, linewidth = 0.3) +
    geom_text(data = sp_info, aes(x = 0.98, y = y, label = label, colour = MAT_mean),
              hjust = 1, fontface = "italic", size = 3.3) +
    scale_colour_gradientn(colours = col_affinity, limits = mat_limits, oob = scales::squish, guide = "none") +
    coord_cartesian(xlim = c(0.30, 1), ylim = c(0.5, n + 0.5), expand = FALSE, clip = "off") +
    labs(title = strata_titles[[st]]) +
    theme_void() +
    theme(plot.title = element_text(colour = "black", size = 19, face = "bold", margin = margin(b = 14)),
          plot.margin = margin(12, 10, 10, 10))
}

# -----------------------------------------------------------------------------
# 3. Heat maps (shared elevation axis and colour limits per response)
# -----------------------------------------------------------------------------
fold_labels <- function(x) paste0("×", signif(10^x, 2))
fill_labs <- c(occurrence = "Change in occurrence probability",
               log10_abundance = "Change in relative abundance")
heat <- function(st, resp, show_x) {
  dd <- filter(strata[[st]]$d, response == resp); n <- nrow(strata[[st]]$sp)
  p <- ggplot(mapping = aes(x = elevation, y = y)) +
    # f >= 0.90 full colour, 0.80 <= f < 0.90 faded, f < 0.80 light grey
    geom_tile(data = filter(dd, support_class == "none"), fill = "grey93", width = steps[[st]], height = 0.92) +
    geom_tile(data = filter(dd, support_class != "none"), aes(fill = value, alpha = support_class),
              width = steps[[st]], height = 0.92) +
    scale_alpha_manual(values = c(strong = 1, weak = 0.45), guide = "none") +
    scale_fill_gradientn(colours = col_change, limits = fill_lim[[resp]] * c(-1, 1),
                         labels = if (resp == "log10_abundance") fold_labels else waiver(),
                         name = fill_labs[[resp]],
                         guide = guide_colourbar(barwidth = unit(5, "cm"), title.position = "top")) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = elev_lim, ylim = c(0.5, n + 0.5), expand = FALSE) +
    labs(x = if (show_x) "Elevation (m)" else NULL, y = NULL,
         title = if (st == names(model_files)[1]) c(occurrence = "Occurrence", log10_abundance = "Abundance")[[resp]] else NULL) +
    theme_mora() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

rows <- map(names(model_files), function(st) {
  last <- st == tail(names(model_files), 1)
  names_panel(st) + heat(st, "occurrence", last) + heat(st, "log10_abundance", last) +
    plot_layout(widths = c(0.9, 1, 1))
})
n_sp <- map_int(strata, ~ nrow(.x$sp))
# One legend per response at the top (colour limits are shared between strata)
fig <- wrap_plots(rows, ncol = 1, heights = n_sp + 3) +   # +3 rows for titles
  plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.justification = "right")
fig_file <- file.path(fig_dir, paste0(fig_prefix, "fig2b_species_change_elevation_heatmap_", tag, ".png"))
ggsave(fig_file, fig,
       device = "png", dpi = 300, width = 10, height = 12, scale = 1.3, bg = "white")
cat("Figure written to", fig_dir, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
