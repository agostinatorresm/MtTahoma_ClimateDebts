# =============================================================================
# 05g_species_change_range_margins.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the laptop, after
#           04_hmsc_results.R and 25_species_change_elevation_centered.R)
# Author: Agostina Torres
#
# Change between past and present of each species along its own elevational
# range, to evaluate upslope shifts relative to each species' distribution
# rather than to absolute elevation. An upslope shift is expected as losses at
# the lower margin of a species' range and gains (or no loss) at its upper
# margin, whatever the species' climatic affinity.
#
# Elevational range of each species: quantiles range_q of the elevations of the
# plots where it occurs (both surveys combined); the outer quantiles are used
# instead of the minimum and maximum so that a few outlying plots do not define
# the margins. Elevation is rescaled within each species' range
# (0 = lower margin, 1 = upper margin).
#
# Change at each elevation as in 05d (fixed effects, ~ time * elevation;
# delta_elevation_species.csv from 04); understory abundance is centered on the
# mean change across species at each elevation (25). Because the model is
# linear in elevation, the contrast between margins reflects each species'
# time x elevation coefficient evaluated over its own range.
#
# Classification of each species and response from the change at the lowest and
# highest grid cell within its range (supported when f >= 0.90):
#   Upslope    : supported loss at the lower margin without loss at the upper
#                margin, or supported gain at the upper margin without gain at
#                the lower margin
#   Downslope  : the reverse
#   Decline / Increase across range : supported change of the same sign at both margins
#   No supported change at margins
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/<stratum>_..._rel_res_delta_elevation_species.csv (04)
#          pipeline/02_analysis/02_hmsc/outputs/centered/understory_abundance_change_elevation_centered.csv (25)
#          pipeline/01_data_prep/01_cleaning/outputs/plot_info.csv, <stratum>_long.csv
# Outputs: pipeline/02_analysis/02_hmsc/outputs/species_change_range_margins.csv
#          pipeline/02_analysis/02_hmsc/figures/m1_final/ms1__figS6_4_species_change_range_margins_m1_final.png
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
centered_strata <- "understory"     # abundance centred across species only here (as in 05b and 05d)
strata_titles <- c(overstory = "\U0001F332 Overstory", understory = "\U0001F33F Understory")   # as in the other HMSC figures
range_q    <- c(0.10, 0.90)         # quantiles of occupied-plot elevations defining each species' range
sig_level  <- 0.90                  # full colour and classification
weak_level <- 0.80                  # faded colour (0.80 <= f < 0.90); below, light grey
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
          axis.title.x = element_text(colour = "black", size = 15, hjust = 0.5),
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
    filter(elev_sc >= elev_min, elev_sc <= elev_max) %>%
    mutate(response = recode(response, abundance = "log10_abundance"), stratum = st)
  
  # Understory abundance: centered change (25), as in 05d
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
  
  # Elevational range of each species from the plots where it occurs
  sp_range <- read_csv(cleaning_dir("outputs", paste0(st, "_long.csv")), show_col_types = FALSE) %>%
    filter(occurrence > 0) %>% mutate(plot_id = as.character(plot_id)) %>%
    left_join(plot_elev, by = "plot_id") %>%
    group_by(species) %>%
    summarise(mean_elevation = mean(elevation, na.rm = TRUE),
              range_lo = quantile(elevation, range_q[1], na.rm = TRUE),
              range_hi = quantile(elevation, range_q[2], na.rm = TRUE), .groups = "drop") %>%
    mutate(species = str_replace_all(species, " ", "_"))
  
  step <- median(diff(sort(unique(d$elevation))))
  d <- d %>%
    left_join(sp_range, by = "species") %>%
    group_by(species, response) %>%
    # keep the cells within the quantile range; if fewer than two remain, keep all cells
    mutate(in_q = elevation >= range_lo - step / 2 & elevation <= range_hi + step / 2,
           in_q = if (sum(in_q) >= 2) in_q else TRUE) %>%
    filter(in_q) %>%
    mutate(lo = min(elevation), hi = max(elevation),
           rel_pos = (elevation - lo) / pmax(hi - lo, step),
           width   = step / pmax(hi - lo, step)) %>%
    ungroup() %>%
    mutate(value = mean,
           f = pmax(support, 1 - support),
           support_class = case_when(f >= sig_level ~ "strong", f >= weak_level ~ "weak", TRUE ~ "none"))
  
  sp_info <- d %>% distinct(species, MAT_mean, mean_elevation, lo, hi) %>%
    group_by(species) %>% slice(1) %>% ungroup() %>%
    arrange(mean_elevation) %>%
    mutate(y = row_number(), label = sp_label(species),
           range_label = paste0(round(lo, -1), "–", round(hi, -1), " m"))
  list(d = left_join(d, select(sp_info, species, y), by = "species"), sp = sp_info)
}
strata <- set_names(names(model_files)) %>% map(read_stratum)
all_d  <- map_dfr(strata, "d")
fill_lim <- all_d %>% group_by(response) %>% summarise(m = max(abs(value))) %>% deframe()

# -----------------------------------------------------------------------------
# 2. Change at the lower and upper margins and classification
# -----------------------------------------------------------------------------
margins <- all_d %>%
  group_by(stratum, species, response) %>%
  summarise(MAT_mean = first(MAT_mean),
            range_lo_m = first(lo), range_hi_m = first(hi),
            lower_change = value[which.min(elevation)], lower_f = f[which.min(elevation)],
            upper_change = value[which.max(elevation)], upper_f = f[which.max(elevation)],
            .groups = "drop") %>%
  mutate(lo_dec = lower_change < 0 & lower_f >= sig_level, lo_inc = lower_change > 0 & lower_f >= sig_level,
         up_dec = upper_change < 0 & upper_f >= sig_level, up_inc = upper_change > 0 & upper_f >= sig_level,
         pattern = case_when(
           lo_dec & !up_dec ~ "Upslope (lower-margin loss)",
           up_inc & !lo_inc ~ "Upslope (upper-margin gain)",
           up_dec & !lo_dec ~ "Downslope (upper-margin loss)",
           lo_inc & !up_inc ~ "Downslope (lower-margin gain)",
           lo_dec & up_dec  ~ "Decline across range",
           lo_inc & up_inc  ~ "Increase across range",
           TRUE             ~ "No supported change at margins")) %>%
  select(-lo_dec, -lo_inc, -up_dec, -up_inc) %>%
  mutate(across(c(lower_change, upper_change), ~ round(.x, 3)),
         across(c(lower_f, upper_f), ~ round(.x, 2))) %>%
  arrange(stratum, response, pattern, species)

write_csv(margins, file.path(module_dir, "outputs", "species_change_range_margins.csv"))
cat("\nPatterns by stratum and response:\n")
print(count(margins, stratum, response, pattern), n = Inf)
cat("\nSpecies with upslope or downslope patterns:\n")
print(filter(margins, str_detect(pattern, "slope")) %>%
        select(stratum, response, species, pattern, lower_change, lower_f, upper_change, upper_f), n = Inf)

# -----------------------------------------------------------------------------
# 3. Figure: change along each species' range (0 = lower margin, 1 = upper margin)
# -----------------------------------------------------------------------------
# Mountain profile behind the species names, as in Figure S6.2 (05b, 05d):
# species are ordered by mean elevation of occurrence, so the profile marks the
# ordering from low (bottom) to high (top) elevation
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
range_panel <- function(st) {
  sp_info <- strata[[st]]$sp; n <- nrow(sp_info)
  ggplot(sp_info) +
    geom_text(aes(x = 0, y = y, label = range_label), hjust = 0, size = 3, colour = "grey30") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.5, n + 0.5), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(12, 0, 10, 0))
}

fold_labels <- function(x) paste0("×", signif(10^x, 2))
fill_labs <- c(occurrence = "Change in occurrence probability",
               log10_abundance = "Change in relative abundance")
heat <- function(st, resp, show_x) {
  dd <- filter(strata[[st]]$d, response == resp); n <- nrow(strata[[st]]$sp)
  p <- ggplot(mapping = aes(x = rel_pos, y = y, width = width)) +
    geom_tile(data = filter(dd, support_class == "none"), fill = "grey93", height = 0.92) +
    geom_tile(data = filter(dd, support_class != "none"), aes(fill = value, alpha = support_class), height = 0.92) +
    scale_alpha_manual(values = c(strong = 1, weak = 0.45), guide = "none") +
    scale_fill_gradientn(colours = col_change, limits = fill_lim[[resp]] * c(-1, 1),
                         labels = if (resp == "log10_abundance") fold_labels else waiver(),
                         name = fill_labs[[resp]],
                         guide = guide_colourbar(barwidth = unit(5, "cm"), title.position = "top")) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("Lower\nmargin", "", "Upper\nmargin")) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(0.5, n + 0.5), expand = FALSE) +
    labs(x = if (show_x) "Position within the species' elevational range" else NULL, y = NULL,
         title = if (st == names(model_files)[1]) c(occurrence = "Occurrence", log10_abundance = "Abundance")[[resp]] else NULL) +
    theme_mora() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

rows <- map(names(model_files), function(st) {
  last <- st == tail(names(model_files), 1)
  names_panel(st) + heat(st, "occurrence", last) + heat(st, "log10_abundance", last) + range_panel(st) +
    plot_layout(widths = c(0.8, 1, 1, 0.3))
})
n_sp <- map_int(strata, ~ nrow(.x$sp))
fig <- wrap_plots(rows, ncol = 1, heights = n_sp + 3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.justification = "right")
fig_file <- file.path(fig_dir, paste0(fig_prefix, "figS6_4_species_change_range_margins_", tag, ".png"))
ggsave(fig_file, fig, device = "png", dpi = 300, width = 10, height = 12, scale = 1.3, bg = "white")
cat("Figure written to", fig_dir, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)