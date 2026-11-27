# ------------------------------------------------------------------------------
# 07_species_replacement_chord.R
# Module: pipeline/02_analysis/01_thermophilization/01_overstory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: show which overstory species replaced which between the past
# (~1975) and present (~2015) surveys, as a chord diagram.
#
# Method: in each plot, the changes in relative basal area (p_present -
# p_past, same pool and weights as 01 and 05) sum to zero, so the share lost
# by declining species equals the share gained by increasing ones. The loss
# of each declining species i is allocated to the increasing species j in
# proportion to their gains:
#     flow(i -> j) = loss_i * gain_j / total gain in the plot
# Flows are summed over plots and divided by the number of plots, so they are
# expressed as the mean share of plot basal area transferred from i to j.
# This is an allocation rule, not a record of individual trees: it assumes
# that, within a plot, every gaining species took over the losses of all
# declining species in proportion to its gain.
#
# Sectors are ordered by thermal affinity (MAT_mean, coldest first) and
# coloured with the same plasma scale as Figure 3, so flows from cold to warm
# species run clockwise. A second version orders and colours species by
# moisture affinity (MAP_mean). Links take the colour of the species that
# lost basal area; flows smaller than `min_flow_share` of the total are
# drawn but not labelled.
#
# Gross flows run in both directions between two species when one replaced
# the other in some plots and the reverse happened in others. The net version
# keeps, for each pair, only the difference between the two directions, so
# each chord shows the dominant direction of replacement.
#
# Output:
#   - outputs/overstory_species_replacement_flows.csv (from, to, flow, net)
#   - figures/overstory_species_replacement_chord_{MAT,MAP}.png (gross)
#   - figures/overstory_species_replacement_chord_{MAT,MAP}_net.png (net)
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(circlize)      # install.packages("circlize") if missing
library(viridisLite)

# 0a. Module directories --------------------------------------------------------
this_module_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)
cleaning_dir     <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
affinity_dir     <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
species_pool_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)
dir.create(this_module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

# Palettes and limits shared with Figure 3 (02_hmsc/05b_species_change_centered_figure.R)
COL_MAT    <- viridisLite::plasma(256, end = 0.80)
MAT_LIMITS <- c(0, 12.5)
COL_MAP    <- c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D")
MAP_LIMITS <- c(800, 2700)   # must equal map_lims in 05b_species_change_centered_figure.R
min_flow_share <- 0.01       # flows below 1% of the total are drawn thinner and lighter

# 1. Load (as in 05_species_contribution.R) -------------------------------------
overstory_long <- read_csv(cleaning_dir("outputs", "overstory_long.csv"), show_col_types = FALSE)
species_affinity <- read_csv(affinity_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)
overstory_pool <- read_csv(species_pool_dir("outputs", "overstory_species_pool_broad_full_clim_cov.csv"),
                           show_col_types = FALSE) %>% pull(species)

# 2. Change in relative basal area per plot and species --------------------------
dp <- overstory_long %>%
  filter(occurrence == 1) %>%
  select(plot_id, period, species, weight = basal_area_spp_perc) %>%
  complete(nesting(plot_id, species), period, fill = list(weight = 0)) %>%
  filter(species %in% species_affinity$species, species %in% overstory_pool) %>%
  pivot_wider(names_from = period, values_from = weight, values_fill = 0) %>%
  group_by(plot_id) %>%
  filter(sum(past) > 0, sum(present) > 0) %>%
  mutate(dp = present / sum(present) - past / sum(past)) %>%
  ungroup()
n_plots <- n_distinct(dp$plot_id)
cat("Plots with basal area in both surveys:", n_plots, "\n")

# 3. Flows from declining to increasing species within each plot ------------------
flows <- dp %>%
  group_by(plot_id) %>%
  group_modify(~ {
    loss <- .x %>% filter(dp < 0) %>% transmute(from = species, loss = -dp)
    gain <- .x %>% filter(dp > 0) %>% transmute(to = species, gain = dp)
    if (nrow(loss) == 0 || nrow(gain) == 0) return(tibble(from = character(), to = character(), flow = numeric()))
    tidyr::crossing(loss, gain) %>% mutate(flow = loss * gain / sum(gain[!duplicated(to)])) %>% select(from, to, flow)
  }) %>%
  ungroup() %>%
  group_by(from, to) %>%
  summarise(flow = sum(flow) / n_plots, .groups = "drop") %>%
  arrange(desc(flow))
# Net flow per pair: flow(i -> j) - flow(j -> i), kept where positive
flows_net <- flows %>%
  left_join(flows %>% select(from = to, to = from, back = flow), by = c("from", "to")) %>%
  mutate(flow = flow - coalesce(back, 0)) %>%
  filter(flow > 0) %>% select(from, to, flow) %>% arrange(desc(flow))
write_csv(flows %>% left_join(flows_net %>% rename(net = flow), by = c("from", "to")),
          this_module_dir("outputs", "overstory_species_replacement_flows.csv"))

# Check: the total flow equals the mean share of basal area lost per plot
cat("Total flow:", signif(sum(flows$flow), 3), "| mean loss per plot:",
    signif(sum(pmax(-dp$dp, 0)) / n_plots, 3), "\n")

# Direction of replacement relative to the affinities of the two species
dir_sum <- flows %>%
  left_join(species_affinity, by = c("from" = "species")) %>%
  left_join(species_affinity, by = c("to" = "species"), suffix = c("_from", "_to")) %>%
  summarise(to_warmer = sum(flow[MAT_mean_to > MAT_mean_from]) / sum(flow),
            to_wetter = sum(flow[MAP_mean_to > MAP_mean_from]) / sum(flow))
cat("Share of the replaced basal area that went to a warmer-affiliated species:", round(dir_sum$to_warmer, 2),
    "| to a wetter-affiliated species:", round(dir_sum$to_wetter, 2), "\n")
cat("\nLargest flows (mean share of plot basal area):\n"); print(head(flows, 12))
dir_net <- flows_net %>%
  left_join(species_affinity, by = c("from" = "species")) %>%
  left_join(species_affinity, by = c("to" = "species"), suffix = c("_from", "_to")) %>%
  summarise(to_warmer = sum(flow[MAT_mean_to > MAT_mean_from]) / sum(flow),
            to_wetter = sum(flow[MAP_mean_to > MAP_mean_from]) / sum(flow))
cat("\nNet flows: share to a warmer-affiliated species:", round(dir_net$to_warmer, 2),
    "| to a wetter-affiliated species:", round(dir_net$to_wetter, 2), "\n")
cat("Largest net flows:\n"); print(head(flows_net, 12))

# 4. Chord diagram -----------------------------------------------------------------
short_name <- function(x) sub("^(\\w)\\w*\\s", "\\1. ", x)   # "Abies amabilis" -> "A. amabilis"

draw_chord <- function(flows, affinity_col, palette, limits, legend_title, file) {
  sp <- species_affinity %>% filter(species %in% c(flows$from, flows$to)) %>% arrange(.data[[affinity_col]])
  col_fun <- circlize::colorRamp2(seq(limits[1], limits[2], length.out = length(palette)), palette)
  grid_col <- setNames(col_fun(pmin(pmax(sp[[affinity_col]], limits[1]), limits[2])), sp$species)
  small <- flows$flow < min_flow_share * sum(flows$flow)
  link_col <- adjustcolor(grid_col[flows$from], alpha.f = 0.65)
  link_col[small] <- adjustcolor(grid_col[flows$from][small], alpha.f = 0.25)

  png(file, width = 9, height = 9, units = "in", res = 300, bg = "white")
  par(mar = c(1, 1, 1, 1), family = "sans")
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 3, track.margin = c(0.01, 0.01), points.overflow.warning = FALSE)
  chordDiagram(flows %>% as.data.frame(), order = sp$species, grid.col = grid_col, col = link_col,
               directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow",
               diffHeight = -mm_h(2), annotationTrack = "grid", preAllocateTracks = list(track.height = 0.18),
               link.sort = TRUE, link.largest.ontop = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], short_name(CELL_META$sector.index),
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.95, font = 3)
  }, bg.border = NA)
  # Colour bar for the affinity scale
  # Colour bar in the top-left corner, clear of the species labels
  xs <- seq(-1.08, -0.62, length.out = 100)
  rect(xs[-100], 1.02, xs[-1], 1.06, col = col_fun(seq(limits[1], limits[2], length.out = 99)), border = NA, xpd = NA)
  text(c(xs[1], xs[100]), 0.98, limits, cex = 0.8, xpd = NA)
  text(mean(xs), 1.10, legend_title, cex = 0.85, xpd = NA)
  circos.clear()
  dev.off()
  file
}

f_mat <- draw_chord(flows, "MAT_mean", COL_MAT, MAT_LIMITS, "Thermal affinity (MAT, °C)",
                    this_module_dir("figures", "overstory_species_replacement_chord_MAT.png"))
f_map <- draw_chord(flows, "MAP_mean", COL_MAP, MAP_LIMITS, "Moisture affinity (MAP, mm)",
                    this_module_dir("figures", "overstory_species_replacement_chord_MAP.png"))
f_mat_net <- draw_chord(flows_net, "MAT_mean", COL_MAT, MAT_LIMITS, "Thermal affinity (MAT, °C)",
                        this_module_dir("figures", "overstory_species_replacement_chord_MAT_net.png"))
f_map_net <- draw_chord(flows_net, "MAP_mean", COL_MAP, MAP_LIMITS, "Moisture affinity (MAP, mm)",
                        this_module_dir("figures", "overstory_species_replacement_chord_MAP_net.png"))

# 5. Open the figures -----------------------------------------------------------------
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(f_mat_net); open_fig(f_map_net)
