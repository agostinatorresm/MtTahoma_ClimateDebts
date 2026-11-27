# ------------------------------------------------------------------------------
# 08_species_replacement_chord.R
# Module: pipeline/02_analysis/01_thermophilization/02_understory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: understory counterpart of
# 01_overstory/code/07_species_replacement_chord.R: which species (or groups
# of species) replaced which between the past (~1975) and present (~2015)
# surveys, as a chord diagram.
#
# Method: as in the overstory script, with relative cover (same pool and
# weights as 01 and 05). In each plot the changes in relative cover sum to
# zero; the loss of each declining species is allocated to the increasing
# species in proportion to their gains, and flows are summed over plots and
# divided by the number of plots (mean share of plot cover transferred).
# This is an allocation rule, not a record of individual plants. Net flows
# keep, for each pair, only the difference between the two directions.
#
# With over a hundred species, flows are aggregated into groups before
# drawing (both versions are produced):
#   "top_species": the `n_top` species with the largest total flow, and the
#                  remaining species grouped by life form ("Other forbs", ...)
#   "lifeform_thermal": life form x thermal class (cold, intermediate, warm:
#                  terciles of MAT_mean among the pool species)
# Groups are ordered and coloured by their mean thermal (or moisture)
# affinity, weighted by flow, with the palettes of Figure 3.
#
# Caution: present relative cover is inflated because fewer species had cover
# per plot in the present survey (Table S3.1), so gross gains and losses are
# larger than the true turnover; the net version is less affected.
#
# Output:
#   - outputs/understory_species_replacement_flows.csv (species-level flows)
#   - figures/understory_replacement_chord_<grouping>_{MAT,MAP}_net.png, for
#     both groupings; net flows below `draw_min_share` of the total are not
#     drawn. Sedges (one species) are excluded, as in Figure S4.5.
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(circlize)      # install.packages("circlize") if missing
library(viridisLite)

# 0a. Module directories --------------------------------------------------------
this_module_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "02_understory", ...)
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
groupings <- c("top_species", "lifeform_thermal")   # both versions are drawn
n_top    <- 15
draw_min_share <- 0.015      # net flows below 1.5% of the total are not drawn (they are kept in the csv)
LIFEFORMS_EXCLUDED <- c("Sedge")   # one species in five plots, as in Figure S4.5

# 1. Load (as in 05_species_contribution.R) -------------------------------------
understory_long <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE)
spp_info <- read_csv(cleaning_dir("outputs", "spp_info.csv"), show_col_types = FALSE) %>% select(species, life_form)
species_affinity <- read_csv(affinity_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)
understory_pool <- read_csv(species_pool_dir("outputs", "understory_species_pool_broad_full_clim_cov.csv"),
                            show_col_types = FALSE) %>% pull(species)

# 2. Change in relative cover per plot and species ---------------------------------
# Species present only in the 50-m2 census have no cover and weight 0 (as in 05)
dp <- understory_long %>%
  filter(occurrence == 1) %>%
  select(plot_id, period, species, weight = cover) %>%
  mutate(weight = coalesce(weight, 0)) %>%
  complete(nesting(plot_id, species), period, fill = list(weight = 0)) %>%
  filter(species %in% species_affinity$species, species %in% understory_pool,
         !species %in% spp_info$species[spp_info$life_form %in% LIFEFORMS_EXCLUDED]) %>%
  pivot_wider(names_from = period, values_from = weight, values_fill = 0) %>%
  group_by(plot_id) %>%
  filter(sum(past) > 0, sum(present) > 0) %>%
  mutate(dp = present / sum(present) - past / sum(past)) %>%
  ungroup()
n_plots <- n_distinct(dp$plot_id)
cat("Plots with cover in both surveys:", n_plots, "\n")

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
          this_module_dir("outputs", "understory_species_replacement_flows.csv"))

# Check: the total flow equals the mean share of cover lost per plot
cat("Total flow:", signif(sum(flows$flow), 3), "| mean loss per plot:",
    signif(sum(pmax(-dp$dp, 0)) / n_plots, 3), "\n")

# Direction of replacement relative to the affinities of the two species
dir_sum <- flows %>%
  left_join(species_affinity, by = c("from" = "species")) %>%
  left_join(species_affinity, by = c("to" = "species"), suffix = c("_from", "_to")) %>%
  summarise(to_warmer = sum(flow[MAT_mean_to > MAT_mean_from]) / sum(flow),
            to_wetter = sum(flow[MAP_mean_to > MAP_mean_from]) / sum(flow))
cat("Share of the replaced cover that went to a warmer-affiliated species:", round(dir_sum$to_warmer, 2),
    "| to a wetter-affiliated species:", round(dir_sum$to_wetter, 2), "\n")
cat("\nLargest flows (mean share of plot cover):\n"); print(head(flows, 12))
dir_net <- flows_net %>%
  left_join(species_affinity, by = c("from" = "species")) %>%
  left_join(species_affinity, by = c("to" = "species"), suffix = c("_from", "_to")) %>%
  summarise(to_warmer = sum(flow[MAT_mean_to > MAT_mean_from]) / sum(flow),
            to_wetter = sum(flow[MAP_mean_to > MAP_mean_from]) / sum(flow))
cat("\nNet flows: share to a warmer-affiliated species:", round(dir_net$to_warmer, 2),
    "| to a wetter-affiliated species:", round(dir_net$to_wetter, 2), "\n")
cat("Largest net flows:\n"); print(head(flows_net, 12))

# 4. Groups -----------------------------------------------------------------------
sp_info <- species_affinity %>% left_join(spp_info, by = "species") %>%
  mutate(life_form = coalesce(life_form, "Other"),
         thermal = cut(MAT_mean, quantile(MAT_mean[species %in% dp$species], c(0, 1/3, 2/3, 1)),
                       labels = c("cold", "intermediate", "warm"), include.lowest = TRUE))
for (grouping in groupings) {
total_by_sp <- bind_rows(flows %>% select(sp = from, flow), flows %>% select(sp = to, flow)) %>%
  group_by(sp) %>% summarise(total = sum(flow), .groups = "drop") %>% arrange(desc(total))
plural <- function(x) ifelse(grepl("s$", x), paste0(x, "es"), paste0(x, "s"))
group_of <- function(x) {
  i <- match(x, sp_info$species)
  if (grouping == "top_species") ifelse(x %in% head(total_by_sp$sp, n_top), x, paste("Other", tolower(plural(sp_info$life_form[i]))))
  else paste(sp_info$life_form[i], sp_info$thermal[i], sep = ", ")
}
group_flows <- function(fl) {
  fl %>% mutate(from_g = group_of(from), to_g = group_of(to)) %>%
    group_by(from = from_g, to = to_g) %>% summarise(flow = sum(flow), .groups = "drop") %>%
    filter(from != to)
}
# Flow-weighted mean affinity of each group
group_aff <- bind_rows(flows %>% transmute(species = from, flow), flows %>% transmute(species = to, flow)) %>%
  left_join(sp_info, by = "species") %>% mutate(g = group_of(species)) %>%
  group_by(g) %>% summarise(MAT_mean = weighted.mean(MAT_mean, flow), MAP_mean = weighted.mean(MAP_mean, flow), .groups = "drop")
gnet <- group_flows(flows) %>%
  left_join(group_flows(flows) %>% select(from = to, to = from, back = flow), by = c("from", "to")) %>%
  mutate(flow = flow - coalesce(back, 0)) %>% filter(flow > 0) %>% select(from, to, flow) %>% arrange(desc(flow))
cat("\nLargest net flows between groups (", grouping, "):\n"); print(head(gnet, 15))

# 5. Chord diagram -------------------------------------------------------------------
short_name <- function(x) ifelse(grepl("^Other|, ", x), x, sub("^(\\w)\\w*\\s", "\\1. ", x))
draw_chord <- function(fl, affinity_col, palette, limits, legend_title, file) {
  aff <- group_aff %>% filter(g %in% c(fl$from, fl$to)) %>% arrange(.data[[affinity_col]])
  col_fun <- circlize::colorRamp2(seq(limits[1], limits[2], length.out = length(palette)), palette)
  grid_col <- setNames(col_fun(pmin(pmax(aff[[affinity_col]], limits[1]), limits[2])), aff$g)
  small <- fl$flow < min_flow_share * sum(fl$flow)
  link_col <- adjustcolor(grid_col[fl$from], alpha.f = 0.65)
  link_col[small] <- adjustcolor(grid_col[fl$from][small], alpha.f = 0.25)
  png(file, width = 9, height = 9, units = "in", res = 300, bg = "white")
  par(mar = c(1, 1, 1, 1), family = "sans")
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2.5, track.margin = c(0.01, 0.01), points.overflow.warning = FALSE)
  chordDiagram(as.data.frame(fl), order = aff$g, grid.col = grid_col, col = link_col,
               directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow",
               diffHeight = -mm_h(2), annotationTrack = "grid", preAllocateTracks = list(track.height = 0.2),
               link.sort = TRUE, link.largest.ontop = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    lab <- short_name(CELL_META$sector.index)
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], lab, facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.85, font = ifelse(grepl("^Other|, ", lab), 1, 3))
  }, bg.border = NA)
  xs <- seq(-1.08, -0.62, length.out = 100)
  rect(xs[-100], 1.02, xs[-1], 1.06, col = col_fun(seq(limits[1], limits[2], length.out = 99)), border = NA, xpd = NA)
  text(c(xs[1], xs[100]), 0.98, limits, cex = 0.8, xpd = NA)
  text(mean(xs), 1.10, legend_title, cex = 0.85, xpd = NA)
  circos.clear(); dev.off(); file
}
gdraw <- gnet %>% filter(flow >= draw_min_share * sum(flow))
f_mat_net <- draw_chord(gdraw, "MAT_mean", COL_MAT, MAT_LIMITS, "Thermal affinity (MAT, °C)",
                        this_module_dir("figures", paste0("understory_replacement_chord_", grouping, "_MAT_net.png")))
f_map_net <- draw_chord(gdraw, "MAP_mean", COL_MAP, MAP_LIMITS, "Moisture affinity (MAP, mm)",
                        this_module_dir("figures", paste0("understory_replacement_chord_", grouping, "_MAP_net.png")))

}

# 6. Open the figures -----------------------------------------------------------------
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(f_mat_net); open_fig(f_map_net)
