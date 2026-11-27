# ------------------------------------------------------------------------------
# 05_species_contribution.R
# Module: pipeline/02_analysis/01_thermophilization/01_overstory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: rank overstory species by how much each one contributes to the
# community-level CTI and CMI shifts between periods -- i.e. which species
# are actually driving (or opposing) the thermophilization/hygrophilization
# signal reported in 01_compute_cti_cmi.R, rather than looking only at the
# aggregate rates.
#
# METHOD: exact decomposition of the plot-level change in
# the community-weighted mean, applied in parallel to both variables. CTI/CMI
# are weighted means, so they depend on RELATIVE dominance, p_i =
# basal_area_spp_perc_i / sum of basal_area_spp_perc over the pool species in
# that plot and period (same pool and weights as 01_compute_cti_cmi.R). With
# t_i the species' affinity (MAT_mean for CTI, MAP_mean for CMI; single
# truncated-mean values, not bootstrap resamples),
#     index = sum_i p_i * t_i
#     change in index = sum_i (p_i,present - p_i,past) * (t_i - index_past)
# where index_past is the plot's own past CTI (or CMI). Subtracting the
# constant index_past is allowed because the changes in relative dominance
# sum to zero within a plot, and it makes each term interpretable: a species
# contributes positively when it gains relative dominance and is warmer
# (wetter) than the plot's initial community, or loses relative dominance
# and is colder (drier) than it. The terms sum exactly to the plot's change
# in index, which is verified in the console (section 3).
# The decomposition uses point affinities, whereas the rates in 01 use the
# median across affinity resamples, so the summed contributions are close to,
# but not identical with, the reported change in the median index. Plots with
# no pool species in one of the two periods have an undefined index and are
# excluded from the decomposition.
#
# Per-species contributions are then averaged across all plots where that
# species occurs in either period, giving one number per species per
# variable: positive means the species is, on average, pushing plot-level
# CTI (or CMI) up (either because it is a warm/wet-affinity species gaining
# basal area, or a cold/dry-affinity species losing it); negative means the
# opposite.
#
# SPECIES POOL: restricted to the same BROAD species pool
# with full climatic coverage used in 01_compute_cti_cmi.R
# (overstory_species_pool_broad_full_clim_cov.csv), in place of the plain
# "has a climatic affinity" filter -- so the species ranked here are exactly
# the species entering the CTI/CMI calculation itself, not a superset of it.
#
# Output:
#   - outputs/overstory_species_contribution_CTI.csv
#   - outputs/overstory_species_contribution_CMI.csv
#   - figures/overstory_species_contribution_CTI.png
#   - figures/overstory_species_contribution_CMI.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)

# 0a. Module directories --------------------------------------------------------
this_module_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)
cleaning_dir     <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
affinity_dir     <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
species_pool_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)

dir.create(this_module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

# 1. Load -------------------------------------------------------------------------
overstory_long <- read_csv(cleaning_dir("outputs", "overstory_long.csv"), show_col_types = FALSE)

species_affinity <- read_csv(affinity_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)

# BROAD species pool, full climatic coverage -- same pool used in
# 01_compute_cti_cmi.R, see header note above.
overstory_pool <- read_csv(
  species_pool_dir("outputs", "overstory_species_pool_broad_full_clim_cov.csv"),
  show_col_types = FALSE
) %>%
  pull(species)

# 2. Per plot x species x period dominance, filled to 0 where a species is
#    absent in one of the two periods (needed so a species gained/lost
#    entirely between periods still gets a non-zero weight_change) -----------
overstory_weight <- overstory_long %>%
  filter(occurrence == 1) %>%
  select(plot_id, period, species, weight = basal_area_spp_perc) %>%
  complete(nesting(plot_id, species), period, fill = list(weight = 0))

cat("Plot x species combinations (both periods filled):", nrow(overstory_weight), "\n")

# 3. Restrict to species with a climatic affinity AND in the broad species -----
#    pool, compute each species' contribution to CTI and CMI change at each
#    plot where it occurs in either period
contribution_by_plot_all <- overstory_weight %>%
  filter(species %in% species_affinity$species, species %in% overstory_pool) %>%
  pivot_wider(names_from = period, values_from = weight, values_fill = 0) %>%
  left_join(species_affinity, by = "species")

# Plots where the index is defined in both periods (some pool species with
# weight > 0 in each period)
plots_defined <- contribution_by_plot_all %>%
  group_by(plot_id) %>%
  summarise(total_past = sum(past), total_present = sum(present), .groups = "drop") %>%
  filter(total_past > 0, total_present > 0)

cat("Plots with the index defined in both periods:", nrow(plots_defined),
    "of", n_distinct(contribution_by_plot_all$plot_id), "\n")

contribution_by_plot <- contribution_by_plot_all %>%
  filter(plot_id %in% plots_defined$plot_id) %>%
  group_by(plot_id) %>%
  mutate(
    p_past          = past / sum(past),
    p_present       = present / sum(present),
    CTI_past        = sum(p_past * MAT_mean),
    CMI_past        = sum(p_past * MAP_mean),
    CTI_present     = sum(p_present * MAT_mean),
    CMI_present     = sum(p_present * MAP_mean),
    weight_change   = present - past,                    # percentage points of basal area (absolute)
    rel_weight_change   = p_present - p_past,            # change in relative dominance
    contribution_to_CTI = rel_weight_change * (MAT_mean - CTI_past),
    contribution_to_CMI = rel_weight_change * (MAP_mean - CMI_past)
  ) %>%
  ungroup()

# Check: contributions must sum, per plot, to the plot's change in index
decomposition_check <- contribution_by_plot %>%
  group_by(plot_id) %>%
  summarise(
    err_CTI = abs(sum(contribution_to_CTI) - (first(CTI_present) - first(CTI_past))),
    err_CMI = abs(sum(contribution_to_CMI) - (first(CMI_present) - first(CMI_past))),
    .groups = "drop"
  )
cat("Decomposition check -- max abs discrepancy per plot: CTI",
    signif(max(decomposition_check$err_CTI), 3), ", CMI",
    signif(max(decomposition_check$err_CMI), 3), "(should be ~0)\n")

cat("Species entering the contribution ranking (affinity + broad pool):",
    n_distinct(contribution_by_plot$species), "\n")

# 4. Aggregate to one value per species, per variable -----------------------------
# Two summaries per species:
#   mean_contribution: mean over the plots where the species occurs (size of
#     its effect where present);
#   community_contribution: sum over plots divided by all plots with a
#     defined index (contribution to the mean change in index; the species'
#     values sum to the mean plot-level change, so frequent species weigh
#     more). Species are ranked and plotted by community_contribution.
n_plots_index <- n_distinct(contribution_by_plot$plot_id)
summarise_contribution <- function(df, contribution_col, affinity_col) {
  df %>%
    group_by(species) %>%
    summarise(
      affinity            = first(.data[[affinity_col]]),
      n_plots             = n(),
      mean_contribution   = mean(.data[[contribution_col]], na.rm = TRUE),
      total_contribution  = sum(.data[[contribution_col]], na.rm = TRUE),
      community_contribution = total_contribution / n_plots_index,
      mean_weight_change  = mean(weight_change, na.rm = TRUE),
      mean_rel_weight_change = mean(rel_weight_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(community_contribution))
}

species_contribution_CTI <- summarise_contribution(contribution_by_plot, "contribution_to_CTI", "MAT_mean") %>%
  rename(MAT_mean = affinity)

species_contribution_CMI <- summarise_contribution(contribution_by_plot, "contribution_to_CMI", "MAP_mean") %>%
  rename(MAP_mean = affinity)

write_csv(species_contribution_CTI, this_module_dir("outputs", "overstory_species_contribution_CTI.csv"))
write_csv(species_contribution_CMI, this_module_dir("outputs", "overstory_species_contribution_CMI.csv"))

cat("\n-- CTI --\n")
cat("Top 5 species pushing CTI up (warm-affinity species gaining dominance,",
    "or cold-affinity species losing it):\n")
print(species_contribution_CTI %>% slice_head(n = 5) %>% select(species, MAT_mean, n_plots, community_contribution, mean_contribution))
cat("\nTop 5 species pushing CTI down:\n")
print(species_contribution_CTI %>% slice_tail(n = 5) %>% select(species, MAT_mean, n_plots, community_contribution, mean_contribution))

cat("\n-- CMI --\n")
cat("Top 5 species pushing CMI up (wet-affinity species gaining dominance,",
    "or dry-affinity species losing it):\n")
print(species_contribution_CMI %>% slice_head(n = 5) %>% select(species, MAP_mean, n_plots, community_contribution, mean_contribution))
cat("\nTop 5 species pushing CMI down:\n")
print(species_contribution_CMI %>% slice_tail(n = 5) %>% select(species, MAP_mean, n_plots, community_contribution, mean_contribution))

# 5. Figures: all species ranked by mean contribution, colored by their own affinity
plot_contribution <- function(df, affinity_col, color_scale, x_lab, title_lab) {
  df %>%
    mutate(species = fct_reorder(species, community_contribution)) %>%
    ggplot(aes(x = community_contribution, y = species, color = .data[[affinity_col]])) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = community_contribution, y = species, yend = species), linewidth = 0.6) +
    geom_point(size = 2.5) +
    color_scale +
    labs(x = x_lab, y = NULL, title = title_lab) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "plain"))
}

p_contribution_CTI <- plot_contribution(
  species_contribution_CTI, "MAT_mean",
  scale_color_viridis_c(option = "plasma", direction = -1, name = "MAT affinity (°C)"),
  "Contribution to mean CTI change (°C)",
  "Overstory species ranked by contribution to community CTI change"
)
p_contribution_CTI

p_contribution_CMI <- plot_contribution(
  species_contribution_CMI, "MAP_mean",
  scale_color_viridis_c(option = "viridis", direction = -1, name = "MAP affinity (mm)"),
  "Contribution to mean CMI change (mm)",
  "Overstory species ranked by contribution to community CMI change"
)
p_contribution_CMI

ggsave(this_module_dir("figures", "overstory_species_contribution_CTI.png"), p_contribution_CTI,
       width = 8, height = max(4, nrow(species_contribution_CTI) * 0.28), dpi = 300, bg = "white")

ggsave(this_module_dir("figures", "overstory_species_contribution_CMI.png"), p_contribution_CMI,
       width = 8, height = max(4, nrow(species_contribution_CMI) * 0.28), dpi = 300, bg = "white")

# 6. Summary ----------------------------------------------------------------------
cat("\n========================================\n")
cat("SPECIES CONTRIBUTION TO CTI AND CMI CHANGE -- COMPLETE\n")
cat("========================================\n")
cat("Saved:\n")
cat(" - outputs/overstory_species_contribution_CTI.csv\n")
cat(" - outputs/overstory_species_contribution_CMI.csv\n")
cat(" - figures/overstory_species_contribution_CTI.png\n")
cat(" - figures/overstory_species_contribution_CMI.png\n")
