# ------------------------------------------------------------------------------
# 05_species_contribution.R
# Module: pipeline/02_analysis/01_thermophilization/02_understory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: rank understory species by how much each one contributes to the
# community-level CTI and CMI shifts between periods -- i.e. which species
# are actually driving (or opposing) the thermophilization/hygrophilization
# signal reported in 01_compute_cti_cmi_understory.R, rather than looking only at the
# aggregate rates. Understory counterpart of
# 01_overstory/code/05_species_contribution.R.
#
# Method: exact decomposition of the plot-level change in the community-
# weighted mean, for CTI and CMI. The indices depend on relative cover,
# p_i = cover_i / total cover of the pool species in the plot and survey.
# With t_i the species' affinity (MAT_mean for CTI, MAP_mean for CMI; mean of
# the truncated distribution),
#     index = sum_i p_i * t_i
#     change in index = sum_i (p_i,present - p_i,past) * (t_i - index_past)
# where index_past is the plot's past index. Subtracting index_past is
# possible because the changes in relative cover sum to zero within a plot;
# a species contributes positively when it gains relative cover and is warmer
# (wetter) than the initial community, or loses relative cover and is colder
# (drier) than it. The terms sum exactly to the plot's change in index
# (checked in section 3). The decomposition uses point affinities, whereas the
# rates in 01 use the median across affinity resamples, so the summed
# contributions are close to, but not identical with, the reported change in
# the median index. Plots with no pool species with cover in one survey have
# an undefined index and are excluded. Species present only in the 50-m2
# census have no cover and weight 0.
#
# Species: shared pool with climatic affinity
# (understory_species_pool_broad_full_clim_cov.csv), as in 01.
#
# Per-species contributions are then averaged across all plots where that
# species occurs in either period, giving one number per species per
# variable: positive means the species is, on average, pushing plot-level
# CTI (or CMI) up (either because it is a warm/wet-affinity species gaining
# cover, or a cold/dry-affinity species losing it); negative means the
# opposite.
#
# Output:
#   - outputs/understory_species_contribution_CTI.csv
#   - outputs/understory_species_contribution_CMI.csv
#   - figures/understory_species_contribution_CTI.png
#   - figures/understory_species_contribution_CMI.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)

# 0a. Module directories --------------------------------------------------------
this_module_dir  <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "02_understory", ...)
cleaning_dir     <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
affinity_dir     <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
species_pool_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)

dir.create(this_module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

# 1. Load -------------------------------------------------------------------------
understory_long <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE)

species_affinity <- read_csv(affinity_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)

# Shared pool with climatic affinity (same as 01_compute_cti_cmi_understory.R)
understory_pool <- read_csv(
  species_pool_dir("outputs", "understory_species_pool_broad_full_clim_cov.csv"),
  show_col_types = FALSE
) %>%
  pull(species)


# 2. Per plot x species x period cover, filled to 0 where a species is
#    absent in one of the two periods (needed so a species gained/lost
#    entirely between periods still gets a non-zero weight_change) -----------
understory_weight <- understory_long %>%
  filter(occurrence == 1) %>%
  select(plot_id, period, species, weight = cover) %>%
  # species present only in the 50-m2 census have no cover: weight 0
  mutate(weight = coalesce(weight, 0)) %>%
  complete(nesting(plot_id, species), period, fill = list(weight = 0))

cat("Plot x species combinations (both periods filled):", nrow(understory_weight), "\n")

# 3. Restrict to species with a climatic affinity AND in the broad species -----
#    pool, compute each species' contribution to CTI and CMI change at each
#    plot where it occurs in either period
contribution_by_plot_all <- understory_weight %>%
  filter(species %in% species_affinity$species, species %in% understory_pool) %>%
  pivot_wider(names_from = period, values_from = weight, values_fill = 0) %>%
  left_join(species_affinity, by = "species")

# Plots where the index is defined in both periods (some pool species with
# cover > 0 in each period)
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
    weight_change   = present - past,                    # percentage points of cover (absolute)
    rel_weight_change   = p_present - p_past,            # change in relative cover
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

write_csv(species_contribution_CTI, this_module_dir("outputs", "understory_species_contribution_CTI.csv"))
write_csv(species_contribution_CMI, this_module_dir("outputs", "understory_species_contribution_CMI.csv"))

cat("\n-- CTI --\n")
cat("Top 5 species pushing CTI up (warm-affinity species gaining cover,",
    "or cold-affinity species losing it):\n")
print(species_contribution_CTI %>% slice_head(n = 5) %>% select(species, MAT_mean, n_plots, community_contribution, mean_contribution))
cat("\nTop 5 species pushing CTI down:\n")
print(species_contribution_CTI %>% slice_tail(n = 5) %>% select(species, MAT_mean, n_plots, community_contribution, mean_contribution))

cat("\n-- CMI --\n")
cat("Top 5 species pushing CMI up (wet-affinity species gaining cover,",
    "or dry-affinity species losing it):\n")
print(species_contribution_CMI %>% slice_head(n = 5) %>% select(species, MAP_mean, n_plots, community_contribution, mean_contribution))
cat("\nTop 5 species pushing CMI down:\n")
print(species_contribution_CMI %>% slice_tail(n = 5) %>% select(species, MAP_mean, n_plots, community_contribution, mean_contribution))

# 5. Figures: top/bottom N species by contribution, colored by their own affinity.
#    The full ranking is in the CSV (section 4); with ~117 species carrying a
#    climatic affinity (the broad pool), a one-row-per-species figure exceeds
#    a legible (and a ggsave-plottable) height, so the figure itself is
#    restricted to the two tails of the ranking, where the pattern actually is.
N_SPECIES_SHOWN <- 20

plot_contribution <- function(df, affinity_col, color_scale, x_lab, title_lab, n_show = N_SPECIES_SHOWN) {
  df_show <- bind_rows(
    slice_max(df, order_by = community_contribution, n = n_show),
    slice_min(df, order_by = community_contribution, n = n_show)
  ) %>%
    distinct(species, .keep_all = TRUE)
  
  df_show %>%
    mutate(species = fct_reorder(species, community_contribution)) %>%
    ggplot(aes(x = community_contribution, y = species, color = .data[[affinity_col]])) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = community_contribution, y = species, yend = species), linewidth = 0.6) +
    geom_point(size = 2.5) +
    color_scale +
    labs(x = x_lab, y = NULL, title = paste0(title_lab, " (top/bottom ", n_show, ")")) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "plain"))
}

p_contribution_CTI <- plot_contribution(
  species_contribution_CTI, "MAT_mean",
  scale_color_viridis_c(option = "plasma", direction = -1, name = "MAT affinity (°C)"),
  "Contribution to mean CTI change (°C)",
  "Understory species ranked by contribution to community CTI change"
)
p_contribution_CTI

p_contribution_CMI <- plot_contribution(
  species_contribution_CMI, "MAP_mean",
  scale_color_viridis_c(option = "viridis", direction = -1, name = "MAP affinity (mm)"),
  "Contribution to mean CMI change (mm)",
  "Understory species ranked by contribution to community CMI change"
)
p_contribution_CMI

ggsave(this_module_dir("figures", "understory_species_contribution_CTI.png"), p_contribution_CTI,
       width = 8, height = max(4, 2 * N_SPECIES_SHOWN * 0.28), dpi = 300, bg = "white")

ggsave(this_module_dir("figures", "understory_species_contribution_CMI.png"), p_contribution_CMI,
       width = 8, height = max(4, 2 * N_SPECIES_SHOWN * 0.28), dpi = 300, bg = "white")

# 6. Summary ----------------------------------------------------------------------
cat("\n========================================\n")
cat("SPECIES CONTRIBUTION TO CTI AND CMI CHANGE -- COMPLETE\n")
cat("========================================\n")
cat("Saved:\n")
cat(" - outputs/understory_species_contribution_CTI.csv\n")
cat(" - outputs/understory_species_contribution_CMI.csv\n")
cat(" - figures/understory_species_contribution_CTI.png\n")
cat(" - figures/understory_species_contribution_CMI.png\n")
