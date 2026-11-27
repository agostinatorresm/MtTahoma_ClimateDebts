# ------------------------------------------------------------------------------
# 07_species_contribution_by_lifeform.R
# Module: pipeline/02_analysis/01_thermophilization/02_understory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: aggregate the same abundance-shift decomposition used in
# 05_species_contribution.R up to life-form categories (e.g. forb,
# graminoid, shrub, fern), to check whether the community-level CTI/CMI
# shift is concentrated in a particular growth form rather than spread
# evenly across the understory flora.
#
# Method: identical to 05 at the plot x species level (exact decomposition
# of the plot's change in CTI/CMI), so group contributions of a plot sum
# exactly to that plot's change in index. Each species' contribution is
# summed within its group at each plot, and the per-plot totals are averaged
# across the plots where the group has at least one pool species.
#
# Groups:
#   life form: spp_info.csv (01_cleaning/outputs), column life_form; species
#     without an entry are grouped as "Unknown";
#   woodiness: 04_traits/outputs/understory_traits.csv (species-level
#     records, genus imputation or literature).
#
# Species: shared pool with climatic affinity, as in 01 and 05. Species
# present only in the 50-m2 census have no cover and weight 0.
#
# Output:
#   - outputs/understory_lifeform_contribution.csv: per life form -- number
#     of plots and species involved, mean contribution to CTI and CMI
#   - outputs/understory_woodiness_contribution.csv: same, by woodiness
#   - outputs/understory_lifeform_contribution_by_process.csv,
#     outputs/understory_woodiness_contribution_by_process.csv and
#     outputs/understory_total_contribution_by_process.csv: contributions
#     split into gains and losses of species warmer/colder (wetter/drier) than
#     the community
#   - figures/understory_lifeform_contribution_CTI.png
#   - figures/understory_lifeform_contribution_CMI.png
#   - figures/understory_woodiness_contribution_CTI.png
#   - figures/understory_woodiness_contribution_CMI.png
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

# 1. Load ---------------------------------------------------------------------
understory_long <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE)

species_affinity <- read_csv(affinity_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)

# ASSUMPTION: species-name column "species" -- see header note.
spp_info <- read_csv(cleaning_dir("outputs", "spp_info.csv"), show_col_types = FALSE) %>%
  select(species, life_form)

# WOODINESS -- see header note. Loaded from the
# trait-imputation module's output rather than spp_info.csv, since that is
# where the genus-level TRY+GIFT imputed woodiness values live.
traits_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

understory_traits <- read_csv(traits_dir("outputs", "understory_traits.csv"), show_col_types = FALSE) %>%
  select(species, woodiness)

# Shared pool with climatic affinity (same as 01_compute_cti_cmi_understory.R)
understory_pool <- read_csv(
  species_pool_dir("outputs", "understory_species_pool_broad_full_clim_cov.csv"),
  show_col_types = FALSE
) %>%
  pull(species)


# 2. Per plot x species x period cover, filled to 0 where a species is absent
#    in one of the two periods -- identical to 05, section 2 -------------------
understory_weight <- understory_long %>%
  filter(occurrence == 1) %>%
  select(plot_id, period, species, weight = cover) %>%
  # species present only in the 50-m2 census have no cover: weight 0
  mutate(weight = coalesce(weight, 0)) %>%
  complete(nesting(plot_id, species), period, fill = list(weight = 0))

# 3. Restrict to species with a climatic affinity AND in the broad species ----
#    pool (see header note), attach life form, compute each species' per-plot
#    contribution
contribution_by_plot_all <- understory_weight %>%
  filter(species %in% species_affinity$species, species %in% understory_pool) %>%
  pivot_wider(names_from = period, values_from = weight, values_fill = 0) %>%
  left_join(species_affinity, by = "species") %>%
  left_join(spp_info, by = "species") %>%
  left_join(understory_traits, by = "species") %>%
  mutate(
    life_form = coalesce(life_form, "Unknown"),
    woodiness = coalesce(woodiness, "Unknown")
  )

# Plots where the index is defined in both periods (see 05's header)
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
    weight_change   = present - past,
    rel_weight_change   = p_present - p_past,
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

cat("Species entering the life-form/woodiness aggregation (affinity + broad pool):",
    n_distinct(contribution_by_plot$species), "\n")

n_no_lifeform <- contribution_by_plot %>% filter(life_form == "Unknown") %>% distinct(species) %>% nrow()
cat("Species with no life_form entry in spp_info.csv (grouped as 'Unknown'):", n_no_lifeform, "\n")

n_no_woodiness <- contribution_by_plot %>% filter(woodiness == "Unknown") %>% distinct(species) %>% nrow()
cat("Species with no woodiness entry in understory_traits.csv (grouped as 'Unknown'):", n_no_woodiness, "\n")
cat("-> coverage comparison: life_form 'Unknown' is the pre-existing, already-flagged gap (see header);",
    "woodiness is shown here as the better-covered alternative categorical breakdown.\n")

# 4. Sum contribution within life form, per plot -------------------------------
contribution_by_plot_lifeform <- contribution_by_plot %>%
  group_by(plot_id, life_form) %>%
  summarise(
    n_species_lifeform        = n_distinct(species),
    total_contribution_to_CTI = sum(contribution_to_CTI, na.rm = TRUE),
    total_contribution_to_CMI = sum(contribution_to_CMI, na.rm = TRUE),
    .groups = "drop"
  )

# 5. Average across plots where the life form is present (same convention as
#    05's species-level average -- across plots where it occurs, not zero-
#    filled for plots where it is absent) --------------------------------------
# mean_contribution_*: mean over the plots where the group occurs;
# community_contribution_*: sum over plots divided by all plots with a
# defined index (contribution to the mean change in index; groups sum to the
# mean plot-level change)
n_plots_index <- n_distinct(contribution_by_plot_lifeform$plot_id)
lifeform_summary <- contribution_by_plot_lifeform %>%
  group_by(life_form) %>%
  summarise(
    n_plots               = n(),
    n_species_total        = sum(n_species_lifeform),
    mean_contribution_CTI  = mean(total_contribution_to_CTI, na.rm = TRUE),
    mean_contribution_CMI  = mean(total_contribution_to_CMI, na.rm = TRUE),
    community_contribution_CTI = sum(total_contribution_to_CTI, na.rm = TRUE) / n_plots_index,
    community_contribution_CMI = sum(total_contribution_to_CMI, na.rm = TRUE) / n_plots_index,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_contribution_CTI))

write_csv(lifeform_summary, this_module_dir("outputs", "understory_lifeform_contribution.csv"))

cat("\nLife-form contribution to CTI and CMI change:\n")
print(lifeform_summary)

# 4b/5b. Same aggregation repeated with woodiness in place of life_form -------
contribution_by_plot_woodiness <- contribution_by_plot %>%
  group_by(plot_id, woodiness) %>%
  summarise(
    n_species_woodiness       = n_distinct(species),
    total_contribution_to_CTI = sum(contribution_to_CTI, na.rm = TRUE),
    total_contribution_to_CMI = sum(contribution_to_CMI, na.rm = TRUE),
    .groups = "drop"
  )

woodiness_summary <- contribution_by_plot_woodiness %>%
  group_by(woodiness) %>%
  summarise(
    n_plots               = n(),
    n_species_total        = sum(n_species_woodiness),
    mean_contribution_CTI  = mean(total_contribution_to_CTI, na.rm = TRUE),
    mean_contribution_CMI  = mean(total_contribution_to_CMI, na.rm = TRUE),
    community_contribution_CTI = sum(total_contribution_to_CTI, na.rm = TRUE) / n_plots_index,
    community_contribution_CMI = sum(total_contribution_to_CMI, na.rm = TRUE) / n_plots_index,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_contribution_CTI))

write_csv(woodiness_summary, this_module_dir("outputs", "understory_woodiness_contribution.csv"))

cat("\nWoodiness contribution to CTI and CMI change:\n")
print(woodiness_summary)

# 5c. Contributions split into four processes ---------------------------------
# The sign of each plot x species contribution is set by two signs: whether
# the species' relative cover increased or decreased in that plot, and whether
# its affinity is above or below the plot's past index. The four combinations
# are the processes that shift a community index:
#   CTI: gain of warm-affiliated (+), loss of cold-affiliated (+),
#        gain of cold-affiliated (-), loss of warm-affiliated (-)
#   CMI: gain of wet-affiliated (+),  loss of dry-affiliated (+),
#        gain of dry-affiliated (-),  loss of wet-affiliated (-)
# "Warm-affiliated" means warmer than the community of that plot (MAT_mean >
# CTI_past), not warm in absolute terms. For each group the four parts add up
# to community_contribution_CTI/CMI above. Relative cover increases are
# shifted upward in the present survey, because fewer species had cover per
# plot (Table S3.1), so the gain/loss split is read as a relative indication.
summarise_by_process <- function(df, group_col) {
  classify <- function(d, affinity_col, index_past_col, contribution_col, index, hi, lo) {
    d %>%
      filter(rel_weight_change != 0) %>%
      mutate(
        index    = index,
        change   = if_else(rel_weight_change > 0, "Gain", "Loss"),
        affinity = if_else(.data[[affinity_col]] > .data[[index_past_col]], hi, lo),
        process  = paste0(change, " of ", affinity, "-affiliated"),
        contribution = .data[[contribution_col]]
      )
  }
  bind_rows(
    classify(df, "MAT_mean", "CTI_past", "contribution_to_CTI", "CTI", "warm", "cold"),
    classify(df, "MAP_mean", "CMI_past", "contribution_to_CMI", "CMI", "wet",  "dry")
  ) %>%
    group_by(.data[[group_col]], index, process) %>%
    summarise(
      n_species              = n_distinct(species),
      community_contribution = sum(contribution, na.rm = TRUE) / n_plots_index,
      .groups = "drop"
    )
}

lifeform_process  <- summarise_by_process(contribution_by_plot, "life_form")
woodiness_process <- summarise_by_process(contribution_by_plot, "woodiness")
# Whole understory, for reference
total_process <- summarise_by_process(mutate(contribution_by_plot, group = "All understory"), "group")

write_csv(lifeform_process,  this_module_dir("outputs", "understory_lifeform_contribution_by_process.csv"))
write_csv(woodiness_process, this_module_dir("outputs", "understory_woodiness_contribution_by_process.csv"))
write_csv(total_process,     this_module_dir("outputs", "understory_total_contribution_by_process.csv"))

cat("\nContribution by process, whole understory:\n")
print(total_process, n = Inf)
cat("\nContribution by process and life form:\n")
print(lifeform_process %>% arrange(index, life_form, process), n = Inf)

# 6. Figures: life forms ranked by mean contribution ---------------------------
plot_lifeform <- function(df, group_col, contribution_col, x_lab, title_lab) {
  df %>%
    mutate(.group = fct_reorder(.data[[group_col]], .data[[contribution_col]])) %>%
    ggplot(aes(x = .data[[contribution_col]], y = .group, fill = .data[[contribution_col]] > 0)) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#2E8B57", "FALSE" = "#B03A2E"), guide = "none") +
    labs(x = x_lab, y = NULL, title = title_lab) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 14, face = "plain"))
}

p_lifeform_CTI <- plot_lifeform(
  lifeform_summary, "life_form", "community_contribution_CTI",
  "Contribution to mean CTI change (°C)",
  "Understory life forms ranked by contribution to community CTI change"
)
p_lifeform_CTI

p_lifeform_CMI <- plot_lifeform(
  lifeform_summary, "life_form", "community_contribution_CMI",
  "Contribution to mean CMI change (mm)",
  "Understory life forms ranked by contribution to community CMI change"
)
p_lifeform_CMI

ggsave(this_module_dir("figures", "understory_lifeform_contribution_CTI.png"), p_lifeform_CTI,
       width = 8, height = max(3, nrow(lifeform_summary) * 0.5), dpi = 300, bg = "white")

ggsave(this_module_dir("figures", "understory_lifeform_contribution_CMI.png"), p_lifeform_CMI,
       width = 8, height = max(3, nrow(lifeform_summary) * 0.5), dpi = 300, bg = "white")

# 6b. Figures: same, by woodiness ----------------------------------------------
p_woodiness_CTI <- plot_lifeform(
  woodiness_summary, "woodiness", "community_contribution_CTI",
  "Contribution to mean CTI change (°C)",
  "Understory woodiness classes ranked by contribution to community CTI change"
)
p_woodiness_CTI

p_woodiness_CMI <- plot_lifeform(
  woodiness_summary, "woodiness", "community_contribution_CMI",
  "Contribution to mean CMI change (mm)",
  "Understory woodiness classes ranked by contribution to community CMI change"
)
p_woodiness_CMI

ggsave(this_module_dir("figures", "understory_woodiness_contribution_CTI.png"), p_woodiness_CTI,
       width = 8, height = max(3, nrow(woodiness_summary) * 0.5), dpi = 300, bg = "white")

ggsave(this_module_dir("figures", "understory_woodiness_contribution_CMI.png"), p_woodiness_CMI,
       width = 8, height = max(3, nrow(woodiness_summary) * 0.5), dpi = 300, bg = "white")

# 7. Summary --------------------------------------------------------------------
cat("\n========================================\n")
cat("SPECIES CONTRIBUTION BY LIFE FORM / WOODINESS -- COMPLETE\n")
cat("========================================\n")
cat("Saved:\n")
cat(" - outputs/understory_lifeform_contribution.csv\n")
cat(" - outputs/understory_woodiness_contribution.csv\n")
cat(" - figures/understory_lifeform_contribution_CTI.png\n")
cat(" - figures/understory_lifeform_contribution_CMI.png\n")
cat(" - figures/understory_woodiness_contribution_CTI.png\n")
cat(" - figures/understory_woodiness_contribution_CMI.png\n")
