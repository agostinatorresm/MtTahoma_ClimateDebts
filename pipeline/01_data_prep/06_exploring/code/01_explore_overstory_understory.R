# -------------------------------------------------------------------------
# Author: Agostina Torres
# Project: MORAclimshift
# Script: 01_explore_overstory_understory.R
# Description:
# Exploratory summaries for BOTH overstory (Tree) and understory (Veg),
# comparing past (1975) and present. Every table is saved both as a .csv
# (for further analysis) and as a .png (for quick viewing).
#
# Sections:
#   A - Species coverage: how many species are shared between periods vs
#       found in only one (lollipop, both strata)
#   B - Plot coverage: same idea, at the plot level (lollipop, both strata)
#   C - Final scope summary: plots/species kept for analysis,
#       after requiring real presence in both periods
#   D - General cover/abundance by period, with past-present differences
#   E - Plot coverage by resurvey year (present) and by elevation (past,
#       used as a fixed per-plot attribute, since present-day elevation
#       was not recorded separately)
#   F - Turnover: colonization/extinction/persistence between periods, as
#       a bar chart + heatmap (both strata). An alluvial version
#       (create_alluvial_plot) is defined below but not used, because of
#       an incompatibility between ggalluvial and recent ggplot2 versions
#       ("could not find function gg_par").
#       calculate_turnover/print_turnover_summary/create_turnover_heatmap
#       follow the turnover functions of the earlier MORA and NOCA
#       cleaning scripts.
#   G - Rank change: how each species' abundance rank shifts from past to
#       present (bump/slope chart), both layers
#   H - Elevation shift: mean elevation of occurrence per species, past vs
#       present - are species occurring higher up now? (both layers,
#       elevation from the past survey used as a fixed per-plot proxy)
#   I - Elevation DISTRIBUTION per species (not just the mean): min/max/
#       mean/median/sd/skewness/kurtosis, past vs present, plus a small-
#       multiples density plot per species (both layers; understory
#       limited to top 24 species by occurrence for readability)
#   J - Understory grouped by LIFE FORM (Forb/Shrub/Grass/Fern/Sedge/
#       Rush) instead of by species - covers all 221 species pooled into
#       their life-form category
#
# Outputs (pipeline/01_data_prep/02_exploring/outputs/):
#   *.csv for every table, *.png for every table and figure
# -------------------------------------------------------------------------

library(tidyverse)
library(here)
library(gridExtra)
# library(hrbrthemes)  # no longer needed - switched to theme_minimal() everywhere
# library(ggalluvial)  # not used - see Section F note
library(viridis)

cleaning_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
explore_dir  <- here("pipeline", "01_data_prep", "02_exploring", "outputs")
dir.create(explore_dir, recursive = TRUE, showWarnings = FALSE)

col_past    <- "#7B5EA7"
col_present <- "#E8821A"
col_both    <- "#4A9B8E"

save_table_png <- function(df, filename, title = NULL) {
  g <- tableGrob(df, rows = NULL)
  if (!is.null(title)) {
    title_grob <- grid::textGrob(title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    g <- gtable::gtable_add_rows(g, heights = grid::grobHeight(title_grob) + unit(6, "mm"), pos = 0)
    g <- gtable::gtable_add_grob(g, title_grob, 1, 1, 1, ncol(g), name = "title")
  }
  ggsave(file.path(explore_dir, filename), g, width = min(14, 2 + ncol(df) * 1.4),
         height = 0.4 + 0.3 * nrow(df), dpi = 200, limitsize = FALSE)
}

# Saves every actual chart/plot in BOTH .png and .pdf - filename_base has
# no extension (e.g. "A_species_coverage_lollipop"). Uses theme_minimal()
# throughout (not hrbrthemes::theme_ipsum()) specifically so the standard
# PDF device works everywhere - theme_ipsum's custom font broke the
# default pdf() device, and cairo_pdf() isn't available on this machine's
# R build either ("unable to start device 'cairo_pdf'").
save_plot <- function(plot, filename_base, width, height) {
  ggsave(file.path(explore_dir, paste0(filename_base, ".png")), plot, width = width, height = height, dpi = 200)
  ggsave(file.path(explore_dir, paste0(filename_base, ".pdf")), plot, width = width, height = height)
}


# ==============================================================================
# LOAD DATA
# ==============================================================================

overstory_long     <- read_csv(file.path(cleaning_dir, "overstory_long.csv"), show_col_types = FALSE)
understory_long    <- read_csv(file.path(cleaning_dir, "understory_long.csv"), show_col_types = FALSE)
spp_info           <- read_csv(file.path(cleaning_dir, "spp_info.csv"), show_col_types = FALSE)
overstory_plot_cov <- read_csv(file.path(cleaning_dir, "overstory_plot_coverage.csv"), show_col_types = FALSE)
understory_plot_cov<- read_csv(file.path(cleaning_dir, "understory_plot_coverage.csv"), show_col_types = FALSE)
plot_info          <- read_csv(file.path(cleaning_dir, "plot_info.csv"), show_col_types = FALSE)
plot_env_by_year   <- read_csv(file.path(cleaning_dir, "plot_environment_by_year.csv"), show_col_types = FALSE)


# ==============================================================================
# SECTION A — Species coverage: shared vs unique-to-one-period, both layers
# ==============================================================================

species_coverage <- spp_info |>
  filter(canopy_type %in% c("over", "under")) |>
  mutate(
    stratum = if_else(canopy_type == "over", "Overstory", "Understory"),
    coverage = case_when(
      presence_coverage == "both"         ~ "Both periods",
      presence_coverage == "only_present" ~ "Present only",
      presence_coverage == "only_past"    ~ "Past only",
      TRUE                                 ~ "No real presence"
    )
  ) |>
  count(stratum, coverage)

write_csv(species_coverage, file.path(explore_dir, "A_species_coverage.csv"))
save_table_png(species_coverage, "A_species_coverage.png", "Species coverage by period")

coverage_order <- c("No real presence", "Past only", "Present only", "Both periods")

pA <- species_coverage |>
  mutate(coverage = factor(coverage, levels = coverage_order)) |>
  ggplot(aes(x = n, y = coverage, color = stratum)) +
  geom_segment(aes(x = 0, xend = n, y = coverage, yend = coverage), linewidth = 1) +
  geom_point(size = 4) +
  geom_text(aes(label = n), hjust = -0.6, size = 3.5, color = "black") +
  facet_wrap(~ stratum, scales = "free_x") +
  scale_color_manual(values = c("Overstory" = col_past, "Understory" = col_present)) +
  labs(title = "Species coverage across periods", x = "Number of species", y = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

save_plot(pA, "A_species_coverage_lollipop", width = 10, height = 5)


# ==============================================================================
# SECTION B — Plot coverage: shared vs unique-to-one-period, both layers
# ==============================================================================

plot_coverage <- bind_rows(
  overstory_plot_cov  |> mutate(stratum = "Overstory"),
  understory_plot_cov |> mutate(stratum = "Understory")
) |>
  mutate(coverage = case_when(
    coverage == "both"          ~ "Both periods",
    coverage == "only_present"  ~ "Present only",
    coverage == "only_past"     ~ "Past only"
  )) |>
  count(stratum, coverage)

write_csv(plot_coverage, file.path(explore_dir, "B_plot_coverage.csv"))
save_table_png(plot_coverage, "B_plot_coverage.png", "Plot coverage by period")

pB <- plot_coverage |>
  mutate(coverage = factor(coverage, levels = coverage_order)) |>
  ggplot(aes(x = n, y = coverage, color = stratum)) +
  geom_segment(aes(x = 0, xend = n, y = coverage, yend = coverage), linewidth = 1) +
  geom_point(size = 4) +
  geom_text(aes(label = n), hjust = -0.6, size = 3.5, color = "black") +
  facet_wrap(~ stratum, scales = "free_x") +
  scale_color_manual(values = c("Overstory" = col_past, "Understory" = col_present)) +
  labs(title = "Plot coverage across periods", x = "Number of plots", y = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

save_plot(pB, "B_plot_coverage_lollipop", width = 10, height = 4)


# ==============================================================================
# SECTION C — Final scope: plots/species kept for analysis
# ==============================================================================

final_scope <- tibble(
  stratum   = c("Overstory", "Understory"),
  n_plots   = c(n_distinct(overstory_long$plot_id), n_distinct(understory_long$plot_id)),
  n_species = c(n_distinct(overstory_long$species), n_distinct(understory_long$species))
)

write_csv(final_scope, file.path(explore_dir, "C_final_scope.csv"))
save_table_png(final_scope, "C_final_scope.png", "Final plots/species kept (real presence in both periods)")


# ==============================================================================
# SECTION D — General cover/abundance by period, with past-present
# differences
# ==============================================================================

overstory_summary <- overstory_long |>
  group_by(period) |>
  summarise(
    total_basal_area_m2 = sum(basal_area_spp_m2, na.rm = TRUE),
    mean_count_per_plot  = mean(count, na.rm = TRUE),
    mean_ba_perc_per_plot = mean(basal_area_spp_perc, na.rm = TRUE),
    .groups = "drop"
  )

understory_summary <- understory_long |>
  group_by(period) |>
  summarise(
    mean_cover_per_plot = mean(cover, na.rm = TRUE),
    total_cover         = sum(cover, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(overstory_summary,  file.path(explore_dir, "D_overstory_summary_by_period.csv"))
write_csv(understory_summary, file.path(explore_dir, "D_understory_summary_by_period.csv"))
save_table_png(overstory_summary,  "D_overstory_summary_by_period.png",  "Overstory summary by period")
save_table_png(understory_summary, "D_understory_summary_by_period.png", "Understory summary by period")

# Differences (present - past), per species
overstory_diff <- overstory_long |>
  group_by(species, period) |>
  summarise(total_ba_m2 = sum(basal_area_spp_m2, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = total_ba_m2, values_fill = 0) |>
  mutate(diff_present_minus_past = present - past) |>
  arrange(desc(abs(diff_present_minus_past)))

understory_diff <- understory_long |>
  group_by(species, period) |>
  summarise(total_cover = sum(cover, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = total_cover, values_fill = 0) |>
  mutate(diff_present_minus_past = present - past) |>
  arrange(desc(abs(diff_present_minus_past)))

write_csv(overstory_diff,  file.path(explore_dir, "D_overstory_diff_by_species.csv"))
write_csv(understory_diff, file.path(explore_dir, "D_understory_diff_by_species.csv"))

pD1 <- overstory_diff |>
  slice_max(abs(diff_present_minus_past), n = 15) |>
  mutate(species = fct_reorder(species, diff_present_minus_past)) |>
  ggplot(aes(x = diff_present_minus_past, y = species, fill = diff_present_minus_past > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = col_present, "FALSE" = col_past), guide = "none") +
  labs(title = "Overstory: top 15 species by change in basal area",
       subtitle = "Present minus past (m2), largest absolute change",
       x = "Change in total basal area (m2)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic"))

pD2 <- understory_diff |>
  slice_max(abs(diff_present_minus_past), n = 15) |>
  mutate(species = fct_reorder(species, diff_present_minus_past)) |>
  ggplot(aes(x = diff_present_minus_past, y = species, fill = diff_present_minus_past > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = col_present, "FALSE" = col_past), guide = "none") +
  labs(title = "Understory: top 15 species by change in cover",
       subtitle = "Present minus past (%cover), largest absolute change",
       x = "Change in total cover (%)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic"))

save_plot(pD1, "D_overstory_diff_top15", width = 9, height = 7)
save_plot(pD2, "D_understory_diff_top15", width = 9, height = 7)


# ==============================================================================
# SECTION E — Plot coverage by resurvey year (present) and by elevation
# (past, used as a fixed per-plot proxy - present has no elevation source)
# ==============================================================================

plots_by_year <- plot_env_by_year |>
  filter(period == "present", plot_id %in% unique(c(overstory_long$plot_id, understory_long$plot_id))) |>
  distinct(plot_id, year) |>
  count(year, name = "n_plots")

write_csv(plots_by_year, file.path(explore_dir, "E_plots_by_year.csv"))
save_table_png(plots_by_year, "E_plots_by_year.png", "Plots by resurvey year (present)")

pE1 <- plots_by_year |>
  ggplot(aes(x = factor(year), y = n_plots)) +
  geom_col(fill = col_present, alpha = 0.9) +
  geom_text(aes(label = n_plots), vjust = -0.4, size = 3.5) +
  labs(title = "Plots surveyed by resurvey year", x = "Year", y = "Number of plots") +
  theme_minimal(base_size = 11)

save_plot(pE1, "E_plots_by_year", width = 7, height = 5)

plot_elevation <- plot_info |>
  filter(period == "past", !is.na(elevation)) |>
  distinct(plot_id, elevation) |>
  mutate(elev_bin = cut(elevation, breaks = seq(0, 2000, by = 200),
                        labels = paste0(seq(0, 1800, by = 200), "-", seq(200, 2000, by = 200), "m")))

plots_by_elevation <- plot_elevation |>
  filter(plot_id %in% unique(c(overstory_long$plot_id, understory_long$plot_id))) |>
  count(elev_bin, name = "n_plots")

write_csv(plots_by_elevation, file.path(explore_dir, "E_plots_by_elevation.csv"))
save_table_png(plots_by_elevation, "E_plots_by_elevation.png", "Plots by elevation band (past elevation, fixed proxy)")

pE2 <- plots_by_elevation |>
  filter(!is.na(elev_bin)) |>
  ggplot(aes(x = elev_bin, y = n_plots)) +
  geom_col(fill = col_both, alpha = 0.9) +
  geom_text(aes(label = n_plots), vjust = -0.4, size = 3.5) +
  labs(title = "Plots by elevation band", subtitle = "Elevation from 1975 survey, used as fixed per-plot proxy",
       x = "Elevation band", y = "Number of plots") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(pE2, "E_plots_by_elevation", width = 9, height = 5)


# ==============================================================================
# SECTION F — Turnover: colonization/extinction/persistence between
# periods, at the plot x species level. Functions defined inline below
# (following the earlier MORA and NOCA cleaning scripts; see header).
# ==============================================================================

calculate_turnover <- function(occ_data) {
  
  occ_wide <- occ_data |>
    select(plot_id, species, period, occurrence) |>
    distinct(plot_id, species, period, .keep_all = TRUE) |>
    pivot_wider(names_from = period, values_from = occurrence, values_fill = 0)
  
  turnover_data <- occ_wide |>
    mutate(
      past    = if_else(past    > 0, 1, 0),
      present = if_else(present > 0, 1, 0),
      status = case_when(
        past == 1 & present == 1 ~ "persistence",
        past == 1 & present == 0 ~ "extinction",
        past == 0 & present == 1 ~ "colonization",
        TRUE                      ~ "absent_both"
      )
    ) |>
    filter(status != "absent_both")
  
  summary_counts <- turnover_data |> count(status, name = "n")
  
  species_summary <- turnover_data |>
    count(species, status) |>
    pivot_wider(names_from = status, values_from = n, values_fill = 0)
  
  list(turnover_data = turnover_data, summary = summary_counts, species_summary = species_summary)
}

print_turnover_summary <- function(turnover_results) {
  cat("\n=== Turnover summary ===\n")
  print(turnover_results$summary)
  total <- sum(turnover_results$summary$n)
  cat("\nTotal plot x species combinations with a turnover event:", total, "\n")
  for (i in seq_len(nrow(turnover_results$summary))) {
    s <- turnover_results$summary$status[i]
    n <- turnover_results$summary$n[i]
    cat(sprintf("  %-14s: %5d (%.1f%%)\n", s, n, 100 * n / total))
  }
}

create_turnover_heatmap <- function(turnover_data, top_n_species = 40) {
  
  species_order <- turnover_data |>
    count(species) |>
    slice_max(n, n = top_n_species) |>
    arrange(n) |>
    pull(species)
  
  turnover_data |>
    filter(species %in% species_order) |>
    count(species, status) |>
    mutate(species = factor(species, levels = species_order)) |>
    ggplot(aes(x = status, y = species, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), size = 3, color = "white") +
    scale_fill_viridis_c(name = "N plots") +
    labs(subtitle = paste("Top", length(species_order), "species by total occurrence"), x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(face = "italic", size = 8))
}

create_alluvial_plot <- function(turnover_data,
                                 fill_colors = c(colonization = "#228B22",
                                                 extinction   = "#DC143C",
                                                 persistence  = col_both)) {
  
  alluvial_data <- turnover_data |>
    mutate(
      past_state    = if_else(past    == 1, "Present", "Absent"),
      present_state = if_else(present == 1, "Present", "Absent")
    ) |>
    count(past_state, present_state, status)
  
  ggplot(alluvial_data, aes(axis1 = past_state, axis2 = present_state, y = n)) +
    geom_alluvium(aes(fill = status), width = 1/4, alpha = 0.85, discern = TRUE) +
    geom_stratum(width = 1/4, fill = "grey90", color = "grey40", discern = TRUE) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_x_discrete(limits = c("Past (1975)", "Present"), expand = c(0.15, 0.05)) +
    scale_fill_manual(values = fill_colors, name = "Status") +
    labs(subtitle = "Plot x species combinations, past -> present",
         y = "Number of plot x species combinations") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
}

# Fallback if ggalluvial fails to render (known "could not find function
# gg_par" version-mismatch issue between ggalluvial and newer ggplot2) -
# same data (status counts) shown as a grouped bar chart instead of a
# flow diagram.
create_turnover_barchart_fallback <- function(turnover_data,
                                              fill_colors = c(colonization = "#228B22",
                                                              extinction   = "#DC143C",
                                                              persistence  = col_both)) {
  turnover_data |>
    count(status) |>
    mutate(status = factor(status, levels = c("colonization", "persistence", "extinction"))) |>
    ggplot(aes(x = status, y = n, fill = status)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n), vjust = -0.4, size = 4) +
    scale_fill_manual(values = fill_colors, guide = "none") +
    labs(subtitle = "Plot x species combinations by status (fallback chart - ggalluvial unavailable)",
         x = NULL, y = "Number of plot x species combinations") +
    theme_minimal(base_size = 12)
}

# --- Overstory ---
cat("\n########## OVERSTORY TURNOVER ##########\n")
over_turnover <- calculate_turnover(overstory_long)
print_turnover_summary(over_turnover)

write_csv(over_turnover$summary, file.path(explore_dir, "F_overstory_turnover_summary.csv"))
save_table_png(over_turnover$summary, "F_overstory_turnover_summary.png", "Overstory turnover summary")

# Bar chart of status counts (instead of the alluvial diagram; see header).
over_status_bar <- create_turnover_barchart_fallback(over_turnover$turnover_data) +
  labs(title = "Overstory: colonization and extinction")
over_heatmap  <- create_turnover_heatmap(over_turnover$turnover_data, top_n_species = 23) +
  labs(title = "Overstory: turnover by species")

save_plot(over_status_bar, "F_overstory_status_bar", width = 6, height = 6)
save_plot(over_heatmap, "F_overstory_heatmap", width = 7, height = 8)

# --- Understory ---
cat("\n########## UNDERSTORY TURNOVER ##########\n")
under_turnover <- calculate_turnover(understory_long)
print_turnover_summary(under_turnover)

write_csv(under_turnover$summary, file.path(explore_dir, "F_understory_turnover_summary.csv"))
save_table_png(under_turnover$summary, "F_understory_turnover_summary.png", "Understory turnover summary")

under_status_bar <- create_turnover_barchart_fallback(under_turnover$turnover_data) +
  labs(title = "Understory: colonization and extinction")
under_heatmap  <- create_turnover_heatmap(under_turnover$turnover_data, top_n_species = 40) +
  labs(title = "Understory: turnover by species (top 40)")

save_plot(under_status_bar, "F_understory_status_bar", width = 6, height = 6)
save_plot(under_heatmap, "F_understory_heatmap", width = 7, height = 10)

# --- Combined PDF ---
pdf(file.path(explore_dir, "F_turnover_dynamics.pdf"), width = 10, height = 8, onefile = TRUE)
print(over_status_bar)
print(over_heatmap)
print(under_status_bar)
print(under_heatmap)
dev.off()


# ==============================================================================
# SECTION G — Rank change: how each species' abundance RANK shifts from
# past to present (bump/slope chart). Overstory ranked by total basal
# area (m2); understory ranked by total cover (%). Understory limited to
# the top N species by combined abundance - a bump chart with 221 species
# would be unreadable.
# ==============================================================================

make_rank_chart <- function(long_data, value_col, top_n = NULL, title, y_lab) {
  
  totals <- long_data |>
    group_by(species, period) |>
    summarise(total_value = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  
  if (!is.null(top_n)) {
    keep_species <- totals |>
      group_by(species) |>
      summarise(combined = sum(total_value), .groups = "drop") |>
      slice_max(combined, n = top_n) |>
      pull(species)
    totals <- totals |> filter(species %in% keep_species)
  }
  
  ranked <- totals |>
    group_by(period) |>
    mutate(rank = rank(-total_value, ties.method = "min")) |>
    ungroup() |>
    mutate(period = factor(period, levels = c("past", "present")))
  
  ggplot(ranked, aes(x = period, y = rank, group = species)) +
    geom_line(aes(color = species), linewidth = 0.8, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(color = species), size = 2.5, show.legend = FALSE) +
    geom_text(data = ranked |> filter(period == "past"),
              aes(label = species), hjust = 1.1, size = 2.8, fontface = "italic") +
    geom_text(data = ranked |> filter(period == "present"),
              aes(label = species), hjust = -0.1, size = 2.8, fontface = "italic") +
    scale_y_reverse(breaks = scales::pretty_breaks()) +
    scale_x_discrete(expand = expansion(mult = c(0.6, 0.6))) +
    labs(title = title, subtitle = paste("Ranked by", y_lab), x = NULL, y = "Rank (1 = highest)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), axis.text.y = element_blank())
}

pG1 <- make_rank_chart(overstory_long, "basal_area_spp_m2", top_n = NULL,
                       title = "Overstory: species rank shift (past -> present)",
                       y_lab = "total basal area (m2)")

pG2 <- make_rank_chart(understory_long, "cover", top_n = 25,
                       title = "Understory: species rank shift (past -> present)",
                       y_lab = "total %cover, top 25 species")

save_plot(pG1, "G_overstory_rank_shift",  width = 9,  height = 8)
save_plot(pG2, "G_understory_rank_shift", width = 10, height = 10)


# ==============================================================================
# SECTION H — Elevation shift: mean elevation of occurrence per species,
# past vs present. Elevation from plot_info.csv (past survey only, used
# as a fixed per-plot attribute). Weighted by presence only
# (occurrence == 1), not by every plot in the dataset - answers
# "are species occurring at higher elevations now than before?"
# ==============================================================================

plot_elev_lookup <- plot_info |>
  filter(period == "past", !is.na(elevation)) |>
  distinct(plot_id, elevation)

make_elevation_shift <- function(long_data, min_occurrences = 5, top_n = NULL, title) {
  
  elev_by_species_period <- long_data |>
    filter(occurrence == 1) |>
    left_join(plot_elev_lookup, by = "plot_id") |>
    filter(!is.na(elevation)) |>
    group_by(species, period) |>
    summarise(mean_elev = mean(elevation, na.rm = TRUE), n_plots = n(), .groups = "drop")
  
  wide <- elev_by_species_period |>
    select(species, period, mean_elev, n_plots) |>
    pivot_wider(names_from = period, values_from = c(mean_elev, n_plots)) |>
    filter(n_plots_past >= min_occurrences, n_plots_present >= min_occurrences) |>
    mutate(shift_m = mean_elev_present - mean_elev_past)
  
  if (!is.null(top_n)) {
    wide <- wide |> slice_max(abs(shift_m), n = top_n)
  }
  
  wide <- wide |> mutate(species = fct_reorder(species, shift_m))
  
  ggplot(wide) +
    geom_segment(aes(x = mean_elev_past, xend = mean_elev_present, y = species, yend = species),
                 color = "grey60", linewidth = 0.8) +
    geom_point(aes(x = mean_elev_past, y = species), color = col_past, size = 3) +
    geom_point(aes(x = mean_elev_present, y = species), color = col_present, size = 3) +
    labs(title = title,
         subtitle = paste0("Mean elevation of occurrence, past (purple) -> present (orange). ",
                           "Min ", min_occurrences, " plots/period"),
         x = "Elevation (m)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(face = "italic", size = 8))
}

pH1 <- make_elevation_shift(overstory_long, min_occurrences = 5, top_n = NULL,
                            title = "Overstory: elevation shift by species")
pH2 <- make_elevation_shift(understory_long, min_occurrences = 10, top_n = 40,
                            title = "Understory: elevation shift by species (top 40 by |shift|)")

save_plot(pH1, "H_overstory_elevation_shift",  width = 9,  height = 7)
save_plot(pH2, "H_understory_elevation_shift", width = 9,  height = 10)

# Also save the underlying numbers as CSV, for both layers, unrestricted
overstory_elev_shift  <- overstory_long |>
  filter(occurrence == 1) |> left_join(plot_elev_lookup, by = "plot_id") |>
  filter(!is.na(elevation)) |> group_by(species, period) |>
  summarise(mean_elev = mean(elevation, na.rm = TRUE), n_plots = n(), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = c(mean_elev, n_plots)) |>
  mutate(shift_m = mean_elev_present - mean_elev_past) |>
  arrange(desc(abs(shift_m)))

understory_elev_shift <- understory_long |>
  filter(occurrence == 1) |> left_join(plot_elev_lookup, by = "plot_id") |>
  filter(!is.na(elevation)) |> group_by(species, period) |>
  summarise(mean_elev = mean(elevation, na.rm = TRUE), n_plots = n(), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = c(mean_elev, n_plots)) |>
  mutate(shift_m = mean_elev_present - mean_elev_past) |>
  arrange(desc(abs(shift_m)))

write_csv(overstory_elev_shift,  file.path(explore_dir, "H_overstory_elevation_shift.csv"))
write_csv(understory_elev_shift, file.path(explore_dir, "H_understory_elevation_shift.csv"))

cat("\nMean elevation shift (present - past), overstory:\n")
cat("  Mean shift across species:", round(mean(overstory_elev_shift$shift_m, na.rm = TRUE), 1), "m\n")
cat("  Species shifting upward  :", sum(overstory_elev_shift$shift_m > 0, na.rm = TRUE), "\n")
cat("  Species shifting downward:", sum(overstory_elev_shift$shift_m < 0, na.rm = TRUE), "\n")

cat("\nMean elevation shift (present - past), understory:\n")
cat("  Mean shift across species:", round(mean(understory_elev_shift$shift_m, na.rm = TRUE), 1), "m\n")
cat("  Species shifting upward  :", sum(understory_elev_shift$shift_m > 0, na.rm = TRUE), "\n")
cat("  Species shifting downward:", sum(understory_elev_shift$shift_m < 0, na.rm = TRUE), "\n")


# ==============================================================================
# SECTION I — Elevation DISTRIBUTION per species (not just the mean):
# min/max/mean/median/sd/skewness/kurtosis, past vs present, plus a small-
# multiples density plot per species. Kurtosis/skewness computed manually
# (no extra package needed - moments/e1071 would add another install
# dependency, and package installs have already been a headache today).
# ==============================================================================

# Excess kurtosis (0 = normal-shaped tails, >0 = heavier/pointier tails,
# <0 = flatter). Skewness >0 = long tail toward high elevation, <0 = long
# tail toward low elevation.
compute_shape_stats <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) return(tibble(n = n, min = NA, max = NA, mean = NA, median = NA,
                           sd = NA, skewness = NA, kurtosis = NA))
  m <- mean(x); s <- sd(x)
  tibble(
    n = n, min = min(x), max = max(x), mean = m, median = median(x), sd = s,
    skewness = (sum((x - m)^3) / n) / s^3,
    kurtosis = (sum((x - m)^4) / n) / s^4 - 3
  )
}

get_elev_occurrences <- function(long_data, group_col) {
  long_data |>
    filter(occurrence == 1) |>
    left_join(plot_elev_lookup, by = "plot_id") |>
    filter(!is.na(elevation)) |>
    select(group = all_of(group_col), period, elevation)
}

make_distribution_stats <- function(occ_elev, min_n = 4) {
  occ_elev |>
    group_by(group, period) |>
    group_modify(~ compute_shape_stats(.x$elevation)) |>
    ungroup() |>
    filter(n >= min_n)
}

# --- Overstory: all 23 species ---
over_occ_elev   <- get_elev_occurrences(overstory_long, "species")
over_elev_stats <- make_distribution_stats(over_occ_elev) |> rename(species = group)

write_csv(over_elev_stats, file.path(explore_dir, "I_overstory_elevation_stats.csv"))
save_table_png(over_elev_stats |> mutate(across(where(is.numeric), ~round(., 1))),
               "I_overstory_elevation_stats.png", "Overstory elevation distribution stats")

pI1 <- over_occ_elev |>
  rename(species = group) |>
  ggplot(aes(x = elevation, fill = period, color = period)) +
  geom_density(alpha = 0.4, linewidth = 0.6) +
  geom_vline(data = over_elev_stats, aes(xintercept = mean, color = period),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  facet_wrap(~ species, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(past = col_past, present = col_present)) +
  scale_color_manual(values = c(past = col_past, present = col_present)) +
  labs(title = "Overstory: elevation distribution by species",
       subtitle = "Dashed line = mean elevation. Density of plots with real presence.",
       x = "Elevation (m)", y = "Density", fill = NULL, color = NULL) +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "italic", size = 8), legend.position = "top")

save_plot(pI1, "I_overstory_elevation_distribution", width = 13, height = 14)

# --- Understory: top 24 species by total occurrence (small multiples of
# 221 species would be unreadable - the life-form view in Section J below
# covers the full set at a coarser grouping).
under_occ_elev <- get_elev_occurrences(understory_long, "species")

top_understory_species <- understory_long |>
  filter(occurrence == 1) |>
  count(species, sort = TRUE) |>
  slice_max(n, n = 24) |>
  pull(species)

under_elev_stats_top <- make_distribution_stats(under_occ_elev |> filter(group %in% top_understory_species)) |>
  rename(species = group)

write_csv(under_elev_stats_top, file.path(explore_dir, "I_understory_elevation_stats_top24.csv"))
save_table_png(under_elev_stats_top |> mutate(across(where(is.numeric), ~round(., 1))),
               "I_understory_elevation_stats_top24.png", "Understory elevation stats (top 24 species)")

# Full stats for ALL understory species too, CSV only (too many for a table image)
under_elev_stats_all <- make_distribution_stats(under_occ_elev) |> rename(species = group)
write_csv(under_elev_stats_all, file.path(explore_dir, "I_understory_elevation_stats_all.csv"))

pI2 <- under_occ_elev |>
  filter(group %in% top_understory_species) |>
  rename(species = group) |>
  ggplot(aes(x = elevation, fill = period, color = period)) +
  geom_density(alpha = 0.4, linewidth = 0.6) +
  geom_vline(data = under_elev_stats_top, aes(xintercept = mean, color = period),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  facet_wrap(~ species, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(past = col_past, present = col_present)) +
  scale_color_manual(values = c(past = col_past, present = col_present)) +
  labs(title = "Understory: elevation distribution by species (top 24 by occurrence)",
       subtitle = "Dashed line = mean elevation. Density of plots with real presence.",
       x = "Elevation (m)", y = "Density", fill = NULL, color = NULL) +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "italic", size = 8), legend.position = "top")

save_plot(pI2, "I_understory_elevation_distribution_top24", width = 13, height = 14)


# ==============================================================================
# SECTION J — Understory grouped by LIFE FORM (Forb/Shrub/Grass/Fern/
# Sedge/Rush), instead of by individual species - covers ALL 221 species
# pooled into their life-form category, past vs present.
# ==============================================================================

spp_info <- read_csv(file.path(cleaning_dir, "spp_info.csv"), show_col_types = FALSE)
life_form_lookup <- spp_info |> select(species, life_form) |> filter(!is.na(life_form))

under_occ_elev_lifeform <- understory_long |>
  filter(occurrence == 1) |>
  left_join(plot_elev_lookup, by = "plot_id") |>
  filter(!is.na(elevation)) |>
  left_join(life_form_lookup, by = "species") |>
  filter(!is.na(life_form)) |>
  select(group = life_form, period, elevation)

lifeform_elev_stats <- make_distribution_stats(under_occ_elev_lifeform, min_n = 4) |>
  rename(life_form = group)

write_csv(lifeform_elev_stats, file.path(explore_dir, "J_understory_lifeform_elevation_stats.csv"))
save_table_png(lifeform_elev_stats |> mutate(across(where(is.numeric), ~round(., 1))),
               "J_understory_lifeform_elevation_stats.png", "Understory elevation stats by life form")

pJ1 <- under_occ_elev_lifeform |>
  rename(life_form = group) |>
  ggplot(aes(x = elevation, fill = period, color = period)) +
  geom_density(alpha = 0.45, linewidth = 0.8) +
  geom_vline(data = lifeform_elev_stats, aes(xintercept = mean, color = period),
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  facet_wrap(~ life_form, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(past = col_past, present = col_present)) +
  scale_color_manual(values = c(past = col_past, present = col_present)) +
  labs(title = "Understory: elevation distribution by life form",
       subtitle = "All species pooled within each life form. Dashed line = mean elevation.",
       x = "Elevation (m)", y = "Density", fill = NULL, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 10), legend.position = "top")

save_plot(pJ1, "J_understory_lifeform_elevation_distribution", width = 11, height = 8)

# Life-form mean-elevation shift, dumbbell style (mirrors Section H, but
# at the life-form level - much easier to read at a glance)
lifeform_wide <- lifeform_elev_stats |>
  select(life_form, period, mean) |>
  mutate(life_form = as.character(life_form)) |>
  # Guard against any duplicate life_form x period rows (would otherwise
  # break pivot_wider's one-value-per-cell assumption)
  group_by(life_form, period) |>
  summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = mean) |>
  mutate(diff = present - past)

lifeform_order <- lifeform_wide |> arrange(diff) |> pull(life_form)

pJ2 <- lifeform_wide |>
  mutate(life_form = factor(life_form, levels = lifeform_order)) |>
  ggplot() +
  geom_segment(aes(x = past, xend = present, y = life_form, yend = life_form), color = "grey60", linewidth = 1) +
  geom_point(aes(x = past, y = life_form), color = col_past, size = 4) +
  geom_point(aes(x = present, y = life_form), color = col_present, size = 4) +
  labs(title = "Understory: elevation shift by life form",
       subtitle = "Mean elevation of occurrence, past (purple) -> present (orange)",
       x = "Elevation (m)", y = NULL) +
  theme_minimal(base_size = 12)

save_plot(pJ2, "J_understory_lifeform_elevation_shift", width = 8, height = 5)

cat("\nUnderstory elevation shift by life form (present - past):\n")
print(lifeform_wide |> mutate(shift_m = round(diff, 1)) |> select(life_form, past, present, shift_m) |>
        arrange(desc(abs(shift_m))))

cat("\nAll exploration outputs saved to:", explore_dir, "\n")



# ==============================================================================
# SECTION K — Turnover (colonization/extinction/persistence, from Section
# F) BY ELEVATION BAND. Tests the climate-shift hypothesis directly: more
# colonization at high elevation, more extinction at low elevation?
# Reuses over_turnover/under_turnover from Section F and plot_elev_lookup
# from Section H (both already in the session by this point).
# ==============================================================================

elev_bands <- tibble(
  elev_bin = factor(paste0(seq(0, 1800, by = 200), "-", seq(200, 2000, by = 200), "m"),
                    levels = paste0(seq(0, 1800, by = 200), "-", seq(200, 2000, by = 200), "m"))
)

add_elevation_band <- function(turnover_data) {
  turnover_data |>
    left_join(plot_elev_lookup, by = "plot_id") |>
    filter(!is.na(elevation)) |>
    mutate(elev_bin = cut(elevation, breaks = seq(0, 2000, by = 200),
                          labels = levels(elev_bands$elev_bin), include.lowest = TRUE))
}

make_turnover_by_elevation <- function(turnover_data) {
  add_elevation_band(turnover_data) |>
    count(elev_bin, status) |>
    group_by(elev_bin) |>
    mutate(total = sum(n), rate = n / total) |>
    ungroup() |>
    filter(!is.na(elev_bin))
}

over_turnover_elev  <- make_turnover_by_elevation(over_turnover$turnover_data)
under_turnover_elev <- make_turnover_by_elevation(under_turnover$turnover_data)

write_csv(over_turnover_elev,  file.path(explore_dir, "K_overstory_turnover_by_elevation.csv"))
write_csv(under_turnover_elev, file.path(explore_dir, "K_understory_turnover_by_elevation.csv"))

status_colors <- c(colonization = "#228B22", extinction = "#DC143C", persistence = col_both)

make_turnover_elev_barplot <- function(turnover_elev, title) {
  turnover_elev |>
    mutate(status = factor(status, levels = c("extinction", "persistence", "colonization"))) |>
    ggplot(aes(x = elev_bin, y = rate, fill = status)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = status_colors, name = NULL) +
    labs(title = title, subtitle = "Proportion of plot x species combinations by status, per elevation band",
         x = "Elevation band", y = "Proportion") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
}

pK1 <- make_turnover_elev_barplot(over_turnover_elev,  "Overstory: turnover by elevation band")
pK2 <- make_turnover_elev_barplot(under_turnover_elev, "Understory: turnover by elevation band")

save_plot(pK1, "K_overstory_turnover_by_elevation",  width = 9, height = 6)
save_plot(pK2, "K_understory_turnover_by_elevation", width = 9, height = 6)

# Net turnover (colonization rate - extinction rate) per elevation band -
# directly tests the hypothesis: positive = net gain (colonization-driven,
# expected at high elevation under a climate-shift hypothesis), negative =
# net loss (extinction-driven, expected at low elevation).
make_net_turnover <- function(turnover_elev) {
  turnover_elev |>
    filter(status %in% c("colonization", "extinction")) |>
    select(elev_bin, status, rate) |>
    pivot_wider(names_from = status, values_from = rate, values_fill = 0) |>
    mutate(net_turnover = colonization - extinction)
}

over_net_turnover  <- make_net_turnover(over_turnover_elev)
under_net_turnover <- make_net_turnover(under_turnover_elev)

write_csv(over_net_turnover,  file.path(explore_dir, "K_overstory_net_turnover_by_elevation.csv"))
write_csv(under_net_turnover, file.path(explore_dir, "K_understory_net_turnover_by_elevation.csv"))

make_net_turnover_plot <- function(net_turnover, title) {
  net_turnover |>
    ggplot(aes(x = elev_bin, y = net_turnover, fill = net_turnover > 0)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
    scale_fill_manual(values = c("TRUE" = "#228B22", "FALSE" = "#DC143C"), guide = "none") +
    labs(title = title,
         subtitle = "Colonization rate minus extinction rate. Positive = net gain, negative = net loss.",
         x = "Elevation band", y = "Net turnover (colonization - extinction)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

pK3 <- make_net_turnover_plot(over_net_turnover,  "Overstory: net turnover by elevation band")
pK4 <- make_net_turnover_plot(under_net_turnover, "Understory: net turnover by elevation band")

save_plot(pK3, "K_overstory_net_turnover_by_elevation",  width = 8, height = 5)
save_plot(pK4, "K_understory_net_turnover_by_elevation", width = 8, height = 5)

cat("\nOverstory net turnover by elevation band (positive = net colonization gain):\n")
print(over_net_turnover |> select(elev_bin, colonization, extinction, net_turnover))
cat("\nUnderstory net turnover by elevation band:\n")
print(under_net_turnover |> select(elev_bin, colonization, extinction, net_turnover))


# ==============================================================================
# SECTION L — Overstory density/cover vs elevation, past vs present
# ("treeline advance" plot). Total basal area per plot, aggregated across
# species, plotted against elevation with a smoothed trend line per
# period - is the forest getting denser at higher elevations now?
# ==============================================================================

overstory_plot_totals <- overstory_long |>
  group_by(period, plot_id) |>
  summarise(total_ba_perc = sum(basal_area_spp_perc, na.rm = TRUE),
            total_count   = sum(count, na.rm = TRUE), .groups = "drop") |>
  left_join(plot_elev_lookup, by = "plot_id") |>
  filter(!is.na(elevation))

write_csv(overstory_plot_totals, file.path(explore_dir, "L_overstory_cover_vs_elevation.csv"))

pL1 <- overstory_plot_totals |>
  ggplot(aes(x = elevation, y = total_ba_perc, color = period, fill = period)) +
  geom_point(alpha = 0.35, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2, alpha = 0.2) +
  scale_color_manual(values = c(past = col_past, present = col_present)) +
  scale_fill_manual(values = c(past = col_past, present = col_present)) +
  labs(title = "Overstory: total basal area vs elevation",
       subtitle = "Is the forest getting denser at higher elevations? Loess smooth with 95% CI",
       x = "Elevation (m)", y = "Total basal area (% of plot area)", color = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

save_plot(pL1, "L_overstory_basal_area_vs_elevation", width = 9, height = 6)

# Same idea, binned by elevation band instead of a continuous smooth -
# easier to read exact values off, matches Section K/E's banding.
overstory_totals_by_band <- overstory_plot_totals |>
  mutate(elev_bin = cut(elevation, breaks = seq(0, 2000, by = 200),
                        labels = levels(elev_bands$elev_bin), include.lowest = TRUE)) |>
  filter(!is.na(elev_bin)) |>
  group_by(period, elev_bin) |>
  summarise(mean_ba_perc = mean(total_ba_perc, na.rm = TRUE),
            mean_count   = mean(total_count, na.rm = TRUE),
            n_plots      = n(), .groups = "drop")

write_csv(overstory_totals_by_band, file.path(explore_dir, "L_overstory_cover_by_elevation_band.csv"))
save_table_png(overstory_totals_by_band |> mutate(across(where(is.numeric), ~round(., 1))),
               "L_overstory_cover_by_elevation_band.png", "Overstory: mean cover by elevation band")

pL2 <- overstory_totals_by_band |>
  ggplot(aes(x = elev_bin, y = mean_ba_perc, color = period, group = period)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c(past = col_past, present = col_present)) +
  labs(title = "Overstory: mean basal area by elevation band",
       subtitle = "Mean total basal area (% of plot area) per plot, by elevation band",
       x = "Elevation band", y = "Mean total basal area (%)", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

save_plot(pL2, "L_overstory_cover_by_elevation_band", width = 9, height = 6)

cat("\nOverstory mean basal area (%) by elevation band and period:\n")
print(overstory_totals_by_band |> select(period, elev_bin, mean_ba_perc, n_plots) |>
        pivot_wider(names_from = period, values_from = c(mean_ba_perc, n_plots)))

cat("\nAll exploration outputs saved to:", explore_dir, "\n")