# ------------------------------------------------------------------------------
# 01_climate_trends.R
# Pipeline: pipeline/01_data_prep/03_climate/02_explore_climate/code/
# Author: Agostina Torres 
# ------------------------------------------------------------------------------
# Purpose: describe how climate at the MORA plots has changed over the full
# CHELSA record (1979-2021), and how the environmental variables relate to
# each other and to elevation, using the annual series and past/present
# window summary produced by 02_compute_annual_variables.R.
#
# Produces two deliverables:
#   1. A temporal trend figure: the full annual series (across-plot mean +/-
#      1 SE per year, with a smoothed trend line) for MAT and MAP, each paired
#      with a "zoom" panel comparing the past (1979-1983) and present
#      (2013-2017) window summaries directly (one point per plot per window).
#   2. A PCA of the environmental variables (MAT, MAT_min, MAP, PAS, seasonal
#      temperature and precipitation, and elevation), to see how the
#      variables covary across plots and which ones carry the most distinct
#      information -- useful both for describing the site's environmental
#      gradient and for trimming redundant candidate HMSC covariates later.
#      A second panel projects the past-window values into the same PCA
#      space to show whether the plots' climate has shifted position in
#      that multivariate space between periods.
#
# Only the plots used in the analyses are included (the 324 plots of the
# overstory climatic-affinity analyses, which contain the 316 understory
# plots), so that the climate summaries reported in the Results describe the
# same plots as the community analyses.
#
# Aspect (a circular variable) is not included in the PCA.
#
# Output:
#   - figures/climate_trend_full_and_zoom.png: MAT/MAP annual trend + past-
#     present zoom panels, combined
#   - figures/climate_pca_biplot.png: PCA biplot (present-day variables +
#     elevation, one point per plot) with loading arrows
#   - figures/climate_pca_past_present_shift.png: past vs. present plot
#     scores projected into the same PCA space
#   - outputs/climate_pca_variance_explained.csv: variance explained per PC
#   - outputs/climate_pca_loadings.csv: variable loadings per PC
#   - outputs/climate_window_summary.csv: mean past and present climate across
#     plots and their difference (values reported in the Results)
# ------------------------------------------------------------------------------

# 0. Libraries -------------------------------------------------------------------
library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)
library(ggrepel)

# 0a. Module directory --------------------------------------------------------
this_module_dir <- function(...) here("pipeline", "01_data_prep", "03_climate", "02_explore_climate", ...)

dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)

# Shared plot theme, matching the climatic-affinity figures elsewhere in the
# project: no gridlines, thin black axis lines, viridis-family colors.
climate_theme <- theme_ipsum(grid = FALSE) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    axis.line.y = element_line(color = "black", linewidth = 0.3)
  )

# 1. Load data -------------------------------------------------------------------
climate_annual_long <- read_csv(this_module_dir("data_clean", "climate_annual_long.csv"),
                                show_col_types = FALSE)
climate_past_present_wide <- read_csv(this_module_dir("data_clean", "climate_past_present_wide.csv"),
                                      show_col_types = FALSE)

# Plots used in the analyses
analysis_plots <- read_csv(here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", "outputs",
                                "overstory_thermo_xero_rates.csv"), show_col_types = FALSE) %>%
  distinct(plot_id = as.character(plot_id)) %>% pull(plot_id)
climate_annual_long <- climate_annual_long %>%
  mutate(plot_id = as.character(plot_id)) %>% filter(plot_id %in% analysis_plots)
climate_past_present_wide <- climate_past_present_wide %>%
  mutate(plot_id = as.character(plot_id)) %>% filter(plot_id %in% analysis_plots)
cat("Plots used:", n_distinct(climate_past_present_wide$plot_id), "of", length(analysis_plots), "analysis plots\n")

# Precipitation is a sum of monthly values, so a plot without climate data must
# be NA, never 0; stop if a zero is found
if (any(climate_past_present_wide$MAP_past == 0 | climate_past_present_wide$MAP_present == 0, na.rm = TRUE)) {
  stop("MAP = 0 found in climate_past_present_wide.csv: rerun 02_compute_annual_variables.R")
}

# Mean climate across plots in each window and their difference (Results)
window_summary <- tibble(variable = c("MAT", "MAP", "PAS", "temp_JJA")) %>%
  filter(paste0(variable, "_past") %in% names(climate_past_present_wide)) %>%
  mutate(past    = map_dbl(variable, ~ mean(climate_past_present_wide[[paste0(.x, "_past")]], na.rm = TRUE)),
         present = map_dbl(variable, ~ mean(climate_past_present_wide[[paste0(.x, "_present")]], na.rm = TRUE)),
         change  = present - past)
write_csv(mutate(window_summary, across(where(is.double), ~ signif(.x, 4))),
          this_module_dir("outputs", "climate_window_summary.csv"))
cat("\nMean climate across plots, past and present windows:\n")
print(mutate(window_summary, across(where(is.double), ~ signif(.x, 4))))

# ================================================================================
# 2. TEMPORAL TREND: full annual series + past/present zoom
# ================================================================================

# 2a. Across-plot annual summary (mean +/- SE per year), for the full-series
# panel. na.rm handles both the 1979 tas/tasmin gap and the plots with no
# coordinates (NA in every year) without needing to filter them out by name.
annual_summary <- climate_annual_long %>%
  group_by(year) %>%
  summarise(
    MAT_mean = mean(MAT, na.rm = TRUE),
    MAT_se   = sd(MAT, na.rm = TRUE) / sqrt(sum(!is.na(MAT))),
    MAP_mean = mean(MAP, na.rm = TRUE),
    MAP_se   = sd(MAP, na.rm = TRUE) / sqrt(sum(!is.na(MAP))),
    .groups = "drop"
  )

# 2b. Past/present window values, one row per plot per window, for the zoom
# panel -- reuses the already-computed per-plot window means rather than
# recomputing them, so it stays consistent with climate_past_present_wide.csv.
window_long <- bind_rows(
  climate_past_present_wide %>%
    transmute(plot_id, window = "past", MAT = MAT_past, MAP = MAP_past),
  climate_past_present_wide %>%
    transmute(plot_id, window = "present", MAT = MAT_present, MAP = MAP_present)
) %>%
  mutate(window = factor(window, levels = c("past", "present")))

window_color_scale <- scale_color_viridis_d(option = "viridis", begin = 0.15, end = 0.75, name = NULL)
window_fill_scale  <- scale_fill_viridis_d(option = "viridis", begin = 0.15, end = 0.75, name = NULL)

# 2c. Full-series trend panels (MAT, MAP) -----------------------------------------
p_trend_mat <- ggplot(annual_summary, aes(x = year, y = MAT_mean)) +
  geom_ribbon(aes(ymin = MAT_mean - MAT_se, ymax = MAT_mean + MAT_se),
              fill = "grey70", alpha = 0.4) +
  geom_line(color = "grey30", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "#443A83", linewidth = 1) +
  labs(x = "Year", y = "Mean annual temperature (°C)") +
  climate_theme

p_trend_map <- ggplot(annual_summary, aes(x = year, y = MAP_mean)) +
  geom_ribbon(aes(ymin = MAP_mean - MAP_se, ymax = MAP_mean + MAP_se),
              fill = "grey70", alpha = 0.4) +
  geom_line(color = "grey30", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "#21918C", linewidth = 1) +
  labs(x = "Year", y = "Annual precipitation (mm)") +
  climate_theme

# 2d. Past vs. present zoom panels (one point per plot per window) ----------------
p_zoom_mat <- ggplot(window_long, aes(x = window, y = MAT, color = window, fill = window)) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  window_color_scale +
  window_fill_scale +
  labs(x = NULL, y = NULL) +
  climate_theme +
  theme(legend.position = "none")

p_zoom_map <- ggplot(window_long, aes(x = window, y = MAP, color = window, fill = window)) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  window_color_scale +
  window_fill_scale +
  labs(x = NULL, y = NULL) +
  climate_theme +
  theme(legend.position = "none")

# 2e. Combine: full trend (wide) + zoom (narrow) per variable, stacked --------------
p_climate_trend <- (p_trend_mat + p_zoom_mat + plot_layout(widths = c(3, 1))) /
  (p_trend_map + p_zoom_map + plot_layout(widths = c(3, 1)))

ggsave(this_module_dir("figures", "climate_trend_full_and_zoom.png"),
       p_climate_trend, width = 10, height = 8, dpi = 300, bg = "white")

# ================================================================================
# 3. PCA: how do the environmental variables relate to each other?
# ================================================================================

# 3a. Build the PCA input: one row per plot, present-day climate variables
# (the plots' current environment, as used downstream for HMSC) plus
# elevation. Plots missing any of these (the 9 without coordinates, or any
# plot missing elevation) are dropped from the PCA -- can't be projected in
# a space they have no coordinates for.
pca_vars <- c("MAT_present", "MAT_min_present", "MAP_present", "PAS_present",
              "temp_DJF_present", "temp_MAM_present", "temp_JJA_present", "temp_SON_present",
              "precip_DJF_present", "precip_MAM_present", "precip_JJA_present", "precip_SON_present",
              "elevation")

pca_input <- climate_past_present_wide %>%
  select(plot_id, all_of(pca_vars)) %>%
  drop_na()

cat("PCA: ", nrow(pca_input), "of", nrow(climate_past_present_wide),
    "plots with complete data for the variables", "\n  ", paste(pca_vars, collapse = ", "), "\n")

pca_matrix <- pca_input %>% select(-plot_id) %>% as.matrix()
rownames(pca_matrix) <- pca_input$plot_id

# 3b. Run PCA on standardized variables (different units: °C, mm, m) --------------
pca_fit <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)

var_explained <- tibble(
  PC = paste0("PC", seq_along(pca_fit$sdev)),
  variance = pca_fit$sdev^2,
  pct_variance = 100 * variance / sum(variance)
) %>%
  mutate(pct_cumulative = cumsum(pct_variance))

write_csv(var_explained, this_module_dir("outputs", "climate_pca_variance_explained.csv"))

loadings <- as_tibble(pca_fit$rotation, rownames = "variable")
write_csv(loadings, this_module_dir("outputs", "climate_pca_loadings.csv"))

cat("\nVariance explained -- PC1:", round(var_explained$pct_variance[1], 1),
    "% | PC2:", round(var_explained$pct_variance[2], 1),
    "% | cumulative:", round(var_explained$pct_cumulative[2], 1), "%\n")

# 3c. Biplot: plot scores + variable loading arrows --------------------------------
pc1_lab <- paste0("PC1 (", round(var_explained$pct_variance[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained$pct_variance[2], 1), "%)")

scores <- as_tibble(pca_fit$x, rownames = "plot_id") %>%
  left_join(pca_input %>% select(plot_id, elevation), by = "plot_id")

# Scale loading arrows to a comparable range with the point cloud, a
# standard biplot convention (arrows otherwise sit on their own tiny scale).
arrow_scale <- 0.8 * max(abs(scores$PC1), abs(scores$PC2)) / max(abs(loadings$PC1), abs(loadings$PC2))
loadings_scaled <- loadings %>%
  mutate(PC1_arrow = PC1 * arrow_scale, PC2_arrow = PC2 * arrow_scale)

p_pca_biplot <- ggplot() +
  geom_point(data = scores, aes(x = PC1, y = PC2, color = elevation), size = 1.8, alpha = 0.7) +
  scale_color_viridis_c(option = "viridis", direction = -1, name = "Elevation (m)") +
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1_arrow, yend = PC2_arrow),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.4) +
  geom_text_repel(data = loadings_scaled, aes(x = PC1_arrow, y = PC2_arrow, label = variable),
                  size = 3, color = "black", max.overlaps = 20) +
  labs(x = pc1_lab, y = pc2_lab) +
  climate_theme

ggsave(this_module_dir("figures", "climate_pca_biplot.png"),
       p_pca_biplot, width = 8, height = 7, dpi = 300, bg = "white")

# ================================================================================
# 4. PCA PAST + PRESENT: climatic variables associated with temporal separation
# ================================================================================
# Fit a new PCA jointly to past and present climate. Elevation is NOT an active
# variable here because it does not change between periods. This PCA is used to
# identify climatic variables associated with the multivariate past -> present shift.

climate_stems <- c("MAT", "MAT_min", "MAP", "PAS",
                   "temp_DJF", "temp_MAM", "temp_JJA", "temp_SON",
                   "precip_DJF", "precip_MAM", "precip_JJA", "precip_SON")

pp_long <- bind_rows(
  climate_past_present_wide %>%
    transmute(plot_id, window = "past",
              across(all_of(paste0(climate_stems, "_past")),
                     .names = "{str_remove(.col, '_past$')}")),
  climate_past_present_wide %>%
    transmute(plot_id, window = "present",
              across(all_of(paste0(climate_stems, "_present")),
                     .names = "{str_remove(.col, '_present$')}"))
) %>%
  mutate(window = factor(window, levels = c("past", "present"))) %>%
  drop_na()

pp_matrix <- pp_long %>% select(all_of(climate_stems)) %>% as.matrix()
pca_pp <- prcomp(pp_matrix, center = TRUE, scale. = TRUE)

var_pp <- tibble(
  PC = paste0("PC", seq_along(pca_pp$sdev)),
  variance = pca_pp$sdev^2,
  pct_variance = 100 * variance / sum(variance)
) %>% mutate(pct_cumulative = cumsum(pct_variance))

loadings_pp <- as_tibble(pca_pp$rotation, rownames = "variable")
scores_pp <- as_tibble(pca_pp$x) %>%
  select(PC1, PC2) %>%
  bind_cols(pp_long %>% select(plot_id, window))

write_csv(var_pp, this_module_dir("outputs", "climate_pca_past_present_variance_explained.csv"))
write_csv(loadings_pp, this_module_dir("outputs", "climate_pca_past_present_loadings.csv"))

# Direction of the temporal shift in the displayed PC1-PC2 space.
centroids_pp <- scores_pp %>%
  group_by(window) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

past_cent <- centroids_pp %>% filter(window == "past")
present_cent <- centroids_pp %>% filter(window == "present")

shift_vec <- c(PC1 = present_cent$PC1 - past_cent$PC1,
               PC2 = present_cent$PC2 - past_cent$PC2)
shift_norm <- sqrt(sum(shift_vec^2))

# signed_projection ranks variables by how strongly their loading vector projects
# onto the past -> present centroid vector. Positive = toward present; negative =
# toward past. loading_length_PC12 indicates representation in the displayed plane.
temporal_alignment <- loadings_pp %>%
  transmute(
    variable,
    PC1_loading = PC1,
    PC2_loading = PC2,
    loading_length_PC12 = sqrt(PC1^2 + PC2^2),
    alignment_past_to_present =
      (PC1 * shift_vec["PC1"] + PC2 * shift_vec["PC2"]) /
      (sqrt(PC1^2 + PC2^2) * shift_norm),
    signed_projection =
      (PC1 * shift_vec["PC1"] + PC2 * shift_vec["PC2"]) / shift_norm
  ) %>%
  arrange(desc(abs(signed_projection)))

write_csv(temporal_alignment,
          this_module_dir("outputs", "climate_pca_past_present_alignment.csv"))

cat("\nVariables ranked by association with the past -> present displacement:\n")
print(temporal_alignment, n = Inf)

# Biplot with climate loadings and the past -> present centroid vector.
arrow_scale_pp <- 0.8 * max(abs(scores_pp$PC1), abs(scores_pp$PC2)) /
  max(abs(loadings_pp$PC1), abs(loadings_pp$PC2))

loadings_pp_scaled <- loadings_pp %>%
  mutate(PC1_arrow = PC1 * arrow_scale_pp,
         PC2_arrow = PC2 * arrow_scale_pp)

pc1_lab_pp <- paste0("PC1 (", round(var_pp$pct_variance[1], 1), "%)")
pc2_lab_pp <- paste0("PC2 (", round(var_pp$pct_variance[2], 1), "%)")

p_pca_shift <- ggplot(scores_pp, aes(x = PC1, y = PC2, color = window, fill = window)) +
  geom_point(alpha = 0.35, size = 1.5) +
  stat_ellipse(geom = "polygon", alpha = 0.12, linewidth = 0.6, level = 0.68) +
  geom_segment(data = loadings_pp_scaled,
               aes(x = 0, y = 0, xend = PC1_arrow, yend = PC2_arrow),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", linewidth = 0.4) +
  geom_text_repel(data = loadings_pp_scaled,
                  aes(x = PC1_arrow, y = PC2_arrow, label = variable),
                  inherit.aes = FALSE, size = 3, color = "black", max.overlaps = 20) +
  geom_segment(data = tibble(x = past_cent$PC1, y = past_cent$PC2,
                             xend = present_cent$PC1, yend = present_cent$PC2),
               aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 0.8, color = "black") +
  geom_point(data = centroids_pp, aes(x = PC1, y = PC2),
             inherit.aes = FALSE, shape = 21, size = 3,
             fill = "white", color = "black") +
  window_color_scale + window_fill_scale +
  labs(x = pc1_lab_pp, y = pc2_lab_pp) +
  climate_theme

ggsave(this_module_dir("figures", "climate_pca_past_present_shift.png"),
       p_pca_shift, width = 8, height = 7, dpi = 300, bg = "white")

cat("\nPast/present PCA and climate-variable alignment tables saved.\n")