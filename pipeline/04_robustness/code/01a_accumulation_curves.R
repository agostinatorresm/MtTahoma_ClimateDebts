# =============================================================================
# 01a_accumulation_curves.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Sampling completeness of the present-day understory design, for the two
# responses used in the analyses:
#
#   Occurrence: species accumulation curve. Mean cumulative richness per plot
#     for k = 1-4 subplots (1-4 m2; averaged over every combination of k of
#     the sampled subplots) and for the 50-m2 census added to the four
#     subplots. A Michaelis-Menten (saturating) and a power-law
#     (non-saturating) species-area model are fitted to the mean curve and
#     compared by AIC.
#
#   Relative cover: stability of the relative-cover estimate with sampled
#     area. For each plot x species, relative cover (species cover / total
#     cover of the denominator species in the plot) is computed from every
#     combination of k subplots and compared with the estimate from all
#     subplots of the plot: mean absolute difference and correlation, per k.
#     A small difference at k = 3 indicates that the four-subplot estimate is
#     close to stable.
#
# Plots with exactly four sampled subplots are used for the curves.
#
# Inputs : pipeline/04_robustness/tmp_files/robustness_inputs.rds (00)
# Outputs: pipeline/04_robustness/tmp_files/accumulation_curve/
#            species_accumulation_curve.csv, species_accumulation_fit.csv,
#            relative_cover_stability.csv
#          pipeline/04_robustness/figures/accumulation_curves/
#            accumulation_occurrence_relative_cover.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

tmp_dir <- here("pipeline", "04_robustness", "tmp_files")
out_tmp <- file.path(tmp_dir, "accumulation_curve")
fig_dir <- here("pipeline", "04_robustness", "figures", "accumulation_curves")
walk(c(out_tmp, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(),
          axis.line  = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 11),
          axis.title.y = element_text(colour = "black", size = 11),
          legend.text  = element_text(colour = "black"),
          plot.title   = element_text(colour = "black", size = 11, face = "bold"),
          plot.margin  = margin(12, 10, 10, 10))
}
col_data <- "#9A3F12"; col_mm <- "#E07B39"; col_pow <- "grey45"

inp <- readRDS(file.path(tmp_dir, "robustness_inputs.rds"))
subplot_codes <- c("QN", "QE", "QS", "QW")

plots_4 <- inp$subplot_grid %>% distinct(plot_id, replicate) %>% count(plot_id) %>%
  filter(n == 4) %>% pull(plot_id)
grid <- inp$subplot_grid %>% filter(plot_id %in% plots_4)
cat("Plots with four subplots:", length(plots_4), "\n")

combos <- unlist(map(1:4, ~ combn(subplot_codes, .x, simplify = FALSE)), recursive = FALSE)

# -----------------------------------------------------------------------------
# 1. Species accumulation (occurrence)
# -----------------------------------------------------------------------------
occ_by_plot <- split(grid %>% filter(occurrence == 1), grid$plot_id[grid$occurrence == 1])
richness_k <- map_dfr(names(occ_by_plot), function(p) {
  d <- occ_by_plot[[p]]
  map_dfr(combos, function(cmb) tibble(plot_id = p, k = length(cmb),
                                       richness = n_distinct(d$species[d$replicate %in% cmb])))
}) %>%
  group_by(plot_id, k) %>% summarise(richness = mean(richness), .groups = "drop")

richness_all <- bind_rows(grid %>% filter(occurrence == 1) %>% distinct(plot_id, species),
                          inp$all_census %>% filter(plot_id %in% plots_4)) %>%
  distinct() %>% count(plot_id, name = "richness") %>% mutate(k = 5L)

acc <- bind_rows(richness_k, richness_all) %>%
  complete(plot_id = plots_4, k = 1:5, fill = list(richness = 0)) %>%
  group_by(k) %>%
  summarise(se = sd(richness) / sqrt(n()), richness = mean(richness), .groups = "drop") %>%
  mutate(area_m2 = c(1, 2, 3, 4, 50)[k],
         unit = c("1 subplot", "2 subplots", "3 subplots", "4 subplots", "+ 50-m2 census")[k])
write_csv(acc, file.path(out_tmp, "species_accumulation_curve.csv"))
print(acc)

mm  <- tryCatch(nls(richness ~ SSmicmen(area_m2, Smax, khalf), data = acc), error = function(e) NULL)
pow <- tryCatch(nls(richness ~ cc * area_m2^z, data = acc, start = list(cc = acc$richness[1], z = 0.3)),
                error = function(e) NULL)
fit_tab <- tibble(
  Smax = if (!is.null(mm)) coef(mm)[["Smax"]] else NA_real_,
  khalf = if (!is.null(mm)) coef(mm)[["khalf"]] else NA_real_,
  pct_Smax_4m2  = if (!is.null(mm)) 100 * 4 / (coef(mm)[["khalf"]] + 4) else NA_real_,
  pct_Smax_50m2 = if (!is.null(mm)) 100 * 50 / (coef(mm)[["khalf"]] + 50) else NA_real_,
  power_z = if (!is.null(pow)) coef(pow)[["z"]] else NA_real_,
  AIC_mm = if (!is.null(mm)) AIC(mm) else NA_real_,
  AIC_power = if (!is.null(pow)) AIC(pow) else NA_real_)
write_csv(fit_tab, file.path(out_tmp, "species_accumulation_fit.csv"))
cat("\nSpecies-area fits (five points: estimates are indicative):\n"); print(fit_tab)

area_seq <- exp(seq(log(0.8), log(60), length.out = 200))
fit_lines <- bind_rows(
  if (!is.null(mm))  tibble(area_m2 = area_seq, richness = predict(mm, data.frame(area_m2 = area_seq)), model = "Michaelis-Menten"),
  if (!is.null(pow)) tibble(area_m2 = area_seq, richness = predict(pow, data.frame(area_m2 = area_seq)), model = "Power law"))

p_occ <- ggplot(acc, aes(area_m2, richness)) +
  geom_line(data = fit_lines, aes(colour = model, linetype = model), linewidth = 0.8) +
  geom_errorbar(aes(ymin = richness - 1.96 * se, ymax = richness + 1.96 * se), width = 0, colour = col_data) +
  geom_point(size = 2.5, colour = col_data) +
  scale_x_log10(breaks = c(1, 2, 3, 4, 50)) +
  scale_colour_manual(values = c("Michaelis-Menten" = col_mm, "Power law" = col_pow), name = NULL) +
  scale_linetype_manual(values = c("Michaelis-Menten" = "dashed", "Power law" = "solid"), name = NULL) +
  labs(x = "Sampled area (m², log scale)", y = "Mean species richness per plot",
       title = "Occurrence: species accumulation") +
  theme_mora() + theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# 2. Stability of relative cover
# -----------------------------------------------------------------------------
grid_den <- grid %>% filter(species %in% inp$denominator)
cov_by_plot <- split(grid_den, grid_den$plot_id)
rel_k <- map_dfr(names(cov_by_plot), function(p) {
  d <- cov_by_plot[[p]]
  map_dfr(combos, function(cmb) {
    s <- d %>% filter(replicate %in% cmb) %>% group_by(species) %>%
      summarise(cover = mean(cover), .groups = "drop")
    tot <- sum(s$cover)
    if (tot == 0) return(NULL)
    s %>% transmute(plot_id = p, species, k = length(cmb), combo = paste(cmb, collapse = "+"),
                    rel = cover / tot)
  })
})
rel_full <- rel_k %>% filter(k == 4) %>% select(plot_id, species, rel_4 = rel)
stability <- rel_k %>%
  filter(k < 4) %>%
  inner_join(rel_full, by = c("plot_id", "species")) %>%
  filter(rel_4 > 0) %>%
  group_by(k) %>%
  summarise(mean_abs_diff = mean(abs(rel - rel_4)),
            median_abs_diff = median(abs(rel - rel_4)),
            mean_rel_4 = mean(rel_4),
            cor_log = cor(log10(rel + 1e-3), log10(rel_4 + 1e-3)),
            .groups = "drop") %>%
  bind_rows(tibble(k = 4L, mean_abs_diff = 0, median_abs_diff = 0, cor_log = 1))
write_csv(stability, file.path(out_tmp, "relative_cover_stability.csv"))
cat("\nRelative cover: difference from the four-subplot estimate (species present in the plot):\n")
print(stability)

p_rel <- ggplot(stability, aes(k, mean_abs_diff)) +
  geom_line(colour = col_data, linewidth = 0.8) +
  geom_point(size = 2.5, colour = col_data) +
  scale_x_continuous(breaks = 1:4) +
  labs(x = "Quadrats used (k, 1 m² each)", y = "Mean absolute difference in relative cover\nfrom the four-quadrat estimate",
       title = "Relative cover: estimate stability") +
  theme_mora()

fig <- p_occ + p_rel + plot_annotation(tag_levels = "A")
ggsave(file.path(fig_dir, "accumulation_occurrence_relative_cover.png"), fig,
       width = 11, height = 5, dpi = 300, bg = "white")
cat("\nFigure written to", fig_dir, "\n")

# Panels for the combined diagnosis figures (08_diagnosis_figures.R)
panel_dir <- file.path(here("pipeline", "04_robustness", "tmp_files"), "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(accumulation = p_occ, cover_stability = p_rel), file.path(panel_dir, "panels_01a.rds"))
