# =============================================================================
# 02_bootstrapping.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Is the change in understory cover between surveys distinguishable from the
# sampling variation of the present-day subplot design? Applied to the
# abundance responses of the analyses (cover and relative cover when
# present), for plot x species pairs recorded in both surveys with a
# present-day cover value.
#
# Following the null-model logic of Zhou et al. (2021, Ecological
# Applications 31: e02437): for each plot, the sampled subplots are
# resampled with replacement (same number as sampled), jointly for all
# species, n_boot times. Each replicate gives present-day cover per species
# (mean over the resampled subplots) and relative cover (species cover /
# total cover of the denominator species). The past value is compared with
# the bootstrap distribution: two-sided p-value = proportion of replicates
# whose distance from the bootstrap mean is at least the distance of the past
# value. P-values are adjusted for multiple tests (Benjamini-Hochberg).
# Resampling subplots jointly keeps the dependence among species within a
# subplot, which relative cover requires.
#
# Per species (evaluated set, 00): mean change over its plots, with all changes counted (raw) or
# only FDR-significant changes counted (others set to 0); shrinkage = raw -
# FDR. Species-level reliability flag (used in 04): at least half of its plot
# tests valid and relative shrinkage (|shrinkage| / |raw change|) below 0.75,
# for relative cover.
#
# Inputs : pipeline/04_robustness/tmp_files/robustness_inputs.rds (00)
# Outputs: pipeline/04_robustness/outputs/bootstrapping/
#            bootstrap_plot_species.csv  (one row per plot x species x response)
#            bootstrap_species.csv       (one row per species x response)
#          pipeline/04_robustness/figures/bootstrapping/bootstrap_relative_cover.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

n_boot <- 2000
alpha  <- 0.05
set.seed(20260916)

tmp_dir <- here("pipeline", "04_robustness", "tmp_files")
out_dir <- here("pipeline", "04_robustness", "outputs", "bootstrapping")
fig_dir <- here("pipeline", "04_robustness", "figures", "bootstrapping")
walk(c(out_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

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
col_raw <- "#E07B39"; col_fdr <- "#9A3F12"

inp <- readRDS(file.path(tmp_dir, "robustness_inputs.rds"))

# Past cover and relative cover per plot x species
past <- inp$past %>%
  group_by(plot_id) %>%
  mutate(past_rel = cover / sum(cover[species %in% inp$denominator], na.rm = TRUE)) %>%
  ungroup() %>%
  transmute(plot_id, species, past_cover = cover, past_rel)

# -----------------------------------------------------------------------------
# 1. Bootstrap per plot
# -----------------------------------------------------------------------------
boot_plot <- function(d) {
  M <- d %>% select(replicate, species, cover) %>%
    pivot_wider(names_from = species, values_from = cover, values_fill = 0)
  M <- as.matrix(M[, -1, drop = FALSE])
  n <- nrow(M)
  if (n < 2) return(NULL)
  den <- colnames(M) %in% inp$denominator
  # resampling counts: n_boot x n subplots
  counts <- t(rmultinom(n_boot, n, rep(1 / n, n)))
  cov_b  <- counts %*% M / n                                   # n_boot x species
  tot_b  <- rowSums(cov_b[, den, drop = FALSE])
  rel_b  <- cov_b / ifelse(tot_b > 0, tot_b, NA)
  obs_cov <- colMeans(M)
  obs_rel <- obs_cov / sum(obs_cov[den])
  tibble(species = colnames(M), n_subplots = n,
         present_cover = obs_cov, present_rel = obs_rel,
         cov_draws = asplit(cov_b, 2), rel_draws = asplit(rel_b, 2))
}

grid <- inp$subplot_grid
plots <- split(grid, grid$plot_id)
cat("Bootstrapping", length(plots), "plots (", n_boot, "replicates each)...\n")
boot <- imap_dfr(plots, function(d, p) {
  b <- boot_plot(d)                       # NULL for plots with a single subplot
  if (is.null(b)) NULL else mutate(b, plot_id = p, .before = 1)
})
cat("Plots with a single subplot (not testable):", sum(map_int(plots, ~ n_distinct(.x$replicate)) < 2), "\n")

test_one <- function(draws, past_value) {
  draws <- draws[is.finite(draws)]
  if (length(draws) < 10 || is.na(past_value)) return(NA_real_)
  m <- mean(draws)
  mean(abs(draws - m) >= abs(past_value - m))
}

tests <- boot %>%
  inner_join(past, by = c("plot_id", "species")) %>%
  filter(present_cover > 0) %>%
  mutate(p_cover = map2_dbl(cov_draws, past_cover, test_one),
         p_rel   = map2_dbl(rel_draws, past_rel, test_one)) %>%
  select(-cov_draws, -rel_draws)

results <- bind_rows(
  tests %>% transmute(plot_id, species, response = "cover", n_subplots,
                      past = past_cover, present = present_cover, p_value = p_cover),
  tests %>% transmute(plot_id, species, response = "relative_cover", n_subplots,
                      past = past_rel, present = present_rel, p_value = p_rel)) %>%
  group_by(response) %>%
  mutate(change = present - past,
         p_adjusted = p.adjust(p_value, method = "BH"),
         significant_fdr = !is.na(p_adjusted) & p_adjusted < alpha,
         direction = if_else(change > 0, "increase", "decrease")) %>%
  ungroup() %>%
  left_join(inp$plot_info %>% select(plot_id, elevation), by = "plot_id")

write_csv(results, file.path(out_dir, "bootstrap_plot_species.csv"))
cat("\nPlot x species tests (present in both surveys, with cover):\n")
print(results %>% count(response, significant_fdr, direction))

# -----------------------------------------------------------------------------
# 2. Species summaries and reliability flag
# -----------------------------------------------------------------------------
n_pairs <- inp$understory_long %>%
  filter(occurrence == 1) %>%
  distinct(period, plot_id, species) %>%
  count(plot_id, species) %>% filter(n == 2) %>%
  filter(plot_id %in% inp$plots_subplot) %>%
  count(species, name = "n_plots_both")

species_summary <- results %>%
  filter(species %in% inp$species_set) %>%
  mutate(change_fdr = if_else(significant_fdr, change, 0)) %>%
  group_by(response, species) %>%
  summarise(n_tests = sum(!is.na(p_value)),
            mean_elevation = mean(elevation, na.rm = TRUE),
            mean_change_raw = mean(change[!is.na(p_value)]),
            mean_change_fdr = mean(change_fdr[!is.na(p_value)]),
            .groups = "drop") %>%
  left_join(n_pairs, by = "species") %>%
  mutate(pct_valid = n_tests / n_plots_both,
         shrinkage = mean_change_raw - mean_change_fdr,
         rel_shrinkage = if_else(abs(mean_change_raw) > 1e-9, abs(shrinkage) / abs(mean_change_raw), NA_real_),
         bootstrap_reliable = pct_valid >= 0.5 & (is.na(rel_shrinkage) | rel_shrinkage < 0.75)) %>%
  left_join(inp$life_form, by = "species")

write_csv(species_summary, file.path(out_dir, "bootstrap_species.csv"))
cat("\nSpecies reliability (relative cover):\n")
print(species_summary %>% filter(response == "relative_cover") %>% count(bootstrap_reliable))

# -----------------------------------------------------------------------------
# 3. Figure (relative cover): tests and raw vs FDR change per species
# -----------------------------------------------------------------------------
rel <- results %>% filter(response == "relative_cover", !is.na(p_value))
p_a <- rel %>%
  mutate(result = case_when(significant_fdr & change > 0 ~ "Increase (FDR)",
                            significant_fdr & change < 0 ~ "Decrease (FDR)",
                            TRUE ~ "Not significant"),
         neg_log_p = -log10(pmax(p_value, 1 / n_boot))) %>%
  ggplot(aes(change, neg_log_p, colour = result)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_colour_manual(values = c("Increase (FDR)" = col_raw, "Decrease (FDR)" = col_fdr,
                                 "Not significant" = "grey70"), name = NULL) +
  labs(x = "Change in relative cover (present - past)", y = "-log10 p-value",
       title = "Plot x species tests") +
  theme_mora() + theme(legend.position = "bottom")

sp <- species_summary %>% filter(response == "relative_cover", n_tests > 0)
p_b <- ggplot(sp) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = mean_elevation, xend = mean_elevation, y = mean_change_raw, yend = mean_change_fdr),
               colour = "grey60", linewidth = 0.3) +
  geom_point(aes(mean_elevation, mean_change_raw, shape = "All changes"), colour = col_raw, size = 2) +
  geom_point(aes(mean_elevation, mean_change_fdr, shape = "FDR-significant only"), colour = col_fdr, size = 2) +
  scale_shape_manual(values = c("All changes" = 16, "FDR-significant only" = 17), name = NULL) +
  labs(x = "Mean elevation of occurrence (m)", y = "Mean change in relative cover",
       title = "Species") +
  theme_mora() + theme(legend.position = "bottom")

fig <- p_a + p_b + plot_annotation(tag_levels = "A")
ggsave(file.path(fig_dir, "bootstrap_relative_cover.png"), fig, width = 11, height = 5, dpi = 300, bg = "white")
cat("\nFigure written to", fig_dir, "\n")

# Panels for the combined diagnosis figures (08_diagnosis_figures.R)
panel_dir <- file.path(here("pipeline", "04_robustness", "tmp_files"), "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(bootstrap = p_a), file.path(panel_dir, "panels_02.rds"))
