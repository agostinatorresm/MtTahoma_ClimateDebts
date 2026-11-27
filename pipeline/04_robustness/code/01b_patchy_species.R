# =============================================================================
# 01b_patchy_species.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Species-level subplot occupancy in the present-day understory survey, for
# the evaluated species set (00: more than five plots in each survey). For
# each species, subplot occupancy in a plot is the proportion of sampled
# 1-m2 subplots in which it was recorded, in the plots where it was recorded
# in at least one subplot; mean subplot occupancy is its mean over those
# plots. Plots where the species was recorded only in the 50-m2 census are
# not included: there it has no cover value and does not enter the abundance
# models. Species with mean subplot occupancy below patchy_threshold are
# flagged as patchy: their present-day cover rests on few subplot records and
# is less reliably estimated than that of species recorded in most subplots.
# The flag is used for robustness analyses (models refitted without flagged
# species), not as a filter.
#
# The threshold (0.4) is close to the median of mean subplot occupancy across
# species; the console output gives the number of species flagged under
# other thresholds.
#
# Occupancy is computed for the evaluated set (frequent pool, 00), and the
# flag and figure are reported for the species modelled with HMSC (restricted
# pool: frequent species with complete traits), the set to which the patchy
# flag is applied in the sensitivity analyses (Supporting Information S3.4).
# If the restricted-pool list is not found, the evaluated set is used.
#
# Inputs : pipeline/04_robustness/tmp_files/robustness_inputs.rds (00)
#          pipeline/01_data_prep/02_species_pool/outputs/
#            understory_species_pool_strict_full_func_clim_cov.csv (restricted pool)
# Outputs: pipeline/04_robustness/outputs/accumulation_curve/accumulation_curve_summary.csv
#            (species, n_plots, mean_subplot_occupancy, patchy)
#          pipeline/04_robustness/figures/accumulation_curves/patchy_species.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

patchy_threshold <- 0.4

tmp_dir <- here("pipeline", "04_robustness", "tmp_files")
out_dir <- here("pipeline", "04_robustness", "outputs", "accumulation_curve")
fig_dir <- here("pipeline", "04_robustness", "figures", "accumulation_curves")
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
col_flag <- c("Patchy" = "#9A3F12", "Not patchy" = "grey55")

inp <- readRDS(file.path(tmp_dir, "robustness_inputs.rds"))

occupancy <- inp$subplot_grid %>%
  group_by(plot_id, species) %>%
  summarise(plot_occupancy = mean(occurrence), cover = mean(cover), .groups = "drop") %>%
  filter(plot_occupancy > 0) %>%
  group_by(species) %>%
  summarise(n_plots = n(),
            mean_subplot_occupancy = mean(plot_occupancy),
            mean_cover_present = mean(cover[cover > 0]),
            .groups = "drop") %>%
  filter(species %in% inp$species_set) %>%
  mutate(patchy = mean_subplot_occupancy < patchy_threshold) %>%
  left_join(inp$life_form, by = "species") %>%
  mutate(life_form = replace_na(life_form, "Unknown"),
         flag = if_else(patchy, "Patchy", "Not patchy"))

cat("Frequent pool -- species evaluated:", nrow(occupancy), "| patchy:", sum(occupancy$patchy), "\n")

# Restricted pool (species modelled with HMSC)
restricted_file <- here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                        "understory_species_pool_strict_full_func_clim_cov.csv")
if (file.exists(restricted_file)) {
  restricted <- read_csv(restricted_file, show_col_types = FALSE) %>%
    mutate(species = str_replace_all(species, "_", " ")) %>% pull(species)
  occupancy <- occupancy %>% filter(species %in% restricted)
  cat("Restricted pool -- species evaluated:", nrow(occupancy), "| patchy:", sum(occupancy$patchy), "\n")
} else {
  cat("Restricted-pool list not found; flag reported for the frequent pool\n")
}
print(summary(occupancy$mean_subplot_occupancy))
cat("\nSpecies flagged under other thresholds:\n")
print(tibble(threshold = c(0.25, 0.3, 0.4, 0.5, 0.6)) %>%
        mutate(n_patchy = map_int(threshold, ~ sum(occupancy$mean_subplot_occupancy < .x))))

write_csv(occupancy %>% select(species, n_plots, mean_subplot_occupancy, patchy) %>% arrange(mean_subplot_occupancy),
          file.path(out_dir, "accumulation_curve_summary.csv"))

# -----------------------------------------------------------------------------
# Figure: what the flag tracks
# -----------------------------------------------------------------------------
p_a <- ggplot(occupancy, aes(mean_subplot_occupancy, fill = flag)) +
  geom_histogram(binwidth = 0.05, boundary = 0, colour = "white") +
  geom_vline(xintercept = patchy_threshold, linetype = "dashed") +
  scale_fill_manual(values = col_flag, name = NULL) +
  labs(x = "Mean quadrat occupancy", y = "Number of species", title = "Distribution") +
  theme_mora()
p_b <- ggplot(occupancy, aes(n_plots, mean_subplot_occupancy, colour = flag)) +
  geom_hline(yintercept = patchy_threshold, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  scale_x_log10() +
  scale_colour_manual(values = col_flag, name = NULL) +
  labs(x = "Number of plots occupied (log scale)", y = "Mean quadrat occupancy", title = "Frequency") +
  theme_mora()
p_c <- ggplot(occupancy %>% filter(is.finite(mean_cover_present)),
              aes(mean_cover_present, mean_subplot_occupancy, colour = flag)) +
  geom_hline(yintercept = patchy_threshold, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.8) +
  scale_x_log10() +
  scale_colour_manual(values = col_flag, name = NULL) +
  labs(x = "Mean cover in occupied quadrats (%, log scale)", y = "Mean quadrat occupancy", title = "Abundance") +
  theme_mora()
p_d <- occupancy %>%
  count(life_form, flag) %>%
  group_by(life_form) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  ggplot(aes(pct, life_form, fill = flag)) +
  geom_col() +
  scale_fill_manual(values = col_flag, name = NULL) +
  labs(x = "% of species", y = NULL, title = "Life form") +
  theme_mora()

fig <- (p_a + p_b) / (p_c + p_d) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "patchy_species.png"), fig, width = 11, height = 9, dpi = 300, bg = "white")
cat("\nFigure written to", fig_dir, "\n")

# Panels for the combined diagnosis figures (08_diagnosis_figures.R)
panel_dir <- file.path(here("pipeline", "04_robustness", "tmp_files"), "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(occupancy = p_a), file.path(panel_dir, "panels_01b.rds"))
