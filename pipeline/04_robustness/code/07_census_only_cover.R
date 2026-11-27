# =============================================================================
# 07_census_only_cover.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Understory cover change with abundance paired by plot.
#
# Detection differs between the past (whole plot, 250-1000 m2) and the present
# (four 1-m2 subplots for cover, plus a 50-m2 census for presence). Species
# recorded only in the 50-m2 census are present without a cover value, and
# local losses and gains partly reflect detection in a smaller area. Comparing
# cover only where a species has a cover value in both past and present of the
# same plot removes the effect of detection on the abundance data; detection
# differences remain in the occurrence data, which are treated by reducing
# the past data to the present sampling effort (05-06).
#
# Treatments compared, per species (understory HMSC species):
#   Relative cover in both: species cover / total cover of the modelled
#   species with a cover value in that plot and period (row total of the
#   cover matrix of the modelled species); missing values applied afterwards.
#   current : every presence with cover; census-only pairs missing in the
#             past as well (check responses abuall/relall)
#   paired  : only species with cover in both past and present of the same
#             plot (main responses abu/rel)
# Reported: number of cover values per species and species below the minimum
# for the abundance models (min_abu_obs), and the change in mean log10
# relative cover conditional on presence (present - past).
#
# Inputs : pipeline/04_robustness/tmp_files/robustness_inputs.rds (00)
#          pipeline/01_data_prep/02_species_pool/outputs/understory_species_pool_strict_full_func_clim_cov.csv
# Outputs: pipeline/04_robustness/outputs/census_only_cover/census_only_species.csv
#          pipeline/04_robustness/figures/figSR7_census_only_cover.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

tmp_dir <- here("pipeline", "04_robustness", "tmp_files")
out_dir <- here("pipeline", "04_robustness", "outputs", "census_only_cover")
fig_dir <- here("pipeline", "04_robustness", "figures", "census_only_cover")
sup_dir <- here("pipeline", "04_robustness", "figures")   # supplementary and diagnostic figures
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)
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
col_var <- c(current = "grey55", paired = "#9A3F12")
lab_var <- c(current = "All presences with cover (check)",
             paired = "Paired by plot (main)")
min_abu_obs <- 5   # as in 01_prepare_hmsc_inputs.R

inp <- readRDS(file.path(tmp_dir, "robustness_inputs.rds"))
hmsc_species <- read_csv(here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                              "understory_species_pool_strict_full_func_clim_cov.csv"),
                         show_col_types = FALSE) %>% pull(species) %>% unique()

# -----------------------------------------------------------------------------
# 1. Census-only pairs and relative cover
# -----------------------------------------------------------------------------
occ <- inp$understory_long %>%
  filter(occurrence == 1) %>%
  distinct(period, plot_id, species, .keep_all = TRUE) %>%
  select(period, plot_id, species, cover) %>%
  mutate(has_cover = !is.na(cover) & cover > 0)

census_only <- occ %>% filter(period == "present", !has_cover) %>% distinct(plot_id, species)
cat("Present occurrences (species x plot):", sum(occ$period == "present"),
    "| census-only:", nrow(census_only),
    sprintf("(%.0f%%)", 100 * nrow(census_only) / sum(occ$period == "present")), "\n")

# Relative cover on the cover matrix of the modelled species: species cover /
# total cover of the modelled species with a cover value in that plot and
# period (as in 01_prepare_hmsc_inputs.R). Missing values are applied after
# the totals.
rel_base <- occ %>%
  filter(species %in% hmsc_species) %>%
  group_by(period, plot_id) %>%
  mutate(total = sum(if_else(has_cover, cover, 0)),
         rel = if_else(has_cover & total > 0, cover / total, NA_real_)) %>%
  ungroup()

# All presences with cover; census-only pairs missing in the past as well
rel_current <- rel_base %>%
  left_join(census_only %>% mutate(census_only = TRUE), by = c("plot_id", "species")) %>%
  mutate(census_only = coalesce(census_only, FALSE),
         rel_current = if_else(period == "past" & census_only, NA_real_, rel)) %>%
  select(period, plot_id, species, census_only, rel_current)

# Paired: only species with cover in both times of the plot
paired_pairs <- rel_base %>%
  filter(has_cover) %>%
  distinct(period, plot_id, species) %>%
  count(plot_id, species) %>%
  filter(n == 2) %>%
  select(plot_id, species) %>%
  mutate(paired = TRUE)
rel_paired <- rel_base %>%
  left_join(paired_pairs, by = c("plot_id", "species")) %>%
  mutate(rel_paired = if_else(coalesce(paired, FALSE), rel, NA_real_)) %>%
  select(period, plot_id, species, rel_paired)

rel <- rel_current %>% left_join(rel_paired, by = c("period", "plot_id", "species"))

# -----------------------------------------------------------------------------
# 2. Per species
# -----------------------------------------------------------------------------
sp <- rel %>% filter(species %in% hmsc_species)

counts <- sp %>%
  group_by(species) %>%
  summarise(n_present = sum(period == "present"),
            n_census_only = sum(period == "present" & census_only),
            prop_census_only = n_census_only / n_present,
            n_values_current = sum(!is.na(rel_current)),
            n_values_paired  = sum(!is.na(rel_paired)),
            n_plots_paired   = n_values_paired / 2,
            .groups = "drop")

change <- sp %>%
  pivot_longer(c(rel_current, rel_paired), names_to = "variant", values_to = "r", names_prefix = "rel_") %>%
  filter(!is.na(r), r > 0) %>%
  group_by(variant, species, period) %>%
  summarise(m = mean(log10(r)), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = m) %>%
  mutate(change = present - past)

species_out <- change %>%
  select(variant, species, change) %>%
  pivot_wider(names_from = variant, values_from = change, names_prefix = "change_") %>%
  right_join(counts, by = "species") %>%
  left_join(inp$life_form, by = "species") %>%
  arrange(n_plots_paired)
write_csv(species_out, file.path(out_dir, "census_only_species.csv"))

cat("\nCover values of the HMSC species (both times): current", sum(counts$n_values_current),
    "| paired", sum(counts$n_values_paired), "\n")
cat("Plots per species with cover in both times: min", min(counts$n_plots_paired),
    "| median", median(counts$n_plots_paired), "\n")
cat("Species below", min_abu_obs, "cover values (excluded from the abundance models): current",
    sum(counts$n_values_current < min_abu_obs), "| paired", sum(counts$n_values_paired < min_abu_obs), "\n")
cat("Species with fewer than 10 paired plots:",
    paste(counts$species[counts$n_plots_paired < 10], collapse = ", "), "\n")
cat("\nChange in mean log10 relative cover conditional on presence (present - past):\n")
print(change %>% group_by(variant) %>%
        summarise(mean_change = mean(change, na.rm = TRUE),
                  n_increase = sum(change > 0, na.rm = TRUE),
                  n_decrease = sum(change < 0, na.rm = TRUE), .groups = "drop"))
cat("Correlation of species changes, current vs paired:",
    round(cor(species_out$change_current, species_out$change_paired, use = "complete.obs"), 3), "\n")

# -----------------------------------------------------------------------------
# Figure: (A) change per species; (B) shift against the proportion of
# census-only occurrences
# -----------------------------------------------------------------------------
ord <- species_out %>% arrange(change_paired) %>% pull(species)
pa <- change %>%
  mutate(species = factor(species, levels = ord), variant = factor(variant, levels = names(lab_var))) %>%
  ggplot(aes(change, species, colour = variant)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.8) +
  scale_colour_manual(values = col_var, labels = lab_var, name = NULL) +
  labs(x = "Change in mean log10 relative cover\n(present - past, conditional on presence)", y = NULL,
       title = "Species change") +
  theme_mora() +
  theme(axis.text.y = element_text(face = "italic", colour = "black", size = 8))

pb <- species_out %>%
  mutate(shift = change_paired - change_current) %>%
  ggplot(aes(n_plots_paired, shift)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.8, alpha = 0.85, colour = col_var[["paired"]]) +
  labs(x = "Plots with cover in both times (paired)",
       y = "Shift in change (paired - current)", title = "Effect of pairing") +
  theme_mora()

fig <- (pa | pb) + plot_layout(widths = c(1.1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
ggsave(file.path(sup_dir, "figSR7_census_only_cover.png"), fig,
       width = 12, height = 2 + 0.17 * length(ord), dpi = 300, bg = "white", limitsize = FALSE)
cat("\nFigure written to", sup_dir, "\n")
