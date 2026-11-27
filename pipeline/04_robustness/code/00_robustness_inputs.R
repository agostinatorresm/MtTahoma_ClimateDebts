# =============================================================================
# 00_robustness_inputs.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Builds the understory tables used by the sampling-reliability diagnostics
# (01a, 01b, 02, 03, 04) from the outputs of the cleaning module, so that the
# diagnostics use the same plots, species, present-day census and species
# complexes as the analyses (06_build_understory.R).
#
# Sampling design:
#   past    : one cover estimate per species for the whole plot (250, 500 or
#             1000 m2)
#   present : cover in four 1-m2 subplots (QN, QE, QS, QW) and presence in a
#             nested 50-m2 census (ALL)
#
# Tables (saved in one list):
#   subplot_grid : present-day subplot x species cover with true zeros: every
#                  species recorded in a plot (subplots or ALL) crossed with
#                  every subplot sampled in that plot, cover 0 where the
#                  species was not recorded in that subplot
#   all_census   : species recorded in the 50-m2 census, per plot
#   past         : past cover and occurrence per plot x species, with plot area
#   plot_info    : plot elevation and past plot area
#   denominator  : species used as denominator of relative cover (recorded in
#                  both surveys somewhere in the layer, as in the HMSC
#                  relative-abundance response)
#   species_set  : species evaluated by the diagnostics: recorded in more
#                  than five plots in each survey (frequency criterion of the
#                  frequent pool used in the HMSC models). All species are
#                  used to compute richness, relative cover and rarefaction;
#                  species-level results are reported for this set.
#
# Plots: those of understory_long.csv (both surveys, with cover subplots).
# Plots with two present-day censuses averaged in the analyses are left out
# of the subplot-level diagnostics, because their subplot records cannot be
# assigned to a census in understory_long_subplot.csv; they are listed in the
# console output.
#
# Inputs : pipeline/01_data_prep/01_cleaning/outputs/
#            understory_long.csv, understory_long_subplot.csv, plot_info.csv,
#            present_census_selected.csv, spp_info.csv
# Outputs: pipeline/04_robustness/tmp_files/robustness_inputs.rds
# =============================================================================

library(tidyverse)
library(here)

input_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
tmp_dir   <- here("pipeline", "04_robustness", "tmp_files")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

subplot_codes <- c("QN", "QE", "QS", "QW")
min_plots <- 6      # more than five plots in each survey, as the frequent (strict) pool of the HMSC models

understory_long <- read_csv(file.path(input_dir, "understory_long.csv"), show_col_types = FALSE)
long_subplot    <- read_csv(file.path(input_dir, "understory_long_subplot.csv"), show_col_types = FALSE)
plot_info_raw   <- read_csv(file.path(input_dir, "plot_info.csv"), show_col_types = FALSE)
census          <- read_csv(file.path(input_dir, "present_census_selected.csv"), show_col_types = FALSE)
spp_info        <- read_csv(file.path(input_dir, "spp_info.csv"), show_col_types = FALSE)

# -----------------------------------------------------------------------------
# 1. Plots
# -----------------------------------------------------------------------------
plots_analysis <- unique(understory_long$plot_id)
plots_two_census <- census %>%
  filter(stratum == "understory") %>%
  count(plot_id) %>%
  filter(n > 1) %>%
  pull(plot_id)
plots_subplot <- setdiff(plots_analysis, plots_two_census)

cat("Plots in understory_long.csv:", length(plots_analysis), "\n")
cat("Plots with two present-day censuses (left out of subplot-level diagnostics):",
    paste(intersect(plots_two_census, plots_analysis), collapse = ", "), "\n")

# Elevation (past survey) and past plot area, one row per plot
plot_info <- plot_info_raw %>%
  { if ("period" %in% names(.)) filter(., period == "past") else . } %>%
  distinct(plot_id, .keep_all = TRUE) %>%
  select(plot_id, elevation) %>%
  right_join(understory_long %>% filter(period == "past") %>% distinct(plot_id, plot_area),
             by = "plot_id") %>%
  filter(plot_id %in% plots_analysis)

# -----------------------------------------------------------------------------
# 2. Present-day subplot grid with true zeros, and 50-m2 census
# -----------------------------------------------------------------------------
present_sub <- long_subplot %>%
  filter(period == "present", plot_id %in% plots_subplot)

subplot_records <- present_sub %>%
  filter(replicate %in% subplot_codes) %>%
  group_by(plot_id, replicate, species) %>%
  summarise(cover = sum(cover, na.rm = TRUE), .groups = "drop")

all_census <- present_sub %>%
  filter(replicate == "ALL", occurrence == 1) %>%
  distinct(plot_id, species)

species_per_plot <- bind_rows(subplot_records %>% filter(cover > 0) %>% distinct(plot_id, species),
                              all_census) %>%
  distinct()
subplots_sampled <- present_sub %>%
  filter(replicate %in% subplot_codes) %>%
  distinct(plot_id, replicate)

subplot_grid <- species_per_plot %>%
  inner_join(subplots_sampled, by = "plot_id", relationship = "many-to-many") %>%
  left_join(subplot_records, by = c("plot_id", "replicate", "species")) %>%
  mutate(cover = coalesce(cover, 0), occurrence = as.integer(cover > 0))

cat("Present-day subplot grid:", nrow(subplot_grid), "rows,",
    n_distinct(subplot_grid$plot_id), "plots,", n_distinct(subplot_grid$species), "species\n")
cat("Subplots sampled per plot:\n")
print(subplots_sampled %>% count(plot_id, name = "n_subplots") %>% count(n_subplots, name = "n_plots"))

# -----------------------------------------------------------------------------
# 3. Past survey
# -----------------------------------------------------------------------------
past <- understory_long %>%
  filter(period == "past", occurrence == 1) %>%
  select(plot_id, species, cover, plot_area)

# -----------------------------------------------------------------------------
# 4. Denominator of relative cover
# -----------------------------------------------------------------------------
denominator <- understory_long %>%
  filter(occurrence == 1) %>%
  distinct(period, species) %>%
  count(species) %>%
  filter(n == 2) %>%
  pull(species)
cat("Species recorded in both surveys (relative-cover denominator):", length(denominator), "\n")

n_plots_survey <- understory_long %>%
  filter(occurrence == 1) %>%
  distinct(period, plot_id, species) %>%
  count(species, period) %>%
  pivot_wider(names_from = period, values_from = n, values_fill = 0)
species_set <- n_plots_survey %>% filter(past >= min_plots, present >= min_plots) %>% pull(species)
cat("Species in at least", min_plots, "plots in each survey (evaluated set):", length(species_set), "\n")

life_form <- spp_info %>% distinct(species, .keep_all = TRUE) %>% select(species, life_form)

saveRDS(list(subplot_grid = subplot_grid, all_census = all_census, past = past,
             plot_info = plot_info, denominator = denominator, species_set = species_set,
             life_form = life_form,
             plots_subplot = plots_subplot, understory_long = understory_long),
        file.path(tmp_dir, "robustness_inputs.rds"))
cat("Saved", file.path(tmp_dir, "robustness_inputs.rds"), "\n")
