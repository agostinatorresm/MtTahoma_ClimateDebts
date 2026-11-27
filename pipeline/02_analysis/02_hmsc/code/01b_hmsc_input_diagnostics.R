# =============================================================================
# 01b_hmsc_input_diagnostics.R
# Author: Agostina Torres
#
# Summary of the HMSC inputs built by 01_prepare_hmsc_inputs.R for one layer
# and species pool: sampling design, species occurrence and abundance, traits
# (raw, transformed and scaled), covariates, spatial design and phylogeny, and
# a screening of species prone to quasi-separation (occurrence per survey
# period and elevation tercile).
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/<stratum>_<pool>_hmsc_inputs.rds
# Outputs: pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_diag_*.csv
#          pipeline/02_analysis/02_hmsc/figures/<stratum>_<pool>_diag_*.png
# =============================================================================

library(here)
library(tidyverse)
library(ape)
library(hrbrthemes)
library(viridis)
library(ggpubr)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
stratum <- "overstory"   # "overstory" | "understory"
pool    <- "strict"

# -----------------------------------------------------------------------------
# Plot style (shared across project figures)
# -----------------------------------------------------------------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
}

theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid       = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x      = element_line(colour = "black", linewidth = 0.3),
          axis.line.y      = element_line(colour = "black", linewidth = 0.3),
          plot.margin      = margin(10, 10, 10, 10))
}

# -----------------------------------------------------------------------------
# Paths and inputs
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
run_tag    <- paste(stratum, pool, sep = "_")
out_dir    <- file.path(module_dir, "outputs")
fig_dir    <- file.path(module_dir, "figures")
walk(c(out_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

# File names always start with <stratum>_<pool>
diag_file <- function(name) file.path(out_dir, paste0(run_tag, "_diag_", name))
diag_fig  <- function(name) file.path(fig_dir, paste0(run_tag, "_diag_", name))

inp <- readRDS(file.path(module_dir, "data_clean", paste0(run_tag, "_hmsc_inputs.rds")))

if (is.null(inp$TrData_raw)) {
  stop("hmsc_inputs has no raw trait values (TrData_raw); re-run the current ",
       "version of 01_prepare_hmsc_inputs.R before this script.")
}

Y_pa  <- inp$Y_pa
Y_abu <- inp$Y_abu
XData <- inp$XData
past  <- XData$time == "past"

trait_names <- all.vars(inp$TrFormula)
num_traits  <- trait_names[map_lgl(trait_names, ~ is.numeric(inp$TrData[[.x]]))]

# -----------------------------------------------------------------------------
# 1. Overview
# -----------------------------------------------------------------------------
overview <- tibble(
  item = c("Stratum", "Species pool", "Plots", "Samples (plot x period)",
           "Species", "Traits", "Trait names",
           "Mean prevalence (all samples)", "Abundance observations (non-NA)",
           "Abundance variable", "Elevation range (m)",
           "Spatial extent (km, x / y)", "Minimum distance between plots (m)",
           "Phylogeny tips"),
  value = c(stratum, pool,
            nlevels(inp$studyDesign$plot), nrow(Y_pa),
            ncol(Y_pa), length(trait_names), paste(trait_names, collapse = ", "),
            sprintf("%.3f", mean(Y_pa)), sum(!is.na(Y_abu)),
            inp$abundance_var,
            paste(round(range(inp$plot_covariates$elevation)), collapse = " - "),
            paste(round(diff(range(inp$sData$x)), 1), "/", round(diff(range(inp$sData$y)), 1)),
            if (!is.null(inp$min_dist_km)) round(inp$min_dist_km * 1000, 1) else NA,
            length(inp$tree$tip.label))
)
overview

write_csv(overview, diag_file("overview.csv"))

# -----------------------------------------------------------------------------
# 2. Species: occurrence, abundance and traits (raw, transformed, scaled)
# -----------------------------------------------------------------------------
species_occ <- tibble(species = colnames(Y_pa)) %>%
  mutate(n_past          = colSums(Y_pa[past, , drop = FALSE]),
         n_present       = colSums(Y_pa[!past, , drop = FALSE]),
         prev_past       = round(n_past / sum(past), 3),
         prev_present    = round(n_present / sum(!past), 3),
         delta_prev      = prev_present - prev_past,
         n_abu           = colSums(!is.na(Y_abu)),
         mean_log10_abu  = round(colMeans(Y_abu, na.rm = TRUE), 3))

traits_raw <- inp$TrData_raw %>%
  rename_with(~ paste0(.x, "_raw"), all_of(trait_names))

traits_transformed <- inp$TrData_raw %>%
  mutate(across(any_of(inp$log_traits), log10)) %>%
  rename_with(~ ifelse(.x %in% inp$log_traits, paste0(.x, "_log10"), paste0(.x, "_untransformed")),
              all_of(trait_names)) %>%
  select(species, ends_with("_log10"))

traits_scaled <- inp$TrData %>%
  rownames_to_column("species") %>%
  rename_with(~ paste0(.x, "_scaled"), all_of(trait_names))

species_table <- species_occ %>%
  left_join(traits_raw, by = "species") %>%
  left_join(traits_transformed, by = "species") %>%
  left_join(traits_scaled, by = "species") %>%
  arrange(desc(n_past + n_present))
species_table

write_csv(species_table, diag_file("species_occurrence_traits.csv"))

# -----------------------------------------------------------------------------
# 3. Trait summary
# -----------------------------------------------------------------------------
trait_summary <- map_dfr(num_traits, function(tr) {
  raw    <- inp$TrData_raw[[tr]]
  scaled <- inp$TrData[[tr]]
  tibble(trait          = tr,
         transformation = if (tr %in% inp$log_traits) "log10, then scaled" else "scaled",
         raw_min        = signif(min(raw), 3),
         raw_median     = signif(median(raw), 3),
         raw_max        = signif(max(raw), 3),
         fold_range     = if (all(raw > 0)) round(max(raw) / min(raw), 1) else NA_real_,
         scaled_mean    = round(mean(scaled), 3),
         scaled_sd      = round(sd(scaled), 3),
         scaled_min     = round(min(scaled), 2),
         scaled_max     = round(max(scaled), 2),
         most_extreme   = rownames(inp$TrData)[which.max(abs(scaled))])
})
trait_summary

categorical_summary <- inp$TrData %>%
  select(where(is.factor)) %>%
  imap_dfr(~ tibble(trait = .y, level = names(table(.x)), n_species = as.integer(table(.x))))
categorical_summary

trait_cor <- cor(inp$TrData[, num_traits], method = "spearman")
round(trait_cor, 2)

write_csv(trait_summary, diag_file("trait_summary.csv"))
write_csv(as_tibble(trait_cor, rownames = "trait"), diag_file("trait_correlation_spearman.csv"))

# -----------------------------------------------------------------------------
# 4. Covariates
# -----------------------------------------------------------------------------
covariate_summary <- XData %>%
  select(where(is.numeric)) %>%
  imap_dfr(~ tibble(covariate = .y,
                    mean = round(mean(.x), 3), sd = round(sd(.x), 3),
                    min = round(min(.x), 2), max = round(max(.x), 2)))
covariate_summary

time_counts <- count(XData, time)
time_counts

write_csv(covariate_summary, diag_file("covariate_summary.csv"))

# -----------------------------------------------------------------------------
# 4b. Quasi-separation screening
# -----------------------------------------------------------------------------
# Occurrence of each species per survey period and elevation tercile (plot
# elevation). Cells where a species is never (or always) present indicate
# quasi-complete separation in the occurrence model: the corresponding
# coefficients (time, elevation, time x elevation) are weakly identified by the
# data and rely on the hierarchical prior. Species with few occurrences in
# total are flagged as well.
elev_breaks <- quantile(inp$plot_covariates$elevation, c(0, 1/3, 2/3, 1))

sample_info <- inp$XData_full %>%
  select(sample_id, plot_id, time) %>%
  left_join(select(inp$plot_covariates, plot_id, elevation), by = "plot_id") %>%
  mutate(elev_band = cut(elevation, elev_breaks, include.lowest = TRUE,
                         labels = c("low", "mid", "high")))
stopifnot(identical(sample_info$sample_id, rownames(Y_pa)))

occ_cells <- as_tibble(Y_pa) %>%
  mutate(time = sample_info$time, elev_band = sample_info$elev_band) %>%
  pivot_longer(-c(time, elev_band), names_to = "species", values_to = "presence") %>%
  group_by(species, time, elev_band) %>%
  summarise(n_samples  = n(),
            n_presence = sum(presence),
            prevalence = n_presence / n_samples,
            .groups = "drop")

separation_screen <- occ_cells %>%
  group_by(species) %>%
  summarise(n_presence_total = sum(n_presence),
            n_presence_past    = sum(n_presence[time == "past"]),
            n_presence_present = sum(n_presence[time == "present"]),
            n_cells_zero = sum(n_presence == 0),
            n_cells_full = sum(n_presence == n_samples),
            cells_zero   = paste(paste(time, elev_band, sep = ":")[n_presence == 0], collapse = ", "),
            .groups = "drop") %>%
  mutate(flag_period_absent = n_presence_past == 0 | n_presence_present == 0,
         flag_cells         = n_cells_zero > 0 | n_cells_full > 0,
         flag_rare          = n_presence_total < 10) %>%
  arrange(desc(flag_period_absent), desc(n_cells_zero + n_cells_full), n_presence_total)
separation_screen

count(separation_screen, flag_period_absent, flag_cells, flag_rare)

write_csv(occ_cells, diag_file("occurrence_period_elevation.csv"))
write_csv(separation_screen, diag_file("separation_screen.csv"))

# -----------------------------------------------------------------------------
# 5. Figures
# -----------------------------------------------------------------------------
sp_order <- species_table$species

# Prevalence in each period per species
pa1 <- species_table %>%
  select(species, prev_past, prev_present) %>%
  pivot_longer(-species, names_to = "period", values_to = "prevalence",
               names_prefix = "prev_") %>%
  mutate(species = factor(str_replace_all(species, "_", " "),
                          levels = rev(str_replace_all(sp_order, "_", " ")))) %>%
  ggplot(aes(x = prevalence, y = species)) +
  geom_line(aes(group = species), colour = "grey60") +
  geom_point(aes(colour = period), size = 3) +
  scale_colour_viridis_d(option = "viridis", end = 0.8, name = NULL) +
  labs(x = "Prevalence (proportion of plots)", y = NULL) +
  theme_mora() +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "top")
pa1

# Scaled traits per species
# Numeric traits only (categorical traits such as woodiness are summarised in
# categorical_summary)
pa2 <- traits_scaled %>%
  select(species, all_of(paste0(num_traits, "_scaled"))) %>%
  pivot_longer(-species, names_to = "trait", values_to = "value") %>%
  mutate(trait   = str_remove(trait, "_scaled$"),
         species = factor(str_replace_all(species, "_", " "),
                          levels = rev(str_replace_all(sp_order, "_", " ")))) %>%
  ggplot(aes(x = trait, y = species, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", value)), size = 3,
            colour = "white") +
  scale_fill_viridis_c(option = "viridis", name = "Scaled\nvalue") +
  labs(x = NULL, y = NULL) +
  theme_mora() +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))
pa2

fig_inputs <- ggarrange(pa1, pa2, ncol = 2, widths = c(1, 1), labels = c("A", "B"))
fig_inputs

ggsave(diag_fig("species_prevalence_traits.png"), fig_inputs,
       device = "png", dpi = 300, width = 10, height = 0.35 * length(sp_order) + 2,
       scale = 1.8, bg = "transparent", limitsize = FALSE)

# Prevalence per period and elevation tercile (empty cells outlined in red)
pa3 <- occ_cells %>%
  mutate(species = factor(str_replace_all(species, "_", " "),
                          levels = rev(str_replace_all(sp_order, "_", " "))),
         empty   = n_presence == 0) %>%
  ggplot(aes(x = elev_band, y = species, fill = prevalence)) +
  geom_tile(colour = "white") +
  geom_tile(data = ~ filter(.x, empty), fill = NA, colour = "red", linewidth = 0.5) +
  geom_text(aes(label = n_presence), size = 3, colour = "white") +
  facet_wrap(~ time, ncol = 2) +
  scale_fill_viridis_c(option = "viridis", name = "Prevalence") +
  labs(x = "Elevation tercile", y = NULL) +
  theme_mora() +
  theme(axis.text.y = element_text(face = "italic"))
pa3

ggsave(diag_fig("occurrence_period_elevation.png"), pa3,
       device = "png", dpi = 300, width = 7, height = 0.3 * length(sp_order) + 2,
       scale = 1.8, bg = "transparent", limitsize = FALSE)
