# =============================================================================
# 09_species_affinity_figure.R
#
# Pipeline: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/
# Author: Agostina Torres
#
# Thermal (MAT) and moisture (MAP) affinities of the overstory species of the
# shared pool (all used in the climatic-affinity shift analyses) and of the
# understory species of the restricted pool (those modelled in HMSC; the
# understory shared pool has too many species to show legibly).
#
# For each species, the point is the affinity used as a single value (mean of
# the occurrence values truncated to the 5th-95th percentile) and the bar is
# the 5th-95th percentile of the 100 resampled values used to propagate
# affinity uncertainty into CTI and CMI, i.e. the spread of the climate across
# the species' occurrences. Each panel is ordered by its own affinity. Points
# are coloured with the thermal and moisture scales used in the other figures
# (Figure 3, Figure S4.5). Overstory species of the shared pool that are not
# in the restricted pool (not modelled in HMSC) are drawn as triangles.
#
# Inputs : data_clean/species_climatic_affinity.csv
#          data_clean/species_climatic_affinity_resamples.csv
#          02_species_pool/outputs/<layer>_species_pool_broad_full_clim_cov.csv (shared pool)
#          02_analysis/02_hmsc/outputs/<layer>_strict_diag_species_occurrence_traits.csv
#            (restricted pool)
# Output : figures/species_climatic_affinity.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

affinity_dir     <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
species_pool_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)
hmsc_dir         <- function(...) here("pipeline", "02_analysis", "02_hmsc", ...)
dir.create(affinity_dir("figures"), recursive = TRUE, showWarnings = FALSE)

# Palettes and limits shared with Figure 3 (02_hmsc/05b_species_change_centered_figure.R)
COL_MAT    <- viridisLite::plasma(256, end = 0.80)
MAT_LIMITS <- c(0, 12.5)
COL_MAP    <- c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D")
MAP_LIMITS <- c(800, 2700)

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.x = element_line(colour = "black", linewidth = 0.4),
          axis.text = element_text(colour = "black", size = 10),
          axis.title.x = element_text(colour = "black", size = 11, hjust = 0.5),
          axis.title.y = element_blank(),
          plot.title = element_text(colour = "black", size = 13, face = "bold"),
          plot.margin = margin(4, 8, 4, 4))
}

# -----------------------------------------------------------------------------
# 1. Affinities, their spread and the species pools
# -----------------------------------------------------------------------------
aff <- read_csv(affinity_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean)
spread <- read_csv(affinity_dir("data_clean", "species_climatic_affinity_resamples.csv"), show_col_types = FALSE) %>%
  group_by(species) %>%
  summarise(MAT_lo = quantile(MAT_resample, 0.05), MAT_hi = quantile(MAT_resample, 0.95),
            MAP_lo = quantile(MAP_resample, 0.05), MAP_hi = quantile(MAP_resample, 0.95), .groups = "drop")

layers <- c("overstory", "understory")
shared <- map_dfr(layers, function(l) {
  read_csv(species_pool_dir("outputs", paste0(l, "_species_pool_broad_full_clim_cov.csv")), show_col_types = FALSE) %>%
    distinct(species) %>% mutate(layer = l)
})
restricted <- map_dfr(layers, function(l) {
  f <- hmsc_dir("outputs", paste0(l, "_strict_diag_species_occurrence_traits.csv"))
  if (!file.exists(f)) return(tibble(species = character(), layer = character()))
  d <- read_csv(f, show_col_types = FALSE)
  scol <- intersect(c("species", "Species", "sp"), names(d))[1]
  tibble(species = str_replace_all(d[[scol]], "_", " "), layer = l)
}) %>% distinct()
if (!any(restricted$layer == "understory")) stop("Understory restricted-pool species not found in 02_hmsc/outputs")

sp <- shared %>%
  left_join(aff, by = "species") %>%
  left_join(spread, by = "species") %>%
  mutate(restricted = if (nrow(restricted)) paste(species, layer) %in% paste(restricted$species, restricted$layer) else TRUE)
missing_aff <- filter(sp, is.na(MAT_mean) | is.na(MAP_mean))
if (nrow(missing_aff)) cat("Shared-pool species without affinities (not shown):", paste(missing_aff$species, collapse = ", "), "\n")
sp <- filter(sp, !is.na(MAT_mean), !is.na(MAP_mean)) %>%
  filter(layer == "overstory" | restricted)      # understory: restricted pool only

cat("\nSpecies shown (overstory: shared pool; understory: restricted pool):\n")
print(count(sp, layer, restricted))

# -----------------------------------------------------------------------------
# 2. Figure
# -----------------------------------------------------------------------------
affinity_panel <- function(l, var, cols, limits, xlab, legend_title, title = NULL, name_size = 9) {
  m <- paste0(var, "_mean"); lo <- paste0(var, "_lo"); hi <- paste0(var, "_hi")
  d <- filter(sp, layer == l) %>% mutate(species = fct_reorder(species, .data[[m]]))
  ggplot(d, aes(y = species)) +
    geom_linerange(aes(xmin = .data[[lo]], xmax = .data[[hi]]), colour = "grey60", linewidth = 0.4) +
    geom_point(aes(x = .data[[m]], fill = .data[[m]], shape = restricted), colour = "black", size = 2.1, stroke = 0.3) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), guide = "none") +
    scale_fill_gradientn(colours = cols, limits = limits, oob = scales::squish, name = legend_title,
                         guide = guide_colourbar(barheight = unit(3.5, "cm"), barwidth = unit(0.35, "cm"))) +
    labs(x = xlab, title = title) +
    theme_mora() +
    theme(axis.text.y = element_text(face = "italic", size = name_size), axis.ticks.y = element_blank(),
          legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 9))
}
xlab_mat <- "Thermal affinity (MAT, \u00b0C)"; xlab_map <- "Moisture affinity (MAP, mm)"
pA <- affinity_panel("overstory",  "MAT", COL_MAT, MAT_LIMITS, xlab_mat, "MAT (\u00b0C)", "A. Overstory", 9)
pB <- affinity_panel("overstory",  "MAP", COL_MAP, MAP_LIMITS, xlab_map, "MAP (mm)",      "B. Overstory", 9)
pC <- affinity_panel("understory", "MAT", COL_MAT, MAT_LIMITS, xlab_mat, "MAT (\u00b0C)", "C. Understory", 7.5)
pD <- affinity_panel("understory", "MAP", COL_MAP, MAP_LIMITS, xlab_map, "MAP (mm)",      "D. Understory", 7.5)
n_over  <- sum(sp$layer == "overstory")
n_under <- sum(sp$layer == "understory")
fig <- ((pA | pB) / (pC | pD)) + plot_layout(heights = c(n_over + 4, n_under + 4))

fig_file <- affinity_dir("figures", "species_climatic_affinity.png")
ggsave(fig_file, fig, width = 11, height = 2.5 + 0.17 * (n_over + n_under), dpi = 300, bg = "white",
       limitsize = FALSE)
cat("\nFigure written to", fig_file, "\n")

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)