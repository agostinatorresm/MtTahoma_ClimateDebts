# =============================================================================
# 05c_trait_effects_elevation_figure.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the local computer, after
#           24_trait_effects_elevation.R)
# Author: Agostina Torres
#
# Trait effects on the species' change between past and present along the
# elevational gradient (24_trait_effects_elevation.R): one panel per trait,
# overstory (top row) and understory (bottom row), occurrence and abundance in
# the two shades of orange of the species-change figures. Lines are posterior
# means; ribbons, 90% credible intervals. The posterior support for a change of
# the trait effect with elevation (Gamma for time x elevation) is given in each
# panel.
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/centered/trait_effects_elevation.csv (24)
#          pipeline/02_analysis/02_hmsc/outputs/centered/trait_elevation_interaction.csv (24)
# Outputs: pipeline/02_analysis/02_hmsc/figures/m1_final/ms1__fig3b_trait_effects_elevation_m1_final.png
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)
library(patchwork)

tag        <- "m1_final"
fig_prefix <- "ms1__"
col_resp   <- c(Occurrence = "#E07B39", Abundance = "#9A3F12")
strata_titles <- c(overstory = "\U0001F332 Overstory", understory = "\U0001F33F Understory")
trait_labels  <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity",
                   seed_mass = "Seed mass", wood_density = "Wood density",
                   SLA = "SLA", woodinesswoody = "Woody (vs non-woody)")

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))

theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = "black", linewidth = 0.6),
          axis.line.y = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text = element_text(colour = "black", size = 12),
          axis.title.x = element_text(colour = "black", size = 14, hjust = 0.5),
          axis.title.y = element_text(colour = "black", size = 14, hjust = 0.5),
          strip.text = element_text(colour = "black", face = "bold", size = 13),
          legend.text = element_text(colour = "black", size = 13),
          plot.title = element_text(colour = "black", size = 18, face = "bold"),
          plot.margin = margin(12, 10, 10, 10))
}

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
res_dir    <- file.path(module_dir, "outputs", "centered")
fig_dir    <- file.path(module_dir, "figures", tag)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

prep <- function(d) d %>%
  mutate(trait = factor(recode(trait, !!!trait_labels), levels = unique(trait_labels)),
         response = factor(str_to_sentence(response), levels = c("Occurrence", "Abundance")))
along <- read_csv(file.path(res_dir, "trait_effects_elevation.csv"), show_col_types = FALSE) %>% prep()
inter <- read_csv(file.path(res_dir, "trait_elevation_interaction.csv"), show_col_types = FALSE) %>% prep() %>%
  mutate(label = sprintf("%s: f = %.2f", response, support))

row_plot <- function(s) {
  d <- filter(along, stratum == s)
  lab <- filter(inter, stratum == s) %>% group_by(trait) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  ggplot(d, aes(elevation, mean, colour = response, fill = response)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.2, size = 3.6, lineheight = 0.95) +
    facet_wrap(~ trait, nrow = 1) +
    scale_colour_manual(values = col_resp, name = NULL) +
    scale_fill_manual(values = col_resp, name = NULL) +
    labs(x = "Elevation (m)", y = "Trait effect on change\nbetween past and present (Γ)",
         title = strata_titles[[s]]) +
    theme_mora() +
    theme(legend.position = "top")
}

fig <- (row_plot("overstory") / row_plot("understory")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
fig_file <- file.path(fig_dir, paste0(fig_prefix, "fig3b_trait_effects_elevation_", tag, ".png"))
ggsave(fig_file, fig,
       device = "png", dpi = 300, width = 12, height = 8, scale = 1.1, bg = "white")
cat("Figure written to", fig_dir, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
