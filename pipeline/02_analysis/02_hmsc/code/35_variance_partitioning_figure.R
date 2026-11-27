# =============================================================================
# 35_variance_partitioning_figure.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (after 31_phylogenetic_signal_hmsc.R
#           and 33_variance_partitioning.R; runs from their csv outputs, no fits needed)
# Author: Agostina Torres
#
# Figure S6.1: what the fixed effects, the plot-level random effect and the
# phylogeny explain in the final HMSC models.
#   (A) Mean proportion of the explained variance attributed to time (time and
#       time x elevation), elevation and the plot-level random effect, averaged
#       across species, for the occurrence and abundance models of each layer.
#   (B) Phylogenetic signal in species responses (rho; posterior mean and 90%
#       credible interval).
# In both panels the occurrence model is shown to the left of the abundance
# model within each layer, so that the two panels can be read in the same order.
#
# Inputs : outputs/variance_partitioning/variance_partitioning_summary.csv (33)
#          outputs/phylo_signal/hmsc_rho.csv (31)
# Output : figures/m1_final/ms1__figS6_1_variance_partitioning_m1_final.png
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)
library(patchwork)

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
tag        <- "m1_final"
fig_prefix <- "ms1__"
fig_dir    <- file.path(module_dir, "figures", tag)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

vp_file  <- file.path(module_dir, "outputs", "variance_partitioning", "variance_partitioning_summary.csv")
rho_file <- file.path(module_dir, "outputs", "phylo_signal", "hmsc_rho.csv")
for (f in c(vp_file, rho_file)) if (!file.exists(f)) stop("Run 31 and 33 first (or copy the file): ", f)

layer_labels    <- c(overstory = "Overstory", understory = "Understory")
response_labels <- c(occurrence = "Occurrence", abundance = "Abundance")   # occurrence first
component_labels <- c(Time = "Time (incl. time \u00d7 elevation)", Elevation = "Elevation",
                      Random = "Plot-level random effect")

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
# Mean proportion per component; the random-effect column is named after the
# random level by Hmsc (e.g. "mean_Random: plot"), so it is matched by prefix
vp <- read_csv(vp_file, show_col_types = FALSE) %>%
  select(stratum, response, starts_with("mean_")) %>%
  pivot_longer(starts_with("mean_"), names_to = "component", values_to = "prop") %>%
  mutate(component = str_remove(component, "^mean_"),
         component = case_when(str_detect(component, "^Random") ~ "Random", TRUE ~ component)) %>%
  filter(component %in% names(component_labels)) %>%
  mutate(stratum   = factor(layer_labels[stratum], levels = layer_labels),
         response  = factor(response_labels[response], levels = response_labels),
         component = factor(component_labels[component], levels = rev(component_labels)))

rho <- read_csv(rho_file, show_col_types = FALSE) %>%
  mutate(stratum  = factor(layer_labels[stratum], levels = layer_labels),
         response = factor(response_labels[response], levels = response_labels))

# -----------------------------------------------------------------------------
# 2. Figure
# -----------------------------------------------------------------------------
theme_s6 <- theme_ipsum(base_size = 12) +
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.text = element_text(colour = "black", size = 11),
        axis.title.y = element_text(colour = "black", size = 12),
        strip.text = element_text(colour = "black", face = "bold", size = 12, hjust = 0.5),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "bottom", legend.text = element_text(size = 11))

p_vp <- ggplot(vp, aes(x = response, y = prop, fill = component)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  facet_wrap(~ stratum, nrow = 1) +
  scale_fill_viridis_d(option = "mako", begin = 0.25, end = 0.85, name = NULL,
                       guide = guide_legend(reverse = TRUE, nrow = 1)) +
  # proportions are rounded in the csv, so a stack can add up to slightly more
  # than 1; the range is set with coord_cartesian, which (unlike scale limits)
  # does not drop a bar that goes past 1
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.04))) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Proportion of explained variance") +
  theme_s6
p_vp

p_rho <- ggplot(rho, aes(x = response, y = rho_mean)) +
  geom_linerange(aes(ymin = rho_q05, ymax = rho_q95), linewidth = 0.8, colour = "grey20") +
  geom_point(size = 2.8, colour = "grey20") +
  facet_wrap(~ stratum, nrow = 1) +
  # rho is bounded at 1, so intervals that reach 1 end there; the extra space
  # above 1 shows that they are not cut by the axis
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = expression("Phylogenetic signal ("*rho*")")) +
  theme_s6
p_rho

fig <- (p_vp | p_rho) + plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = "A")
fig_file <- file.path(fig_dir, paste0(fig_prefix, "figS6_1_variance_partitioning_", tag, ".png"))
ggsave(fig_file, fig, width = 12, height = 5.5, dpi = 300, bg = "white")
cat("Figure written:", fig_file, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)