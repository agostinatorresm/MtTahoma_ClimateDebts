# ------------------------------------------------------------------------------
# 02_affinity_edges_integration.R
# Module: pipeline/02_analysis/01_thermophilization/03_integration_over_under
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Overstory and understory community tails in one 2x2 figure: mean rate of
# change of the cold (dry) tail, weighted median and warm (wet) tail of species
# affinities, with 95% confidence intervals across plots (outputs of
# 06_community_affinity_edges.R in each stratum).
#   Row 1: CTI (thermal affinity)  -- overstory (left), understory (right)
#   Row 2: CMI (moisture affinity) -- overstory (left), understory (right)
# The p-values of the paired comparison between the cold/dry and warm/wet tail
# rates are printed to the console and reported in the figure caption.
#
# STYLE: as in 01, 03, 04 and 06 of this module -- theme_ipsum() with the grid
# blanked and a white background, finite black axis lines, colour by variable
# (coral for CTI, ocean blue for CMI, the same in both columns),
# "Overstory"/"Understory" shown once as column headers.
#
# Inputs : 01_overstory/outputs/overstory_affinity_edge_summary.csv
#          02_understory/outputs/understory_affinity_edge_summary.csv
# Outputs: figures/affinity_edges_over_under.png
#          outputs/affinity_edges_over_under.csv
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(patchwork)
library(hrbrthemes)

thermo_dir      <- function(...) here("pipeline", "02_analysis", "01_thermophilization", ...)
this_module_dir <- function(...) thermo_dir("03_integration_over_under", ...)
for (d in c("outputs", "figures")) dir.create(this_module_dir(d), recursive = TRUE, showWarnings = FALSE)

STRATA <- c(overstory = "01_overstory", understory = "02_understory")

THERMO_COLOR <- "#F26B5B"  # coral: temperature / CTI
XERO_COLOR   <- "#1A759F"  # ocean blue: precipitation / CMI

HEADER_OVERSTORY  <- "\U0001F332 Overstory"
HEADER_UNDERSTORY <- "\U0001F33F Understory"

# One shared set of x labels (bottom row only), on two lines so they do not overlap
EDGE_LABELS <- c(q10 = "Cold /\ndry tail", q50 = "Median", q90 = "Warm /\nwet tail")
Y_LABELS    <- c(CTI = "Thermophilization rate (°C/decade)",
                 CMI = "Moisture-affinity shift (mm/decade)")

# 1. Load ---------------------------------------------------------------------
summary_all <- imap_dfr(STRATA, function(dir, st) {
  f <- thermo_dir(dir, "outputs", paste0(st, "_affinity_edge_summary.csv"))
  if (!file.exists(f)) stop("Run ", dir, "/code/06_community_affinity_edges.R first: ", f)
  read_csv(f, show_col_types = FALSE) %>% mutate(stratum = st)
})
write_csv(summary_all, this_module_dir("outputs", "affinity_edges_over_under.csv"))
print(summary_all, width = Inf)
cat("\nPaired test, cold/dry vs warm/wet edge rates (for the caption):\n")
print(summary_all %>% distinct(stratum, index, p_q10_vs_q90))

# 2. Theme ----------------------------------------------------------------------
if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
finite_axis_theme <- theme_ipsum(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line.x      = element_line(color = "black", linewidth = 0.3),
    axis.line.y      = element_line(color = "black", linewidth = 0.3),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 15),
    axis.text        = element_text(size = 13),
    legend.position  = "none",
    plot.title       = element_blank()
  )

# 3. Panel ----------------------------------------------------------------------
# Same y range within a row so the two layers are directly comparable
row_limits <- summary_all %>%
  group_by(index) %>%
  summarise(lo = min(ci_lo, 0), hi = max(ci_hi, 0), .groups = "drop") %>%
  mutate(pad = 0.05 * (hi - lo), lo = lo - 0.03 * (hi - lo), hi = hi + pad)

build_edge_panel <- function(st, ix, base_color) {
  df <- summary_all %>%
    filter(stratum == st, index == ix) %>%
    mutate(edge_lab = factor(EDGE_LABELS[edge], levels = EDGE_LABELS))
  lim <- row_limits %>% filter(index == ix)
  
  p <- ggplot(df, aes(x = edge_lab, y = mean, ymin = ci_lo, ymax = ci_hi)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3, linetype = "dashed") +
    geom_linerange(color = base_color, linewidth = 1) +
    geom_point(color = base_color, size = 3.5) +
    coord_cartesian(ylim = c(lim$lo, lim$hi)) +
    labs(x = NULL, y = Y_LABELS[[ix]]) +
    finite_axis_theme
  # x labels only in the bottom (CMI) row
  if (ix == "CTI") p <- p + theme(axis.text.x = element_blank())
  p
}

p_CTI_over  <- build_edge_panel("overstory",  "CTI", THERMO_COLOR)
p_CTI_under <- build_edge_panel("understory", "CTI", THERMO_COLOR)
p_CMI_over  <- build_edge_panel("overstory",  "CMI", XERO_COLOR)
p_CMI_under <- build_edge_panel("understory", "CMI", XERO_COLOR)

# 4. Column headers + grid ---------------------------------------------------------
column_headers <- wrap_elements(grid::textGrob(HEADER_OVERSTORY,  gp = grid::gpar(fontsize = 20, fontface = "bold"))) |
  wrap_elements(grid::textGrob(HEADER_UNDERSTORY, gp = grid::gpar(fontsize = 20, fontface = "bold")))

p_grid <- (p_CTI_over | p_CTI_under) / (p_CMI_over | p_CMI_under)

fig <- column_headers / p_grid + plot_layout(heights = c(0.07, 1))
print(fig)

ggsave(this_module_dir("figures", "affinity_edges_over_under.png"), fig,
       width = 12, height = 10.5, dpi = 300, bg = "white")
cat("Saved: figures/affinity_edges_over_under.png\n")
