# ------------------------------------------------------------------------------
# 03_cti_cmi_past_present_elevation.R
# Module: pipeline/02_analysis/01_thermophilization/02_understory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: descriptive figures showing the community temperature index (CTI)
# and community moisture index (CMI) directly against elevation, past vs
# present, rather than as a rate of change (that comparison is in
# 02_thermo_xero_elevation_lag.R). This shows the raw elevational pattern of
# each index in both survey periods side by side: a vertical gap between the
# two curves at a given elevation is the community-level shift in thermal or
# moisture affinity at that point on the gradient. Understory counterpart of
# 01_overstory/code/03_CTI_CMI_past_present_elevation.R.
#
# Uses the primary weight variable (cover, see 01_compute_cti_cmi_understory.R) and the
# bootstrap median (CTI_median/CMI_median), for consistency with
# 02_thermo_xero_elevation_lag.R.
#
# Output:
#   - figures/understory_CTI_past_present_elevation.png
#   - figures/understory_CMI_past_present_elevation.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)

# 0a. Module directories --------------------------------------------------------
this_module_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "02_understory", ...)

dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

PRIMARY_WEIGHT <- "cover"

# Colors: past = navy/marine blue, present = coral red (same convention as overstory)
PAST_COLOR    <- "#0A2463"
PRESENT_COLOR <- "#FF6B5B"

# 1. Load ------------------------------------------------------------------------
cti_cmi_quantiles <- read_csv(this_module_dir("data_clean", "understory_CTI_CMI_quantiles.csv"), show_col_types = FALSE)

cq <- cti_cmi_quantiles %>%
  filter(weight_var == PRIMARY_WEIGHT) %>%
  mutate(period = factor(period, levels = c("past", "present"), labels = c("Past (~1975)", "Present (~2015)")))

# 2. Plot function: index vs elevation, past/present LOESS + points + CI --------
plot_index_by_period <- function(df, median_col, q05_col, q95_col, y_lab, title_lab) {
  ggplot(df, aes(x = elevation, y = .data[[median_col]], color = period, fill = period)) +
    geom_point(alpha = 0.35, size = 1.6) +
    geom_smooth(method = "loess", span = 0.75, se = TRUE, linewidth = 1.1, alpha = 0.15) +
    scale_color_manual(values = c("Past (~1975)" = PAST_COLOR, "Present (~2015)" = PRESENT_COLOR), name = NULL) +
    scale_fill_manual(values = c("Past (~1975)" = PAST_COLOR, "Present (~2015)" = PRESENT_COLOR), name = NULL) +
    labs(x = "Elevation (m)", y = y_lab, title = title_lab) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 15, face = "plain"),
      legend.position = "top"
    )
}

p_CTI <- plot_index_by_period(
  cq, "CTI_median", "CTI_q05", "CTI_q95",
  "Community temperature index (°C)",
  "Understory CTI vs elevation, past vs present"
)
p_CTI

p_CMI <- plot_index_by_period(
  cq, "CMI_median", "CMI_q05", "CMI_q95",
  "Community moisture index (mm)",
  "Understory CMI vs elevation, past vs present"
)
p_CMI

# 3. Save (white background, per project convention) -----------------------------
ggsave(this_module_dir("figures", "understory_CTI_past_present_elevation.png"), p_CTI,
       width = 8, height = 6, dpi = 300, bg = "white")

ggsave(this_module_dir("figures", "understory_CMI_past_present_elevation.png"), p_CMI,
       width = 8, height = 6, dpi = 300, bg = "white")

cat("Saved:\n")
cat(" - figures/understory_CTI_past_present_elevation.png\n")
cat(" - figures/understory_CMI_past_present_elevation.png\n")
