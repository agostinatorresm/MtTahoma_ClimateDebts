# ------------------------------------------------------------------------------
# 03_CTI_CMI_past_present_elevation.R
# Module: pipeline/02_analysis/01_thermophilization/01_overstory
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: descriptive figures showing the community temperature index (CTI)
# and community moisture index (CMI) directly against elevation, past vs
# present, rather than as a rate of change (that comparison is in
# 02_thermo_xero_elevation_lag.R). This shows the raw elevational pattern of
# each index in both survey periods side by side: a vertical gap between the
# two curves at a given elevation is the community-level shift in thermal or
# moisture affinity at that point on the gradient.
#
# Uses the primary weight variable (basal_area_spp_perc, see
# 01_compute_cti_cmi.R) and the bootstrap median (CTI_median/CMI_median), for
# consistency with 02_thermo_xero_elevation_lag.R.
#
# Output:
#   - figures/overstory_CTI_past_present_elevation.png
#   - figures/overstory_CMI_past_present_elevation.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)

# 0a. Module directories --------------------------------------------------------
this_module_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", ...)

dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

PRIMARY_WEIGHT <- "basal_area_spp_perc"

# Colors: past = navy/marine blue, present = coral red
PAST_COLOR    <- "#0A2463"
PRESENT_COLOR <- "#FF6B5B"

# 1. Load ------------------------------------------------------------------------
cti_cmi_quantiles <- read_csv(this_module_dir("data_clean", "overstory_CTI_CMI_quantiles.csv"), show_col_types = FALSE)

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
  "Overstory CTI vs elevation, past vs present"
)
p_CTI

p_CMI <- plot_index_by_period(
  cq, "CMI_median", "CMI_q05", "CMI_q95",
  "Community moisture index (mm)",
  "Overstory CMI vs elevation, past vs present"
)
p_CMI

# 3. Save (white background, per project convention) -----------------------------
ggsave(this_module_dir("figures", "overstory_CTI_past_present_elevation.png"), p_CTI,
       width = 8, height = 6, dpi = 300, bg = "white")

ggsave(this_module_dir("figures", "overstory_CMI_past_present_elevation.png"), p_CMI,
       width = 8, height = 6, dpi = 300, bg = "white")

cat("Saved:\n")
cat(" - figures/overstory_CTI_past_present_elevation.png\n")
cat(" - figures/overstory_CMI_past_present_elevation.png\n")
