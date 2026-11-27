# =============================================================================
# 16_baseline_affinity_mismatch.R
#
# Pipeline: pipeline/02_analysis/01_thermophilization/03_integration_over_under/code/
# Author: Agostina Torres
#
# Descriptive comparison of community climatic affinity and local climate at
# the start of the series (past survey). If communities were in equilibrium
# with climate, past CTI would track past MAT one to one (and CMI would track
# MAP). The difference index − climate is shown along the elevational gradient
# together with the 1:1 line and the fitted past relationship (the space-for-
# time relationship used for the expected rates).
#
# The absolute difference is not a direct estimate of an inherited climatic
# debt: species' affinities are means over their whole ranges, so community
# indices are compressed towards the regional climate (slope < 1) even at
# equilibrium. Separating disequilibrium from this compression would require a
# calibration on independent plots assumed to be in equilibrium.
#
# Inputs : <stratum>/outputs/<stratum>_thermo_xero_rates.csv (past indices)
#          01_data_prep/03_climate/02_explore_climate/data_clean/climate_past_present_wide.csv
#          01_data_prep/01_cleaning/outputs/plot_info.csv (elevation)
# Outputs: 03_integration_over_under/outputs/baseline_affinity_mismatch.csv
#          03_integration_over_under/outputs/baseline_affinity_mismatch_summary.csv
#          03_integration_over_under/figures/baseline_affinity_mismatch.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

thermo_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", ...)
out_dir <- thermo_dir("03_integration_over_under", "outputs")
fig_dir <- thermo_dir("03_integration_over_under", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

STRATA <- tibble(stratum = c("overstory", "understory"), dir = c("01_overstory", "02_understory"),
                 weight = c("basal_area_spp_perc", "cover"))
# Thermal panels in coral (CTI), moisture panels in blue (CMI); layers as dark/light shades
col_T <- c(Overstory = "#B8402F", Understory = "#F7A194")
col_P <- c(Overstory = "#0B3C5D", Understory = "#6FB3D2")

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.text = element_text(colour = "black", size = 11),
          axis.title.x = element_text(colour = "black", size = 12, hjust = 0.5),
          axis.title.y = element_text(colour = "black", size = 12, hjust = 0.5),
          plot.title = element_text(colour = "black", size = 13, face = "bold"),
          legend.title = element_blank(), legend.position = "bottom")
}
first_col <- function(d, patterns) {
  for (p in patterns) { hit <- grep(p, names(d), value = TRUE, ignore.case = TRUE); if (length(hit)) return(hit[1]) }
  NA_character_
}

# -----------------------------------------------------------------------------
# 1. Past indices, past climate and elevation
# -----------------------------------------------------------------------------
clim <- read_csv(here("pipeline", "01_data_prep", "03_climate", "02_explore_climate", "data_clean",
                      "climate_past_present_wide.csv"), show_col_types = FALSE) %>%
  transmute(plot_id = as.character(plot_id), MAT_past, MAP_past)
info <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", "plot_info.csv"), show_col_types = FALSE)
elev <- info %>% transmute(plot_id = as.character(.data[[first_col(info, c("^plot_id$", "^plot$"))]]),
                           elevation = .data[[first_col(info, c("^elevation$", "^elev"))]]) %>%
  distinct(plot_id, .keep_all = TRUE)

plots <- pmap_dfr(STRATA, function(stratum, dir, weight) {
  d <- read_csv(thermo_dir(dir, "outputs", paste0(stratum, "_thermo_xero_rates.csv")), show_col_types = FALSE)
  if ("weight_var" %in% names(d)) d <- filter(d, weight_var == weight)
  cti <- first_col(d, c("^CTI_past$", "^CTI\\.past$", "^CTI_(median|q50|med)\\.?_?past$"))
  cmi <- first_col(d, c("^CMI_past$", "^CMI\\.past$", "^CMI_(median|q50|med)\\.?_?past$"))
  if (is.na(cti) || is.na(cmi)) stop("Past CTI/CMI not found for ", stratum, ": ", paste(names(d), collapse = ", "))
  d %>% transmute(plot_id = as.character(plot_id), stratum = str_to_title(stratum),
                  CTI_past = .data[[cti]], CMI_past = .data[[cmi]]) %>% distinct(plot_id, .keep_all = TRUE)
}) %>%
  left_join(clim, by = "plot_id") %>% left_join(elev, by = "plot_id") %>%
  mutate(mismatch_T = CTI_past - MAT_past, mismatch_P = CMI_past - MAP_past)
write_csv(plots, file.path(out_dir, "baseline_affinity_mismatch.csv"))

# -----------------------------------------------------------------------------
# 2. Summary: mean difference, past relationship and elevation where index = climate
# -----------------------------------------------------------------------------
summ <- plots %>% group_by(stratum) %>%
  group_modify(function(d, k) {
    fT <- lm(CTI_past ~ MAT_past, d); fP <- lm(CMI_past ~ MAP_past, d)
    tibble(n = nrow(d),
           mismatch_T_mean = mean(d$mismatch_T, na.rm = TRUE), mismatch_T_sd = sd(d$mismatch_T, na.rm = TRUE),
           slope_T = coef(fT)[2], MAT_equal = coef(fT)[1] / (1 - coef(fT)[2]),      # MAT at which CTI = MAT
           mismatch_P_mean = mean(d$mismatch_P, na.rm = TRUE), mismatch_P_sd = sd(d$mismatch_P, na.rm = TRUE),
           slope_P = coef(fP)[2], MAP_equal = coef(fP)[1] / (1 - coef(fP)[2]),      # MAP at which CMI = MAP
           elev_slope_T_per_km = coef(lm(mismatch_T ~ I(elevation / 1000), d))[2],
           elev_slope_P_per_km = coef(lm(mismatch_P ~ I(elevation / 1000), d))[2])
  }) %>% ungroup() %>% mutate(across(where(is.double), ~ signif(.x, 3)))
write_csv(summ, file.path(out_dir, "baseline_affinity_mismatch_summary.csv"))
cat("\nIndex − climate in the past survey (CTI − MAT, °C; CMI − MAP, mm):\n"); print(summ, width = Inf)

# -----------------------------------------------------------------------------
# 3. Figure
# -----------------------------------------------------------------------------
scatter <- function(x, y, xl, yl, title, cols) {
  ggplot(plots, aes(.data[[x]], .data[[y]], colour = stratum)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey40") +
    geom_point(alpha = 0.45, size = 1.3) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.9) +
    scale_colour_manual(values = cols) +
    labs(title = title, x = xl, y = yl) + theme_mora() +
    theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1))
}
along <- function(y, yl, title, cols) {
  ggplot(plots, aes(elevation, .data[[y]], colour = stratum)) +
    geom_hline(yintercept = 0, linetype = "22", colour = "grey40") +
    geom_point(alpha = 0.45, size = 1.3) +
    geom_smooth(aes(fill = stratum), method = "loess", formula = y ~ x, span = 0.75, se = TRUE, linewidth = 0.9, alpha = 0.15) +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(title = title, x = "Elevation (m)", y = yl) + theme_mora() +
    guides(colour = "none", fill = "none")                                   # legend shown in panels A-B
}
fig <- (scatter("MAT_past", "CTI_past", "Past MAT (°C)", "Past CTI (°C)", "A. Thermal", col_T) |
          scatter("MAP_past", "CMI_past", "Past MAP (mm)", "Past CMI (mm)", "B. Moisture", col_P)) /
  (along("mismatch_T", "CTI − MAT (°C)", "C. Thermal difference", col_T) |
     along("mismatch_P", "CMI − MAP (mm)", "D. Moisture difference", col_P))
fig_file <- file.path(fig_dir, "baseline_affinity_mismatch.png")
ggsave(fig_file, fig, width = 11, height = 8.5, dpi = 300, bg = "white")
cat("\nFigure written to", fig_file, "\n")

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
