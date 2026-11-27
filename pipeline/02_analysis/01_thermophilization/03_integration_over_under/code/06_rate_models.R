# ------------------------------------------------------------------------------
# 06_rate_models.R
# Module: pipeline/02_analysis/01_thermophilization/03_integration_over_under
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: formal tests and distributions of the plot-level rates of change in
# the community temperature index (CTI, thermophilization) and community
# moisture index (CMI, moisture-affinity shift), for overstory and understory
# (primary weight of each layer: basal_area_spp_perc and cover, as in 01-04).
# Complements the descriptive summaries and LOESS curves of 01-04. Reads each
# layer's *_thermo_xero_rates.csv and *_thermo_xero_expected.csv; does not
# recompute indices.
#
# MODELS (Gaussian errors: rates are continuous and approximately symmetric,
# so the Gaussian GLM is a linear model; predictors are centred and scaled, so
# slopes are per 1 SD of the predictor; elevation is scaled once across
# plots, so slopes are comparable between strata):
#   A. Mean rate:                 rate ~ 1                      (per stratum)
#   B. Rate vs elevation:         rate ~ elevation              (per stratum)
#   C. Rate vs measured climate-change rate (MAT trend for CTI, MAP trend for
#      CMI):                      rate ~ climate rate           (per stratum)
#   D. Climatic lag vs elevation: lag ~ elevation               (per stratum)
#   E. Both layers:               rate ~ stratum * elevation + (1 | plot_id)
#      (overstory and understory are surveyed in the same permanent plots, so
#      each plot contributes two non-independent rows)
# Residual spatial autocorrelation is assessed for every model with Moran's I
# on the residuals (inverse-distance weights on projected coordinates, UTM 10N, km).
#
# DISTRIBUTIONS: 2x2 histograms of the plot-level rates (CTI, CMI; overstory,
# understory) with the mean rate and zero (no change) marked. The x range is
# shared by the two layers within a row so distributions are comparable.
#
# STYLE: theme_ipsum() with the grid blanked, finite black axis lines, color by
# variable (coral for CTI, ocean blue for CMI), "Overstory"/"Understory" shown
# once as column headers, one shared legend above the grid.
#
# Output:
#   - outputs/rate_models_coefficients.csv: fixed effects of every model
#   - outputs/rate_models_moran.csv: Moran's I and p-value per model
#   - figures/rate_distributions_overstory_understory.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(lme4)
library(lmerTest)
library(ape)
library(patchwork)
library(hrbrthemes)
library(cowplot)      # get_legend()

thermo_dir      <- function(...) here("pipeline", "02_analysis", "01_thermophilization", ...)
this_module_dir <- function(...) thermo_dir("03_integration_over_under", ...)
dir.create(this_module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

STRATA <- tibble(
  stratum = c("overstory", "understory"),
  dir     = c("01_overstory", "02_understory"),
  weight  = c("basal_area_spp_perc", "cover")
)

# 1. Load per-plot rates and expected-rate tables ------------------------------
load_stratum <- function(st, dir, weight) {
  rates <- read_csv(thermo_dir(dir, "outputs", paste0(st, "_thermo_xero_rates.csv")), show_col_types = FALSE) %>%
    filter(weight_var == weight)
  expected <- read_csv(thermo_dir(dir, "data_clean", paste0(st, "_thermo_xero_expected.csv")), show_col_types = FALSE)
  if ("weight_var" %in% names(expected)) expected <- filter(expected, weight_var == weight)
  
  rates %>%
    left_join(select(expected, plot_id, lag_CTI, lag_CMI), by = "plot_id") %>%
    mutate(stratum = st)
}

plots <- pmap_dfr(STRATA, ~ load_stratum(..1, ..2, ..3))

# Elevation centred and scaled once across plots (each plot counted once), so
# that slopes are comparable between layers and in the joint model
elev_ref <- plots %>% distinct(plot_id, elevation)
elev_mu  <- mean(elev_ref$elevation, na.rm = TRUE)
elev_sd  <- sd(elev_ref$elevation, na.rm = TRUE)
plots <- plots %>% mutate(elev_sc = (elevation - elev_mu) / elev_sd)

# Projected coordinates (UTM zone 10N, km) for the distance weights
xy <- plots %>%
  distinct(plot_id, lat, long) %>%
  filter(!is.na(lat), !is.na(long)) %>%
  sf::st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(32610)
xy <- tibble(plot_id = xy$plot_id,
             x_km = sf::st_coordinates(xy)[, 1] / 1000,
             y_km = sf::st_coordinates(xy)[, 2] / 1000)
plots <- plots %>% left_join(xy, by = "plot_id")
cat("Plots per stratum:\n"); print(count(plots, stratum))

# 2. Helpers -------------------------------------------------------------------
coef_table <- function(model, label) {
  as_tibble(coef(summary(model)), rownames = "term") %>%
    mutate(model = label, .before = 1)
}

# Moran's I on the residuals of a fitted model. Rows dropped by the model
# because of missing values are matched back to their coordinates through the
# residual names (row indices of the data).
moran_resid <- function(fit, data) {
  r   <- residuals(fit)
  idx <- as.integer(names(r))
  if (any(is.na(idx))) idx <- seq_along(r)
  d <- data[idx, ] %>% mutate(r = as.numeric(r)) %>% filter(!is.na(x_km), !is.na(y_km))
  w <- 1 / as.matrix(dist(cbind(d$x_km, d$y_km)))
  diag(w) <- 0
  w[!is.finite(w)] <- 0     # co-located plots (e.g. both layers of one plot) get zero weight
  mi <- Moran.I(d$r, w)
  tibble(moran_I = mi$observed, moran_p = mi$p.value)
}

all_coefs <- list()
all_moran <- list()

run_model <- function(label, fit, data) {
  cat("\n---", label, "---\n")
  print(coef_table(fit, label), width = Inf)
  cat("AIC:", round(AIC(fit), 1), "\n")
  mo <- moran_resid(fit, data)
  cat("Moran's I:", signif(mo$moran_I, 3), "| p-value:", signif(mo$moran_p, 3), "\n")
  all_coefs[[label]] <<- coef_table(fit, label)
  all_moran[[label]] <<- mutate(mo, model = label, .before = 1)
  invisible(fit)
}

# 3. Models per layer and variable -------------------------------------------
for (st in STRATA$stratum) {
  d <- plots %>% filter(stratum == st)
  cat("\n==============================\n", toupper(st), "\n==============================\n", sep = "")
  
  for (v in c("CTI", "CMI")) {
    rate    <- paste0("rate_", v)
    lag     <- paste0("lag_", v)
    climate <- if (v == "CTI") "climate_rate_MAT" else "climate_rate_MAP"
    d$clim_sc <- as.numeric(scale(d[[climate]]))
    
    cat("\n########", v, "########\n")
    
    run_model(paste(st, v, "A mean rate"),
              lm(as.formula(paste(rate, "~ 1")), data = d), d)
    run_model(paste(st, v, "B elevation"),
              lm(as.formula(paste(rate, "~ elev_sc")), data = d), d)
    run_model(paste(st, v, "C climate change rate"),
              lm(as.formula(paste(rate, "~ clim_sc")), data = d), d)
    run_model(paste(st, v, "D lag ~ elevation"),
              lm(as.formula(paste(lag, "~ elev_sc")), data = d), d)
  }
}

# 4. Both layers together ------------------------------------------------------
cat("\n==============================\nBOTH STRATA\n==============================\n")
plots_both <- plots %>% mutate(stratum = factor(stratum, levels = c("overstory", "understory")))

for (v in c("CTI", "CMI")) {
  rate <- paste0("rate_", v)
  m <- lmer(as.formula(paste(rate, "~ stratum * elev_sc + (1 | plot_id)")), data = plots_both)
  run_model(paste("both", v, "E stratum x elevation"), m, plots_both)
}

# 5. Save coefficients and spatial-autocorrelation tests -----------------------
write_csv(bind_rows(all_coefs), this_module_dir("outputs", "rate_models_coefficients.csv"))
write_csv(bind_rows(all_moran), this_module_dir("outputs", "rate_models_moran.csv"))

# 6. Histograms of plot-level rates --------------------------------------------
N_BINS <- 30

HEADER_OVERSTORY  <- "\U0001F332 Overstory"
HEADER_UNDERSTORY <- "\U0001F33F Understory"

THERMO_COLOR <- "#F26B5B"  # coral: temperature / CTI
XERO_COLOR   <- "#1A759F"  # ocean blue: precipitation / CMI

darken_color <- function(col, amount = 0.4) {
  m <- col2rgb(col) * (1 - amount)
  rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

base_theme <- theme_ipsum(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line.x      = element_line(color = "black", linewidth = 0.3),
    axis.line.y      = element_line(color = "black", linewidth = 0.3),
    axis.title.x     = element_text(size = 14, hjust = 0.5),
    axis.title.y     = element_text(size = 14),
    axis.text        = element_text(size = 12),
    legend.position  = "none",
    plot.title       = element_blank()
  )

overstory_rates  <- filter(plots, stratum == "overstory")
understory_rates <- filter(plots, stratum == "understory")

# x_limits is shared by the two layers of a row; bins have the same width.
# coord_cartesian() clips the view without dropping bins.
build_hist_panel <- function(df, rate_col, base_color, x_lab, x_limits, digits) {
  x <- df[[rate_col]]
  x <- x[!is.na(x)]
  bin_width <- diff(x_limits) / N_BINS
  mean_x <- mean(x)
  label  <- sprintf("Mean = %s\n%.0f%% of plots > 0", formatC(mean_x, format = "f", digits = digits), 100 * mean(x > 0))
  
  ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(binwidth = bin_width, boundary = 0, fill = base_color, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = mean_x, color = darken_color(base_color), linewidth = 0.9, linetype = "dashed") +
    annotate("text", x = x_limits[2], y = Inf, label = label, hjust = 1, vjust = 1.3, size = 4.4,
             color = "grey20", lineheight = 0.95) +
    coord_cartesian(xlim = x_limits) +
    labs(x = x_lab, y = "Number of plots") +
    base_theme
}

# Shared x range per row: widest range across the two layers, padded slightly
shared_limits <- function(col) {
  rng <- range(c(overstory_rates[[col]], understory_rates[[col]]), na.rm = TRUE)
  pad <- 0.03 * diff(rng)
  c(rng[1] - pad, rng[2] + pad)
}

CTI_LIMITS <- shared_limits("rate_CTI")
CMI_LIMITS <- shared_limits("rate_CMI")

x_lab_CTI <- "Thermophilization rate (°C/decade)"
x_lab_CMI <- "Moisture-affinity shift rate (mm/decade)"

p_CTI_over  <- build_hist_panel(overstory_rates,  "rate_CTI", THERMO_COLOR, x_lab_CTI, CTI_LIMITS, 3)
p_CTI_under <- build_hist_panel(understory_rates, "rate_CTI", THERMO_COLOR, x_lab_CTI, CTI_LIMITS, 3)
p_CMI_over  <- build_hist_panel(overstory_rates,  "rate_CMI", XERO_COLOR,   x_lab_CMI, CMI_LIMITS, 1)
p_CMI_under <- build_hist_panel(understory_rates, "rate_CMI", XERO_COLOR,   x_lab_CMI, CMI_LIMITS, 1)

# Shared legend: mean (dashed) and no change (solid)
legend_levels <- c("Mean rate", "No change (0)")

legend_dummy_df <- tibble(
  x    = rep(c(1, 2), 2),
  y    = rep(c(1, 2), 2),
  type = factor(rep(legend_levels, each = 2), levels = legend_levels)
)

p_legend_dummy <- ggplot(legend_dummy_df, aes(x = x, y = y, color = type, linetype = type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Mean rate" = "grey20", "No change (0)" = "grey50"), name = NULL) +
  scale_linetype_manual(values = c("Mean rate" = "dashed", "No change (0)" = "solid"), name = NULL) +
  theme_ipsum(base_size = 15) +
  theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 15))

shared_legend <- get_legend(p_legend_dummy)

header_overstory  <- wrap_elements(grid::textGrob(HEADER_OVERSTORY,  gp = grid::gpar(fontsize = 20, fontface = "bold")))
header_understory <- wrap_elements(grid::textGrob(HEADER_UNDERSTORY, gp = grid::gpar(fontsize = 20, fontface = "bold")))
column_headers <- header_overstory | header_understory

p_grid <- (p_CTI_over | p_CTI_under) / (p_CMI_over | p_CMI_under)

p_final <- wrap_elements(shared_legend) / column_headers / p_grid +
  plot_layout(heights = c(0.06, 0.07, 1))
print(p_final)

ggsave(this_module_dir("figures", "rate_distributions_overstory_understory.png"), p_final,
       width = 12, height = 10.5, dpi = 300, bg = "white")

# 7. Summary -------------------------------------------------------------------
cat("\n========================================\n")
cat("RATE MODELS AND DISTRIBUTIONS -- COMPLETE\n")
cat("========================================\n")
cat("Saved:\n")
cat(" - outputs/rate_models_coefficients.csv\n")
cat(" - outputs/rate_models_moran.csv\n")
cat(" - figures/rate_distributions_overstory_understory.png\n")
