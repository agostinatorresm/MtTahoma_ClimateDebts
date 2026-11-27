# =============================================================================
# 15_space_for_time_and_affinity_correlations.R
#
# Pipeline: pipeline/02_analysis/01_thermophilization/03_integration_over_under/code/
# Author: Agostina Torres
#
# Two diagnostics of the climatic-affinity approach:
#
# 1. Space-for-time relationships. The expected rates (and the climatic lags)
#    rely on the relationship between past community indices and past climate
#    across plots (CTI ~ MAT, CMI ~ MAP). A weak relationship makes the
#    expected rate, and therefore the lag, poorly constrained. For each layer,
#    the correlation (Pearson and Spearman), the slope and the R2 of the linear
#    relationship used for the expected rates are reported.
#
# 2. Collinearity of species' affinities. If species' thermal (MAT) and
#    moisture (MAP) affinities are correlated, CTI and CMI are not independent
#    and their shifts (and the HMSC trait effects of the two affinities) cannot
#    be fully separated. The correlation between species' point affinities in
#    each species pool, and the correlations between plot-level CTI and CMI
#    (past) and between their rates of change, are reported.
#
# Inputs : <stratum>/data_clean/<stratum>_thermo_xero_expected.csv (past indices, if present)
#          <stratum>/outputs/<stratum>_thermo_xero_rates.csv
#          01_data_prep/03_climate/02_explore_climate/data_clean/climate_past_present_wide.csv
#          species affinity table (first csv with species, MAT_mean and MAP_mean found
#          under pipeline/01_data_prep; the file used is printed)
#          01_data_prep/02_species_pool outputs (species pools, if found)
# Outputs: 03_integration_over_under/outputs/space_for_time_correlations.csv
#          03_integration_over_under/outputs/affinity_correlations.csv
#          03_integration_over_under/figures/space_for_time_relationships.png
#          03_integration_over_under/figures/affinity_correlation.png
#          03_integration_over_under/figures/cti_cmi_covariation.png
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

STRATA <- tibble(stratum = c("overstory", "understory"),
                 dir     = c("01_overstory", "02_understory"),
                 weight  = c("basal_area_spp_perc", "cover"))
col_cti <- "#F26B5B"; col_cmi <- "#1A759F"

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.text = element_text(colour = "black", size = 11),
          axis.title.x = element_text(colour = "black", size = 12, hjust = 0.5),
          axis.title.y = element_text(colour = "black", size = 12, hjust = 0.5),
          plot.title = element_text(colour = "black", size = 13, face = "bold"))
}

# Correlation summary: Pearson r with 95% CI, Spearman rho, and the linear slope and R2
cor_summary <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  pr <- cor.test(x, y, method = "pearson")
  sp <- suppressWarnings(cor.test(x, y, method = "spearman"))
  fit <- lm(y ~ x)
  tibble(n = length(x), pearson_r = unname(pr$estimate),
         r_lower_95 = pr$conf.int[1], r_upper_95 = pr$conf.int[2], pearson_p = pr$p.value,
         spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
         slope = unname(coef(fit)[2]), r2 = summary(fit)$r.squared)
}
first_col <- function(d, patterns) {
  for (p in patterns) { hit <- grep(p, names(d), value = TRUE, ignore.case = TRUE); if (length(hit)) return(hit[1]) }
  NA_character_
}

# -----------------------------------------------------------------------------
# 1. Plot-level past indices and past climate
# -----------------------------------------------------------------------------
clim <- read_csv(here("pipeline", "01_data_prep", "03_climate", "02_explore_climate", "data_clean",
                      "climate_past_present_wide.csv"), show_col_types = FALSE) %>%
  transmute(plot_id = as.character(plot_id), MAT_past, MAP_past)

load_plots <- function(st, dir, weight) {
  rates <- read_csv(thermo_dir(dir, "outputs", paste0(st, "_thermo_xero_rates.csv")), show_col_types = FALSE)
  if ("weight_var" %in% names(rates)) rates <- filter(rates, weight_var == weight)
  f_exp <- thermo_dir(dir, "data_clean", paste0(st, "_thermo_xero_expected.csv"))
  expd <- if (file.exists(f_exp)) read_csv(f_exp, show_col_types = FALSE) else tibble(plot_id = character())
  if ("weight_var" %in% names(expd)) expd <- filter(expd, weight_var == weight)
  d <- rates %>% mutate(plot_id = as.character(plot_id)) %>%
    left_join(mutate(expd, plot_id = as.character(plot_id)), by = "plot_id", suffix = c("", ".exp"))
  # past indices and rates (column names differ between versions; the names used are printed)
  # the community index is the median over the 100 affinity resamples (not the
  # 5% or 95% quantiles, which are also stored)
  cols <- c(CTI_past = first_col(d, c("^CTI_past$", "^CTI\\.past$", "^CTI_(median|q50|med)\\.?_?past$", "^CTI_(median|q50|med)")),
            CMI_past = first_col(d, c("^CMI_past$", "^CMI\\.past$", "^CMI_(median|q50|med)\\.?_?past$", "^CMI_(median|q50|med)")),
            rate_CTI = first_col(d, c("^rate_CTI$", "CTI.*rate", "rate.*CTI", "^delta_CTI")),
            rate_CMI = first_col(d, c("^rate_CMI$", "CMI.*rate", "rate.*CMI", "^delta_CMI")))
  cat("\n", st, ": columns used\n"); print(cols)
  if (any(is.na(cols[c("CTI_past", "CMI_past")])) || any(grepl("q05|q95|q025|q975", cols[c("CTI_past", "CMI_past")]))) {
    stop("Past CTI/CMI not found for ", st, ". Columns available:\n", paste(names(d), collapse = ", "))
  }
  d %>% transmute(plot_id, stratum = st,
                  CTI_past = .data[[cols["CTI_past"]]], CMI_past = .data[[cols["CMI_past"]]],
                  rate_CTI = if (!is.na(cols["rate_CTI"])) .data[[cols["rate_CTI"]]] else NA_real_,
                  rate_CMI = if (!is.na(cols["rate_CMI"])) .data[[cols["rate_CMI"]]] else NA_real_) %>%
    distinct(plot_id, .keep_all = TRUE) %>%
    left_join(clim, by = "plot_id")
}
plots <- pmap_dfr(STRATA, ~ load_plots(..1, ..2, ..3))

sft <- plots %>%
  group_by(stratum) %>%
  group_modify(~ bind_rows(
    mutate(cor_summary(.x$MAT_past, .x$CTI_past), relationship = "CTI ~ MAT (past)"),
    mutate(cor_summary(.x$MAP_past, .x$CMI_past), relationship = "CMI ~ MAP (past)"),
    mutate(cor_summary(.x$CTI_past, .x$CMI_past), relationship = "CMI ~ CTI (past, plots)"),
    if (all(is.na(.x$rate_CTI))) NULL else
      mutate(cor_summary(.x$rate_CTI, .x$rate_CMI), relationship = "CMI rate ~ CTI rate (plots)"),
    mutate(cor_summary(.x$MAT_past, .x$MAP_past), relationship = "MAP ~ MAT (past climate, plots)"))) %>%
  ungroup() %>%
  relocate(stratum, relationship) %>%
  mutate(across(where(is.numeric) & !n, ~ signif(.x, 3)))

write_csv(sft, file.path(out_dir, "space_for_time_correlations.csv"))
cat("\nSpace-for-time and plot-level correlations:\n")
print(select(sft, stratum, relationship, n, pearson_r, r_lower_95, r_upper_95, spearman_rho, slope, r2), n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 2. Species' thermal and moisture affinities
# -----------------------------------------------------------------------------
cand <- list.files(here("pipeline", "01_data_prep"), pattern = "affinit.*\\.csv$", recursive = TRUE,
                   full.names = TRUE, ignore.case = TRUE)
aff_file <- NA_character_
for (f in cand) {
  h <- names(read_csv(f, n_max = 1, show_col_types = FALSE))
  if (all(c("MAT_mean", "MAP_mean") %in% h) && any(c("species", "Species") %in% h)) { aff_file <- f; break }
}
if (is.na(aff_file)) stop("No species affinity table with species, MAT_mean and MAP_mean found under pipeline/01_data_prep")
cat("\nSpecies affinities read from:", aff_file, "\n")
aff <- read_csv(aff_file, show_col_types = FALSE) %>%
  rename(any_of(c(species = "Species"))) %>%
  mutate(species = str_replace_all(species, "_", " ")) %>%
  distinct(species, .keep_all = TRUE) %>%
  select(species, MAT_mean, MAP_mean)

# Species of each pool: shared pool from the species lists of each layer
# (long tables, species recorded at both times); restricted pool from the HMSC
# trait tables if found
long_sp <- function(st) {
  d <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", paste0(st, "_long.csv")), show_col_types = FALSE)
  tcol <- intersect(c("time", "period", "survey", "census"), names(d))[1]
  d %>% filter(occurrence > 0) %>% group_by(species) %>% filter(n_distinct(.data[[tcol]]) == 2) %>%
    ungroup() %>% distinct(species) %>% mutate(species = str_replace_all(species, "_", " "), stratum = st)
}
pools <- map_dfr(STRATA$stratum, long_sp) %>% mutate(pool = "shared")
# species of the restricted pool (modelled in HMSC, complete traits)
trait_files <- list.files(here("pipeline", "02_analysis", "02_hmsc", "outputs"),
                          pattern = "^(overstory|understory)_strict_diag_species_occurrence_traits\\.csv$", full.names = TRUE)
if (length(trait_files)) {
  cat("Restricted-pool species read from:\n", paste(trait_files, collapse = "\n "), "\n")
  restricted <- map_dfr(trait_files, function(f) {
    d <- read_csv(f, show_col_types = FALSE)
    scol <- intersect(c("species", "Species", "sp"), names(d))[1]
    if (is.na(scol)) return(NULL)
    tibble(species = str_replace_all(d[[scol]], "_", " "),
           stratum = if (str_detect(f, "overstory")) "overstory" else "understory")
  }) %>% distinct() %>% mutate(pool = "restricted")
  pools <- bind_rows(pools, restricted)
}
sp_aff <- pools %>% left_join(aff, by = "species")
missing_aff <- sp_aff %>% filter(is.na(MAT_mean) | is.na(MAP_mean))
if (nrow(missing_aff)) cat("Species without affinities (excluded):", nrow(missing_aff), "\n")
sp_aff <- sp_aff %>% filter(!is.na(MAT_mean), !is.na(MAP_mean))

affc <- sp_aff %>%
  group_by(stratum, pool) %>%
  group_modify(~ cor_summary(.x$MAT_mean, .x$MAP_mean)) %>%
  ungroup() %>%
  mutate(relationship = "MAP affinity ~ MAT affinity (species)") %>%
  relocate(stratum, pool, relationship) %>%
  mutate(across(where(is.numeric) & !n, ~ signif(.x, 3)))
write_csv(affc, file.path(out_dir, "affinity_correlations.csv"))
cat("\nCorrelation between species' thermal and moisture affinities:\n")
print(select(affc, stratum, pool, n, pearson_r, r_lower_95, r_upper_95, spearman_rho), n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 3. Figure
# -----------------------------------------------------------------------------
# Facet labels
layer_lab <- labeller(stratum = c(overstory = "Overstory", understory = "Understory"),
                      pool = c(shared = "Shared pool", restricted = "Restricted pool"))

lab_r <- function(r, r2) sprintf("r = %.2f, R² = %.2f", r, r2)
sft_lab <- sft %>% filter(relationship %in% c("CTI ~ MAT (past)", "CMI ~ MAP (past)"))
pA <- ggplot(plots, aes(MAT_past, CTI_past)) +
  geom_point(colour = col_cti, alpha = 0.5, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.7, se = TRUE, fill = "grey80") +
  geom_text(data = filter(sft_lab, relationship == "CTI ~ MAT (past)"),
            aes(x = -Inf, y = Inf, label = lab_r(pearson_r, r2)), hjust = -0.05, vjust = 1.3, size = 3.6, inherit.aes = FALSE) +
  facet_wrap(~ stratum, scales = "free", labeller = layer_lab) +
  labs(x = "Past MAT (°C)", y = "Past CTI (°C)") + theme_mora()
pB <- ggplot(plots, aes(MAP_past, CMI_past)) +
  geom_point(colour = col_cmi, alpha = 0.5, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.7, se = TRUE, fill = "grey80") +
  geom_text(data = filter(sft_lab, relationship == "CMI ~ MAP (past)"),
            aes(x = -Inf, y = Inf, label = lab_r(pearson_r, r2)), hjust = -0.05, vjust = 1.3, size = 3.6, inherit.aes = FALSE) +
  facet_wrap(~ stratum, scales = "free", labeller = layer_lab) +
  labs(x = "Past MAP (mm)", y = "Past CMI (mm)") + theme_mora()
pC <- ggplot(sp_aff, aes(MAT_mean, MAP_mean)) +
  geom_point(aes(colour = MAT_mean), alpha = 0.8, size = 1.8) +
  scale_colour_gradientn(colours = viridisLite::plasma(256, end = 0.80), limits = c(0, 12.5),
                         oob = scales::squish, guide = "none") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.7, se = TRUE, fill = "grey80") +
  geom_text(data = affc, aes(x = -Inf, y = Inf, label = sprintf("r = %.2f (n = %d)", pearson_r, n)),
            hjust = -0.05, vjust = 1.3, size = 3.6, inherit.aes = FALSE) +
  facet_grid(pool ~ stratum, scales = "free", labeller = layer_lab) +
  labs(x = "Thermal affinity (MAT, °C)", y = "Moisture affinity (MAP, mm)") + theme_mora()

# Plot-level CTI vs CMI: past indices and rates of change
lab_rr <- function(r) sprintf("r = %.2f", r)
pD <- ggplot(plots, aes(CTI_past, CMI_past)) +
  geom_point(colour = "grey40", alpha = 0.5, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.7, se = TRUE, fill = "grey80") +
  geom_text(data = filter(sft, relationship == "CMI ~ CTI (past, plots)"),
            aes(x = Inf, y = Inf, label = lab_rr(pearson_r)), hjust = 1.05, vjust = 1.3, size = 3.6, inherit.aes = FALSE) +
  facet_wrap(~ stratum, scales = "free", labeller = layer_lab) +
  labs(x = "Past CTI (\u00b0C)", y = "Past CMI (mm)") + theme_mora()
pE <- ggplot(plots, aes(rate_CTI, rate_CMI)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(colour = "grey40", alpha = 0.5, size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.7, se = TRUE, fill = "grey80") +
  geom_text(data = filter(sft, relationship == "CMI rate ~ CTI rate (plots)"),
            aes(x = Inf, y = Inf, label = lab_rr(pearson_r)), hjust = 1.05, vjust = 1.3, size = 3.6, inherit.aes = FALSE) +
  facet_wrap(~ stratum, scales = "free", labeller = layer_lab) +
  labs(x = "CTI rate (\u00b0C decade\u207b\u00b9)", y = "CMI rate (mm decade\u207b\u00b9)") + theme_mora()

# Each diagnostic is saved as its own figure, with its own panel letters, so
# that it can be placed next to the section of the Supporting Information it
# supports:
#   - space-for-time relationships (past CTI ~ MAT, CMI ~ MAP), by layer
#   - correlation between species' thermal and moisture affinities, by pool
#   - plot-level covariation of CTI and CMI (past values and rates of change)
figs <- list(
  space_for_time_relationships     = list(plot = (pA / pB) + plot_annotation(tag_levels = "A"),
                                          width = 9, height = 8),
  affinity_correlation             = list(plot = pC, width = 9, height = 7),
  cti_cmi_covariation              = list(plot = (pD / pE) + plot_annotation(tag_levels = "A"),
                                          width = 9, height = 8))
fig_files <- imap_chr(figs, function(x, nm) {
  f <- file.path(fig_dir, paste0(nm, ".png"))
  ggsave(f, x$plot, width = x$width, height = x$height, dpi = 300, bg = "white")
  cat("Figure written to", f, "\n")
  f
})

# Open the figures (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
walk(fig_files, open_fig)