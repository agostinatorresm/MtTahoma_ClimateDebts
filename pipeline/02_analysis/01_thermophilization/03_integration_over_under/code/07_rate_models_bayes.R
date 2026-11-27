# ------------------------------------------------------------------------------
# 07_rate_models_bayes.R
# Module: pipeline/02_analysis/01_thermophilization/03_integration_over_under
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: Bayesian versions of the rate models of 06_rate_models.R, so that
# the climatic-affinity shifts are reported with posterior means, credible
# intervals (95%, 90% and 50%) and the posterior probability of the direction of each effect
# ("support"), as for the HMSC results.
#
# MODELS (Gaussian likelihood; predictors centred and scaled, so slopes are per
# 1 SD; elevation scaled once across plots, as in 06):
#   A.  Mean rate:                    rate ~ 1                  (per stratum)
#   A2. Mean climatic lag:            lag ~ 1                   (per stratum)
#   B.  Rate vs elevation:            rate ~ elevation          (per stratum)
#   C.  Rate vs measured climate-change rate (MAT trend for CTI, MAP trend
#       for CMI):                     rate ~ climate rate       (per stratum)
#   C2. Lag vs climate-change rate:   lag ~ climate rate        (per stratum)
#   D.  Lag vs elevation:             lag ~ elevation           (per stratum)
#   F.  Understory vs overstory canopy: understory rate (CTI, CMI) and thermal
#       lag ~ change in total overstory basal area (present − past), and
#       understory rate ~ past total basal area. A closing canopy creates a
#       cooler, moister understory microclimate that could buffer the
#       understory against regional warming (expected: lower thermophilization
#       and larger lag where basal area increased).
#   G.  Community tails: mean rate of the cold/dry tail (q10) and of the
#       warm/wet tail (q90) of the abundance-weighted affinity distribution,
#       and mean of the per-plot difference q10 − q90 (Bayesian equivalent of
#       the paired t-test of 06_community_affinity_edges.R). A positive
#       difference means the cold/dry tail rose more than the warm/wet tail
#       (loss of cold- or dry-affiliated species); a negative difference, a
#       gain of warm- or wet-affiliated species.
#   E.  Difference between layers (understory − overstory) in the mean rate (A)
#       and in the elevational slope (B), computed from the posterior draws of
#       the per-layer models. This replaces a joint model with a plot random
#       effect and assumes that the rates of the two layers are independent
#       within plots; their within-plot correlation is reported as a check.
#
# Fitting: brms (Stan). Only two model structures are used (y ~ 1 and y ~ x),
# so each is compiled once and refitted to the other responses and predictors
# with update(newdata = ...), which avoids recompiling.
# Each response is standardized (mean 0, SD 1) before fitting, so the same
# weakly informative priors apply to all responses: normal(0, 1) for the
# intercept and slopes, exponential(1) for the residual SD.
# Coefficients are back-transformed to the original units (°C or mm per decade)
# before summarising; the probability of direction is unaffected by the
# scaling. 4 chains of 2000 iterations (1000 warm-up). Convergence: R-hat < 1.01
# and bulk effective sample size > 400 for every coefficient (reported in the
# output table).
#
# Support: for each coefficient, the posterior probability that it is positive
# (p_positive) and the probability of direction (pd = max(p_positive,
# 1 - p_positive)). pd = 0.95 corresponds roughly to a one-sided frequentist
# p = 0.05, i.e. a two-sided p = 0.10.
#
# Between-plot variability: for the mean models, the posterior predictive
# distribution of the rate of a new plot (Normal(intercept, sigma) over the
# posterior draws) gives 90% and 50% prediction intervals.
#
# Residual spatial autocorrelation: Moran's I on the posterior mean residuals
# (inverse-distance weights on projected coordinates, UTM 10N, km).
#
# Inputs : <stratum>_thermo_xero_rates.csv and <stratum>_thermo_xero_expected.csv
#          (same as 06); pipeline/01_data_prep/01_cleaning/outputs/overstory_long.csv
#          (total basal area per plot and time, for F);
#          <stratum>_affinity_edge_rates.csv (06_community_affinity_edges.R, for G)
# Outputs: outputs/rate_models_bayes_coefficients.csv
#          outputs/rate_models_bayes_strata_differences.csv
#          outputs/rate_models_bayes_moran.csv
#          outputs/rate_models_bayes_predictive.csv (posterior predictive
#            intervals of a new plot for the mean models: between-plot
#            variability estimated by the model)
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(brms)
library(ape)

SEED   <- 2026
CHAINS <- 4
ITER   <- 2000
CORES  <- min(CHAINS, parallel::detectCores())

thermo_dir      <- function(...) here("pipeline", "02_analysis", "01_thermophilization", ...)
this_module_dir <- function(...) thermo_dir("03_integration_over_under", ...)
dir.create(this_module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)

STRATA <- tibble(
  stratum = c("overstory", "understory"),
  dir     = c("01_overstory", "02_understory"),
  weight  = c("basal_area_spp_perc", "cover")
)

# 1. Load per-plot rates and expected-rate tables (as in 06) -------------------
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

elev_ref <- plots %>% distinct(plot_id, elevation)
plots <- plots %>%
  mutate(elev_sc = (elevation - mean(elev_ref$elevation, na.rm = TRUE)) / sd(elev_ref$elevation, na.rm = TRUE))

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
# Posterior summary of the fixed effects in the original units: mean, 95%, 90%
# and 50% credible intervals (the 90% and 50% intervals are used in the
# figures, as for the HMSC results), probability of a positive effect, probability of direction
# and convergence. y_mu and y_sd are the mean and SD used to standardize y.
post_table <- function(fit, label, y_mu, y_sd, predictor_name = NULL) {
  draws <- as_draws_df(fit) %>% as_tibble() %>% select(starts_with("b_")) %>%
    mutate(across(everything(), ~ .x * y_sd)) %>%
    mutate(b_Intercept = b_Intercept + y_mu)
  conv  <- posterior::summarise_draws(as_draws_df(fit) %>% select(starts_with("b_")),
                                      "rhat", "ess_bulk")
  out <- tibble(term       = names(draws),
                mean       = map_dbl(draws, mean),
                lower_95   = map_dbl(draws, quantile, 0.025),
                upper_95   = map_dbl(draws, quantile, 0.975),
                lower_90   = map_dbl(draws, quantile, 0.05),
                upper_90   = map_dbl(draws, quantile, 0.95),
                lower_50   = map_dbl(draws, quantile, 0.25),
                upper_50   = map_dbl(draws, quantile, 0.75),
                p_positive = map_dbl(draws, ~ mean(.x > 0))) %>%
    mutate(pd = pmax(p_positive, 1 - p_positive),
           direction = if_else(p_positive >= 0.5, "positive", "negative")) %>%
    left_join(conv, by = c("term" = "variable")) %>%
    mutate(term = str_remove(term, "^b_"))
  if (!is.null(predictor_name)) out <- mutate(out, term = str_replace(term, "^x$", predictor_name))
  names(draws) <- str_remove(names(draws), "^b_")
  list(table = mutate(out, model = label, .before = 1), draws = draws)
}

# Posterior predictive distribution of a new plot for the mean models (y ~ 1):
# for each posterior draw, a value from Normal(intercept, sigma), in the original
# units. Its quantiles describe the variability among plots estimated by the
# model (prediction interval), unlike the credible interval of the mean.
predictive_table <- function(fit, label, y_mu, y_sd) {
  dr <- as_draws_df(fit)
  set.seed(SEED)
  pred <- (dr$b_Intercept + rnorm(nrow(dr), 0, dr$sigma)) * y_sd + y_mu
  tibble(model = label, sigma_mean = mean(dr$sigma) * y_sd,
         pred_lower_90 = quantile(pred, 0.05), pred_upper_90 = quantile(pred, 0.95),
         pred_lower_50 = quantile(pred, 0.25), pred_upper_50 = quantile(pred, 0.75),
         p_new_plot_positive = mean(pred > 0))
}

# Moran's I on the posterior mean residuals
moran_resid <- function(fit, d) {
  r <- residuals(fit)[, "Estimate"]
  d <- d %>% mutate(r = r) %>% filter(!is.na(x_km), !is.na(y_km))
  w <- 1 / as.matrix(dist(cbind(d$x_km, d$y_km)))
  diag(w) <- 0
  w[!is.finite(w)] <- 0     # co-located rows (both layers of one plot) get zero weight
  mi <- Moran.I(d$r, w)
  tibble(moran_I = mi$observed, moran_p = mi$p.value)
}

# Templates of the two model structures, compiled once on the first data
# used and refitted afterwards with update(newdata = ...)
templates <- list()
priors_fixed <- c(prior(normal(0, 1), class = "Intercept"),
                  prior(exponential(1), class = "sigma"))
fit_brm <- function(type, d) {
  f <- switch(type,
              mean  = bf(y ~ 1),
              slope = bf(y ~ x))
  if (is.null(templates[[type]])) {
    pr <- priors_fixed
    if (type != "mean")  pr <- c(pr, prior(normal(0, 1), class = "b"))
    templates[[type]] <<- brm(f, data = d, family = gaussian(), prior = pr,
                              chains = CHAINS, iter = ITER,
                              cores = CORES, seed = SEED, refresh = 0, silent = 2)
    templates[[type]]
  } else {
    update(templates[[type]], newdata = d, chains = CHAINS, iter = ITER,
           cores = CORES, seed = SEED, refresh = 0, silent = 2)
  }
}

all_post  <- list()
all_draws <- list()
all_pred  <- list()
all_moran <- list()

# type: "mean" or "slope"; y: response column; x: predictor column
run_model <- function(label, type, data, y, x = NULL) {
  cat("\n---", label, "---\n")
  d <- data %>% mutate(y = .data[[y]])
  if (!is.null(x)) d <- mutate(d, x = .data[[x]])
  d <- d %>% filter(!is.na(y), if (!is.null(x)) !is.na(x) else TRUE)
  y_mu <- mean(d$y); y_sd <- sd(d$y)
  d <- d %>% mutate(y = (y - y_mu) / y_sd)
  fit <- fit_brm(type, d)
  res <- post_table(fit, label, y_mu, y_sd, predictor_name = x)
  pt  <- res$table
  print(pt %>% select(term, mean, lower_95, upper_95, pd, direction, rhat, ess_bulk), width = Inf)
  mo <- moran_resid(fit, d)
  cat("Moran's I:", signif(mo$moran_I, 3), "| p-value:", signif(mo$moran_p, 3), "\n")
  all_post[[label]]  <<- pt
  all_draws[[label]] <<- res$draws
  # Mean models: posterior predictive distribution of the rate of a new plot
  # (between-plot variability estimated by the model)
  if (type == "mean") all_pred[[label]] <<- predictive_table(fit, label, y_mu, y_sd)
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
    run_model(paste(st, v, "A mean rate"),           "mean",  d, rate)
    run_model(paste(st, v, "A2 mean lag"),           "mean",  d, lag)
    run_model(paste(st, v, "B elevation"),           "slope", d, rate, "elev_sc")
    run_model(paste(st, v, "C climate change rate"), "slope", d, rate, "clim_sc")
    run_model(paste(st, v, "C2 lag ~ climate rate"), "slope", d, lag,  "clim_sc")
    run_model(paste(st, v, "D lag ~ elevation"),     "slope", d, lag,  "elev_sc")
  }
}

# 3b. Understory shifts vs overstory canopy (F) --------------------------------
# Total overstory basal area (% of plot area) per plot and time
ba <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", "overstory_long.csv"),
               show_col_types = FALSE) %>%
  group_by(plot_id, period) %>%
  summarise(ba = sum(basal_area_spp_perc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = ba, names_prefix = "ba_", values_fill = 0) %>%
  mutate(d_ba = ba_present - ba_past)
cat("\nOverstory basal area (% of plot area): past mean", round(mean(ba$ba_past), 2),
    "| present mean", round(mean(ba$ba_present), 2),
    "| plots with an increase:", round(100 * mean(ba$d_ba > 0)), "%\n")

d_u <- plots %>%
  filter(stratum == "understory") %>%
  left_join(ba, by = "plot_id") %>%
  mutate(d_ba_sc = as.numeric(scale(d_ba)), ba_past_sc = as.numeric(scale(ba_past)))

cat("\n==============================\nUNDERSTORY VS OVERSTORY CANOPY\n==============================\n")
run_model("understory CTI F rate ~ change in basal area", "slope", d_u, "rate_CTI", "d_ba_sc")
run_model("understory CTI F lag ~ change in basal area",  "slope", d_u, "lag_CTI",  "d_ba_sc")
run_model("understory CMI F rate ~ change in basal area", "slope", d_u, "rate_CMI", "d_ba_sc")
run_model("understory CTI F rate ~ past basal area",      "slope", d_u, "rate_CTI", "ba_past_sc")
run_model("understory CMI F rate ~ past basal area",      "slope", d_u, "rate_CMI", "ba_past_sc")

# 3c. Community tails (G) -------------------------------------------------------
cat("\n==============================\nCOMMUNITY EDGES\n==============================\n")
for (i in seq_len(nrow(STRATA))) {
  st <- STRATA$stratum[i]
  er <- read_csv(thermo_dir(STRATA$dir[i], "outputs", paste0(st, "_affinity_edge_rates.csv")),
                 show_col_types = FALSE) %>%
    filter(edge %in% c("q10", "q90")) %>%
    select(plot_id, index, edge, rate) %>%
    pivot_wider(names_from = edge, values_from = rate) %>%
    mutate(q10_minus_q90 = q10 - q90) %>%
    left_join(xy, by = "plot_id")
  for (v in c("CTI", "CMI")) {
    d_e <- filter(er, index == v)
    lab <- if (v == "CTI") c("cold edge", "warm edge") else c("dry edge", "wet edge")
    run_model(paste(st, v, "G", lab[1], "(q10) rate"), "mean", d_e, "q10")
    run_model(paste(st, v, "G", lab[2], "(q90) rate"), "mean", d_e, "q90")
    run_model(paste(st, v, "G q10 - q90"),             "mean", d_e, "q10_minus_q90")
  }
}

# 4. Differences between layers ------------------------------------------------
# Posterior of understory − overstory for the mean rate (intercept of A) and the
# elevational slope (B). Draws of the two layers come from independent fits, so
# they are paired by draw index.
cat("\n==============================\nDIFFERENCES BETWEEN STRATA\n==============================\n")

w <- plots %>% select(plot_id, stratum, rate_CTI, rate_CMI) %>%
  pivot_wider(names_from = stratum, values_from = c(rate_CTI, rate_CMI))
cat("Within-plot correlation of overstory and understory rates: CTI r =",
    round(cor(w$rate_CTI_overstory, w$rate_CTI_understory, use = "complete.obs"), 2),
    "| CMI r =", round(cor(w$rate_CMI_overstory, w$rate_CMI_understory, use = "complete.obs"), 2), "\n")

diff_summary <- function(x, label) {
  tibble(comparison = label, mean = mean(x),
         lower_95 = quantile(x, 0.025), upper_95 = quantile(x, 0.975),
         p_positive = mean(x > 0)) %>%
    mutate(pd = pmax(p_positive, 1 - p_positive),
           direction = if_else(p_positive >= 0.5, "positive", "negative"))
}

strata_diff <- map_dfr(c("CTI", "CMI"), function(v) {
  a_o <- all_draws[[paste("overstory",  v, "A mean rate")]]$Intercept
  a_u <- all_draws[[paste("understory", v, "A mean rate")]]$Intercept
  b_o <- all_draws[[paste("overstory",  v, "B elevation")]]$x
  b_u <- all_draws[[paste("understory", v, "B elevation")]]$x
  n <- min(length(a_o), length(a_u), length(b_o), length(b_u))
  bind_rows(diff_summary(a_u[1:n] - a_o[1:n], paste(v, "mean rate (understory - overstory)")),
            diff_summary(b_u[1:n] - b_o[1:n], paste(v, "elevational slope (understory - overstory)")))
})
print(strata_diff, width = Inf)
write_csv(strata_diff %>% mutate(across(c(mean, lower_95, upper_95), ~ signif(.x, 4)),
                                 across(c(p_positive, pd), ~ round(.x, 3))),
          this_module_dir("outputs", "rate_models_bayes_strata_differences.csv"))

# 5. Save ----------------------------------------------------------------------
post_out <- bind_rows(all_post) %>%
  mutate(across(c(mean, starts_with("lower_"), starts_with("upper_")), ~ signif(.x, 4)),
         across(c(p_positive, pd), ~ round(.x, 3)))
write_csv(post_out, this_module_dir("outputs", "rate_models_bayes_coefficients.csv"))
write_csv(bind_rows(all_moran), this_module_dir("outputs", "rate_models_bayes_moran.csv"))
pred_out <- bind_rows(all_pred) %>% mutate(across(where(is.numeric), ~ signif(.x, 4)))
write_csv(pred_out, this_module_dir("outputs", "rate_models_bayes_predictive.csv"))
cat("\nPosterior predictive intervals for a new plot (mean models):\n"); print(pred_out, width = Inf)

cat("\nConvergence: max R-hat =", round(max(post_out$rhat, na.rm = TRUE), 3),
    "| min bulk ESS =", round(min(post_out$ess_bulk, na.rm = TRUE)), "\n")
cat("Saved: outputs/rate_models_bayes_coefficients.csv, outputs/rate_models_bayes_strata_differences.csv,",
    "outputs/rate_models_bayes_moran.csv\n")
