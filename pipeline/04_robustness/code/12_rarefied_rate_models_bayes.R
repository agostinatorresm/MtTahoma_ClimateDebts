# =============================================================================
# 12_rarefied_rate_models_bayes.R
#
# Pipeline: pipeline/04_robustness/code/ (after 09_rarefied_thermophilization.R)
# Author: Agostina Torres
#
# Mean rate of the understory community indices (CTI, CMI) with the past data
# reduced to the present sampling effort, estimated with the same Bayesian
# model as the main results (07_rate_models_bayes.R, model A: rate ~ 1), so that
# observed and reduced rates are reported with comparable credible intervals.
#
#   - The mean-rate model is fitted to the per-plot rates of each of the 100
#     cover-weighted reductions (09), and the posterior draws of the 100 fits
#     are pooled, so that the credible interval includes both the uncertainty
#     of the mean rate within a reduction and the variation among reductions
#     (as for the pooled HMSC occurrence model).
#   - The same model is fitted to the observed data with the same point
#     affinities (09), so that observed and reduced estimates differ only in
#     the past data. The main results (07) use the median of the indices over
#     resamples of the species' climatic affinities; the observed point-
#     affinity estimate is compared with them as a check.
#   - Model as in 07: Gaussian likelihood, response standardized before
#     fitting, priors normal(0, 1) for the intercept and exponential(1) for the
#     residual SD; draws back-transformed to the original units. Compiled once
#     and refitted with update(newdata = ...).
#   - Support: probability of direction f = max(P(> 0), P(< 0)).
#
# Indices: cover-weighted (main) and occurrence-based.
#
# Inputs : pipeline/04_robustness/outputs/rarefaction_thermo/rarefaction_thermo_plot_rates_by_rep.csv (09)
#          pipeline/02_analysis/01_thermophilization/03_integration_over_under/outputs/rate_models_bayes_coefficients.csv (07; check)
# Outputs: pipeline/04_robustness/outputs/rarefaction_thermo/rarefied_rate_models_bayes.csv
#            (index x weight x data: mean, 95/90/50% credible intervals, f, R-hat range)
# =============================================================================

library(tidyverse)
library(here)
library(brms)

SEED      <- 2026
CHAINS    <- 4
ITER      <- 1000           # per chain, half warm-up: 2000 draws per reduction
CORES     <- min(CHAINS, parallel::detectCores())
weighting <- "cover"        # removal used in the main analyses

out_dir <- here("pipeline", "04_robustness", "outputs", "rarefaction_thermo")
rates <- read_csv(file.path(out_dir, "rarefaction_thermo_plot_rates_by_rep.csv"), show_col_types = FALSE) %>%
  filter(weighting %in% c("observed", !!weighting), !is.na(rate))

# -----------------------------------------------------------------------------
# 1. Model (as in 07, model A)
# -----------------------------------------------------------------------------
template <- NULL
fit_mean <- function(y) {
  y_mu <- mean(y); y_sd <- sd(y)
  d <- tibble(y = (y - y_mu) / y_sd)
  fit <- if (is.null(template)) {
    template <<- brm(bf(y ~ 1), data = d, family = gaussian(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     chains = CHAINS, iter = ITER, cores = CORES, seed = SEED, refresh = 0, silent = 2)
  } else {
    update(template, newdata = d, chains = CHAINS, iter = ITER, cores = CORES,
           seed = SEED, refresh = 0, silent = 2)
  }
  dr <- as_draws_df(fit)
  list(draws = dr$b_Intercept * y_sd + y_mu,
       rhat  = posterior::rhat(posterior::extract_variable_matrix(fit, "b_Intercept")))
}

summarise_draws_mean <- function(x) tibble(
  mean = mean(x),
  lower_95 = quantile(x, 0.025), upper_95 = quantile(x, 0.975),
  lower_90 = quantile(x, 0.05),  upper_90 = quantile(x, 0.95),
  lower_50 = quantile(x, 0.25),  upper_50 = quantile(x, 0.75),
  f = max(mean(x > 0), mean(x < 0)))

# -----------------------------------------------------------------------------
# 2. Fits: observed data and each reduction, per index and weight
# -----------------------------------------------------------------------------
groups <- rates %>% distinct(index, weight)
res <- pmap_dfr(groups, function(index, weight) {
  g <- rates %>% filter(index == !!index, weight == !!weight)
  cat("\n---", index, weight, "---\n")
  obs <- fit_mean(g$rate[g$weighting == "observed"])
  reps <- sort(unique(g$rep[g$weighting == weighting]))
  red <- map(reps, function(k) fit_mean(g$rate[g$weighting == weighting & g$rep == k]))
  red_draws <- unlist(map(red, "draws"))
  red_rhat  <- map_dbl(red, "rhat")
  cat(sprintf("Reductions: %d | R-hat max %.3f\n", length(reps), max(red_rhat)))
  bind_rows(
    summarise_draws_mean(obs$draws) %>% mutate(data = "observed", n_fits = 1, rhat_max = obs$rhat),
    summarise_draws_mean(red_draws) %>% mutate(data = paste0("reduced_", weighting, "_pooled"),
                                               n_fits = length(reps), rhat_max = max(red_rhat))
  ) %>% mutate(index = index, weight = weight, .before = 1)
})

write_csv(res %>% mutate(across(where(is.double), ~ signif(.x, 4))),
          file.path(out_dir, "rarefied_rate_models_bayes.csv"))
cat("\nMean rate (CTI: degC per decade; CMI: mm per decade), observed vs reduced (pooled):\n")
print(res %>% select(index, weight, data, mean, lower_95, upper_95, lower_90, upper_90, f, rhat_max) %>%
        mutate(across(where(is.double), ~ signif(.x, 3))), n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 3. Check: observed point-affinity estimate vs main results (07)
# -----------------------------------------------------------------------------
f07 <- here("pipeline", "02_analysis", "01_thermophilization", "03_integration_over_under",
            "outputs", "rate_models_bayes_coefficients.csv")
if (file.exists(f07)) {
  main <- read_csv(f07, show_col_types = FALSE) %>%
    filter(str_detect(model, "^understory (CTI|CMI) A mean rate$"), term == "Intercept") %>%
    transmute(index = str_extract(model, "CTI|CMI"), mean_07 = mean, lower_95_07 = lower_95, upper_95_07 = upper_95)
  cat("\nObserved data, point affinities (this script, cover-weighted) vs median over affinity resamples (07):\n")
  print(res %>% filter(data == "observed", weight == "cover") %>%
          select(index, mean, lower_95, upper_95) %>% left_join(main, by = "index") %>%
          mutate(across(where(is.double), ~ signif(.x, 3))), width = Inf)
}







library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

res <- read_csv(here("pipeline", "04_robustness", "outputs", "rarefaction_thermo",
                     "rarefied_rate_models_bayes.csv"), show_col_types = FALSE)

fig_dir <- here("pipeline", "04_robustness", "figures")
col_index <- c(CTI = "#F26B5B", CMI = "#1A759F")     # temperature, moisture (as in the main figures)
y_lab <- c(CTI = "Thermal affinity shift\n(\u00b0C per decade)", CMI = "Moisture affinity shift\n(mm per decade)")
lab_weight <- c(cover = "Cover-weighted\n(main)", occurrence = "Occurrence-\nbased")
lab_data <- c(observed = "Observed past data (main)", reduced = "Reduced past data (pooled)")

fig_d <- res %>%
  mutate(data = factor(if_else(data == "observed", "observed", "reduced"), levels = names(lab_data)),
         weight = factor(weight, levels = names(lab_weight)),
         supported = f >= 0.90)

panel <- function(ix) {
  dodge <- position_dodge(width = 0.5)
  ggplot(filter(fig_d, index == ix), aes(weight, mean, colour = data, alpha = supported)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(ymin = lower_90, ymax = upper_90), linewidth = 0.8, lineend = "round", position = dodge) +
    geom_linerange(aes(ymin = lower_50, ymax = upper_50), linewidth = 2.2, lineend = "round", position = dodge) +
    geom_point(size = 2.8, position = dodge) +
    scale_colour_manual(values = c(observed = col_index[[ix]], reduced = "black"), labels = lab_data, name = NULL) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
    scale_x_discrete(labels = lab_weight) +
    labs(x = NULL, y = y_lab[[ix]]) +
    theme_ipsum(base_size = 13) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.text = element_text(colour = "black", size = 12),
          axis.title.y = element_text(colour = "black", size = 13, hjust = 1),
          legend.position = "bottom", legend.text = element_text(size = 12))
}

fig <- (panel("CTI") | panel("CMI")) +
  plot_annotation(title = "\U0001F33F Understory", tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", margin = margin(b = 0))))
fig_file <- file.path(fig_dir, "figS4_1_rarefied_rate_models_bayes.png")
ggsave(fig_file, fig, width = 10, height = 5.5, dpi = 300, bg = "white")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
