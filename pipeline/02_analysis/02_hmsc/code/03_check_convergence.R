# =============================================================================
# 03_check_convergence.R
# Author: Agostina Torres
#
# MCMC convergence diagnostics for HMSC models fitted by 02_fit_hmsc.R.
#
# Diagnostics per model and parameter group:
#   - potential scale reduction factor (PSRF, Gelman-Rubin); target < 1.1
#   - effective sample size (ESS) across chains
#   Parameter groups: Beta (species responses), Gamma (trait effects on
#   responses), Rho (phylogenetic signal), Alpha (spatial scale of the plot
#   random level), Omega (residual species associations, subsampled).
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_<model_id>.rds
#          pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_run_log.csv
# Outputs: pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_conv_*.csv
#            (parameters, summary, psrf_above_target, beta_by_covariate,
#             beta_by_species, alpha_by_chain)
#          pipeline/02_analysis/02_hmsc/figures/<stratum>_<pool>_<run_mode>_conv_*.png
# =============================================================================

library(here)
library(tidyverse)
library(Hmsc)
library(coda)
library(hrbrthemes)
library(viridis)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
stratum  <- "overstory"
pool     <- "strict"
run_mode <- "short"

# Folder of the fitted models inside 02_hmsc/outputs ("" = outputs itself)
fit_subdir <- ""
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) >= 1) stratum    <- cmd_args[1]
if (length(cmd_args) >= 2) run_mode   <- cmd_args[2]
if (length(cmd_args) >= 3) fit_subdir <- cmd_args[3]

n_omega_check <- 200   # Omega elements subsampled for diagnostics
psrf_target   <- 1.1
seed          <- 1

# -----------------------------------------------------------------------------
# Plot style (shared across project figures)
# -----------------------------------------------------------------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
}

# In non-interactive runs (Rscript) printed plots go to a temporary png device,
# which uses the Windows fonts registered above
if (!interactive()) png(tempfile(fileext = ".png"))

theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid       = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x      = element_line(colour = "black", linewidth = 0.6),
          axis.line.y      = element_line(colour = "black", linewidth = 0.6),
          axis.ticks       = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          axis.text        = element_text(colour = "black"),
          axis.title.x     = element_text(colour = "black", size = 11),
          axis.title.y     = element_text(colour = "black", size = 11),
          legend.text      = element_text(colour = "black"),
          legend.title     = element_text(colour = "black"),
          strip.text       = element_text(colour = "black", face = "bold", size = 10),
          plot.title       = element_text(colour = "black", size = 11, face = "bold"),
          plot.subtitle    = element_text(colour = "black"),
          plot.caption     = element_text(colour = "black"),
          plot.margin      = margin(12, 10, 10, 10))
}

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
run_tag    <- paste(stratum, pool, sep = "_")
model_dir  <- file.path(module_dir, "outputs", fit_subdir)
conv_dir   <- file.path(module_dir, "outputs")
# Figures go to figures/ or to a subfolder set with options(mora.fig_subdir = "...")
fig_dir    <- file.path(module_dir, "figures", getOption("mora.fig_subdir", ""))
walk(c(conv_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

# File names always start with <stratum>_<pool>_<run_mode>
file_prefix <- paste0(run_tag, "_", run_mode, "_")
# Diagnostics of fits kept in a subfolder carry the folder name, so they do not
# overwrite those of the fits in outputs/
out_prefix  <- if (nzchar(fit_subdir)) paste0(file_prefix, fit_subdir, "_") else file_prefix
conv_file   <- function(name) file.path(conv_dir, paste0(out_prefix, "conv_", name))
conv_fig    <- function(name) file.path(fig_dir, paste0(out_prefix, "conv_", name))

run_log <- read_csv(file.path(model_dir, paste0(file_prefix, "run_log.csv")),
                    show_col_types = FALSE) %>%
  group_by(model_id) %>%
  slice_tail(n = 1) %>%          # latest fit of each model if refitted
  ungroup()
run_log

model_files <- file.path(model_dir, run_log$file)
names(model_files) <- run_log$model_id
model_files <- model_files[file.exists(model_files)]
if (length(model_files) == 0) stop("No fitted models found for prefix ", file_prefix)

# -----------------------------------------------------------------------------
# 1. Diagnostics
# -----------------------------------------------------------------------------
diag_group <- function(mcmc_list, group) {
  if (is.null(mcmc_list)) return(NULL)
  psrf <- gelman.diag(mcmc_list, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]
  ess  <- effectiveSize(mcmc_list)
  tibble(group = group, parameter = names(ess), psrf = unname(psrf), ess = unname(ess))
}

subsample_mcmc <- function(mcmc_list, n) {
  n_par <- ncol(mcmc_list[[1]])
  if (n_par <= n) return(mcmc_list)
  set.seed(seed)
  keep <- sort(sample(n_par, n))
  as.mcmc.list(lapply(mcmc_list, function(x) x[, keep, drop = FALSE]))
}

# Posterior of the spatial scale (Alpha, km) of each latent factor of the plot
# level, per chain. Alpha = 0 means no spatial structure in that factor;
# chains that settle on different Alpha values (some at 0, some > 0) indicate
# that the spatial scale is weakly identified.
alpha_by_chain <- function(alpha_mcmc) {
  if (is.null(alpha_mcmc)) return(NULL)
  imap_dfr(alpha_mcmc, function(ch, k) {
    as_tibble(as.matrix(ch)) %>%
      pivot_longer(everything(), names_to = "factor", values_to = "alpha_km") %>%
      group_by(factor) %>%
      summarise(chain         = k,
                alpha_median  = median(alpha_km),
                alpha_q025    = quantile(alpha_km, 0.025),
                alpha_q975    = quantile(alpha_km, 0.975),
                prop_zero     = mean(alpha_km == 0),
                n_distinct    = n_distinct(alpha_km),
                .groups = "drop")
  })
}

alpha_list <- list()

convergence <- imap_dfr(model_files, function(f, id) {
  m_fit <- readRDS(f)
  post  <- convertToCodaObject(m_fit, spNamesNumbers = c(TRUE, FALSE),
                               covNamesNumbers = c(TRUE, FALSE))
  rl_names <- names(m_fit$ranLevels)
  i_plot   <- which(rl_names == "plot")

  if (length(i_plot) == 1 && !is.null(post$Alpha[[i_plot]])) {
    alpha_list[[id]] <<- alpha_by_chain(post$Alpha[[i_plot]]) %>%
      mutate(model_id = id, .before = 1)
  }

  bind_rows(
    diag_group(post$Beta,  "Beta"),
    diag_group(post$Gamma, "Gamma"),
    if (!is.null(post$Rho)) diag_group(post$Rho, "Rho"),
    if (length(i_plot) == 1 && !is.null(post$Alpha[[i_plot]]))
      diag_group(post$Alpha[[i_plot]], "Alpha_plot"),
    map_dfr(seq_along(rl_names), function(k) {
      diag_group(subsample_mcmc(post$Omega[[k]], n_omega_check),
                 paste0("Omega_", rl_names[k]))
    })
  ) %>%
    mutate(model_id = id, .before = 1)
})

posterior_size <- run_log %>%
  distinct(model_id, samples, n_chains) %>%
  mutate(posterior_size = samples * n_chains) %>%
  select(model_id, posterior_size)

# Parameters that never change across the posterior (e.g. a spatial scale fixed
# at 0 in every sample) have no defined PSRF; they are counted separately
# (n_constant) and excluded from the PSRF and ESS summaries.
safe_stat <- function(x, f) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)

convergence_summary <- convergence %>%
  group_by(model_id, group) %>%
  summarise(n_parameters   = n(),
            n_constant     = sum(is.na(psrf)),
            psrf_median    = round(safe_stat(psrf, median), 3),
            psrf_max       = round(safe_stat(psrf, max), 3),
            prop_psrf_high = round(safe_stat(psrf > psrf_target, mean), 3),
            ess_median     = round(safe_stat(if_else(is.na(psrf), NA_real_, ess), median)),
            ess_min        = round(safe_stat(if_else(is.na(psrf), NA_real_, ess), min)),
            .groups = "drop") %>%
  left_join(posterior_size, by = "model_id")
convergence_summary

# Compact comparison across models (one row per model): proportion of
# parameters above the PSRF target and median ESS per parameter group, and run
# time. Used to compare fixed-effect and random-level structures fitted with
# the same MCMC settings.
model_comparison <- convergence_summary %>%
  select(model_id, group, prop_psrf_high, ess_median) %>%
  pivot_wider(names_from = group, values_from = c(prop_psrf_high, ess_median),
              names_glue = "{group}_{.value}") %>%
  left_join(select(run_log, model_id, elapsed_min), by = "model_id") %>%
  select(model_id, elapsed_min, ends_with("prop_psrf_high"), ends_with("ess_median"))
model_comparison

# Parameters above the PSRF target (to identify problematic species or terms)
psrf_high <- convergence %>%
  filter(psrf > psrf_target) %>%
  arrange(desc(psrf))
psrf_high

# Beta convergence by covariate and by species (to separate term-specific from
# species-specific problems). Parameter names (names without numbers):
# B[<covariate>, <species>]
beta_conv <- convergence %>%
  filter(group == "Beta") %>%
  mutate(covariate = str_match(parameter, "^B\\[(.+), ([^,]+)\\]$")[, 2],
         species   = str_match(parameter, "^B\\[(.+), ([^,]+)\\]$")[, 3])
stopifnot("Beta parameter names not parsed" = !anyNA(beta_conv$species))

beta_by_covariate <- beta_conv %>%
  group_by(model_id, covariate) %>%
  summarise(psrf_median = round(median(psrf), 3), psrf_max = round(max(psrf), 3),
            n_psrf_high = sum(psrf > psrf_target), ess_min = round(min(ess)),
            .groups = "drop")
beta_by_covariate

beta_by_species <- beta_conv %>%
  group_by(model_id, species) %>%
  summarise(psrf_max = round(max(psrf), 3), n_psrf_high = sum(psrf > psrf_target),
            ess_min = round(min(ess)), .groups = "drop") %>%
  arrange(model_id, desc(psrf_max))
beta_by_species

alpha_summary <- bind_rows(alpha_list)
alpha_summary

# Spatial scale per latent factor, summarised over chains: mean proportion of
# samples at 0 and range of the chain medians (a wide range = chains disagree)
alpha_overview <- alpha_summary %>%
  group_by(model_id, factor) %>%
  summarise(prop_zero_mean   = round(mean(prop_zero), 3),
            chain_median_min = round(min(alpha_median), 3),
            chain_median_max = round(max(alpha_median), 3),
            .groups = "drop")
alpha_overview

write_csv(convergence, conv_file("parameters.csv"))
write_csv(convergence_summary, conv_file("summary.csv"))
write_csv(psrf_high, conv_file("psrf_above_target.csv"))
write_csv(beta_by_covariate, conv_file("beta_by_covariate.csv"))
write_csv(beta_by_species, conv_file("beta_by_species.csv"))
write_csv(model_comparison, conv_file("model_comparison.csv"))
if (nrow(alpha_summary) > 0) {
  write_csv(alpha_summary, conv_file("alpha_by_chain.csv"))
  write_csv(alpha_overview, conv_file("alpha_overview.csv"))
}

# -----------------------------------------------------------------------------
# 2. Figures
# -----------------------------------------------------------------------------
conv_plot_data <- convergence %>%
  mutate(group = factor(group, levels = unique(group)))

p_psrf <- ggplot(conv_plot_data, aes(x = group, y = psrf, fill = group)) +
  geom_hline(yintercept = psrf_target, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~ model_id, scales = "free_x") +
  scale_fill_viridis_d(option = "viridis", end = 0.9) +
  labs(x = NULL, y = "PSRF",
       title = paste0("Potential scale reduction factor (", run_mode, " run)")) +
  theme_mora() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_psrf

p_ess <- ggplot(conv_plot_data, aes(x = group, y = ess, fill = group)) +
  geom_hline(data = posterior_size, aes(yintercept = posterior_size),
             linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~ model_id, scales = "free_x") +
  scale_fill_viridis_d(option = "viridis", end = 0.9) +
  labs(x = NULL, y = "Effective sample size",
       title = paste0("Effective sample size (", run_mode, " run; dashed = posterior size)")) +
  theme_mora() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_ess

ggsave(conv_fig("psrf.png"), p_psrf,
       device = "png", dpi = 300, width = 7, height = 4, scale = 1.8, bg = "transparent")
ggsave(conv_fig("ess.png"), p_ess,
       device = "png", dpi = 300, width = 7, height = 4, scale = 1.8, bg = "transparent")
