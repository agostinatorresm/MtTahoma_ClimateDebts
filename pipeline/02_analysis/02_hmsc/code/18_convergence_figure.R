# =============================================================================
# 18_convergence_figure.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine
#           after 03_check_convergence.R, or on the local computer after copying the
#           *_conv_parameters.csv tables listed below)
# Author: Agostina Torres
#
# MCMC convergence of the final HMSC models (Supporting Information S5):
# distribution of the potential scale reduction factor (PSRF) and of the
# effective sample size (ESS) per parameter, by parameter group, for
#   overstory  occurrence and relative basal area (medium run)
#   understory occurrence (10 fits to reduced past datasets, pooled for the
#              results) and paired relative cover (long runs)
# Parameter groups: species responses (Beta), trait effects (Gamma),
# phylogenetic signal (Rho) and plot-level residual associations (Omega).
# Constant parameters (PSRF undefined) are excluded.
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/overstory_strict_medium_conv_parameters.csv (03)
#          pipeline/02_analysis/02_hmsc/outputs/understory_strict_long_run_oct09_conv_parameters.csv (03)
# Outputs: pipeline/02_analysis/02_hmsc/outputs/convergence_final_models_summary.csv
#          pipeline/02_analysis/02_hmsc/figures/m1_final/ms1__figS5_1_convergence_m1_final.png
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
res_dir    <- file.path(module_dir, "outputs")
fig_dir    <- file.path(module_dir, "figures", "m1_final")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
psrf_target <- 1.1

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line  = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 11),
          axis.title.y = element_text(colour = "black", size = 11),
          strip.text   = element_text(colour = "black", face = "bold", size = 10),
          plot.margin  = margin(12, 10, 10, 10))
}

# Final models (model_id as written by 03)
models <- tribble(
  ~file,                                          ~pattern,                                   ~panel,
  "overstory_strict_medium_conv_parameters.csv",  "^m1_time_elev_nf2_nsp_nos_pa$",            "Overstory:\noccurrence",
  "overstory_strict_medium_conv_parameters.csv",  "^m1_time_elev_nf2_nsp_nos_rel$",           "Overstory: relative\nbasal area",
  "understory_strict_long_run_oct09_conv_parameters.csv", "^m1_time_elev_nsp_nos_red\\d+_pa$", "Understory: occurrence\n(10 reduced fits)",
  "understory_strict_long_run_oct09_conv_parameters.csv", "^m1_time_elev_nsp_nos_rel$",       "Understory:\nrelative cover"
)
group_labels <- c(Beta = "β", Gamma = "Γ", Rho = "ρ", Omega_plot = "Ω")

conv <- pmap_dfr(models, function(file, pattern, panel) {
  f <- file.path(res_dir, file)
  if (!file.exists(f)) stop("Run 03_check_convergence.R first (or copy the file): ", f)
  read_csv(f, show_col_types = FALSE) %>%
    filter(str_detect(model_id, pattern), group %in% names(group_labels), !is.na(psrf)) %>%
    mutate(panel = panel)
}) %>%
  mutate(panel = factor(panel, levels = models$panel),
         group = factor(group_labels[group], levels = group_labels))

summary_tbl <- conv %>%
  group_by(panel, group) %>%
  summarise(n_parameters = n(), psrf_median = round(median(psrf), 3), psrf_max = round(max(psrf), 3),
            prop_psrf_high = round(mean(psrf > psrf_target), 3),
            ess_median = round(median(ess)), ess_min = round(min(ess)), .groups = "drop")
print(summary_tbl, n = Inf, width = Inf)
write_csv(summary_tbl, file.path(res_dir, "convergence_final_models_summary.csv"))

# Figure: PSRF (top) and ESS (bottom) per parameter, by group and model
long <- conv %>%
  select(panel, group, PSRF = psrf, ESS = ess) %>%
  pivot_longer(c(PSRF, ESS), names_to = "statistic") %>%
  mutate(statistic = factor(statistic, levels = c("PSRF", "ESS")))
ref_lines <- tibble(statistic = factor(c("PSRF", "ESS"), levels = c("PSRF", "ESS")),
                    value = c(psrf_target, 1000))

fig <- ggplot(long, aes(group, value)) +
  geom_hline(data = ref_lines, aes(yintercept = value), linetype = "dashed", colour = "grey50") +
  geom_jitter(width = 0.2, height = 0, size = 0.6, alpha = 0.35, colour = "grey30") +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5, colour = "#9A3F12") +
  facet_grid(statistic ~ panel, scales = "free_y", switch = "y") +
  labs(x = "Parameter group", y = NULL) +
  theme_mora() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))
fig_file <- file.path(fig_dir, "ms1__figS5_1_convergence_m1_final.png")
ggsave(fig_file, fig, width = 11, height = 6, dpi = 300, bg = "white")
cat("Figure written to", fig_file, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
