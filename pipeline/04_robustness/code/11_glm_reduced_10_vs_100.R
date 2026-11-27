# =============================================================================
# 11_glm_reduced_10_vs_100.R
#
# Pipeline: pipeline/04_robustness/code/ (run after 05 and 06)
# Author: Agostina Torres
#
# Number of reduced past datasets: the understory occurrence HMSC model was
# fitted to 10 of the 100 past datasets reduced to the present sampling effort
# (05; cover-weighted removal), and its posterior pooled across them. To check
# that 10 datasets give the same species changes as 100, the single-species
# GLMs of 06 (probit, presence ~ time * elevation) are fitted to each of the
# 100 cover-weighted reduced datasets, and the results averaged over (i) all
# 100 and (ii) the 10 datasets used in HMSC are compared. The model is the same
# in both sets, so the only difference is the number of reduced datasets.
#
# Change between surveys per species and dataset: predicted occurrence
# probability in the present minus the past at each plot's elevation, averaged
# over all plots (as in 06). Per set: mean change across datasets and its
# 2.5-97.5% range. Compared: correlation of the mean changes (raw and relative
# to the mean across species), mean absolute difference, and agreement in the
# direction of change (range above zero, below zero, or overlapping zero).
#
# Inputs : pipeline/04_robustness/tmp_files/robustness_inputs.rds (00)
#          pipeline/04_robustness/outputs/rarefaction/rarefied_past_occurrence.rds (05)
#          pipeline/04_robustness/outputs/rarefaction/rarefaction_richness.csv (05)
#          pipeline/02_analysis/02_hmsc/data_clean/understory_reduced_datasets_oct09.csv
#            (the 10 datasets of the final HMSC models, run_oct09)
#          pipeline/01_data_prep/02_species_pool/outputs/understory_species_pool_strict_full_func_clim_cov.csv
# Outputs: pipeline/04_robustness/outputs/rarefaction/glm_reduced_10_vs_100.csv
#          pipeline/04_robustness/outputs/rarefaction/glm_reduced_10_vs_100_summary.csv
#          pipeline/04_robustness/figures/rarefaction/glm_reduced_10_vs_100.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)

tmp_dir <- here("pipeline", "04_robustness", "tmp_files")
out_dir <- here("pipeline", "04_robustness", "outputs", "rarefaction")
fig_dir <- here("pipeline", "04_robustness", "figures", "rarefaction")
walk(c(out_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line.x = element_line(colour = "black", linewidth = 0.6),
          axis.line.y = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black", size = 13),
          axis.title.x = element_text(colour = "black", size = 15, hjust = 1),
          axis.title.y = element_text(colour = "black", size = 15, hjust = 1),
          strip.text = element_text(colour = "black", face = "bold", size = 12),
          plot.margin  = margin(12, 10, 10, 10))
}

# -----------------------------------------------------------------------------
# 1. Inputs
# -----------------------------------------------------------------------------
inp      <- readRDS(file.path(tmp_dir, "robustness_inputs.rds"))
rarefied <- readRDS(file.path(out_dir, "rarefied_past_occurrence.rds"))
if (!"weighting" %in% names(rarefied)) rarefied$weighting <- "cover"
rarefied <- filter(rarefied, weighting == "cover")       # weighting used in HMSC
rarefied_survey <- read_csv(file.path(out_dir, "rarefaction_richness.csv"), show_col_types = FALSE)$rarefied_survey[1]
stopifnot("The HMSC reduction replaced the past survey" = rarefied_survey == "past")

reps_all <- sort(unique(rarefied$rep))
# The 10 datasets of the HMSC models: read from the file written by
# the HMSC launch script; if it is not available,
# reproduced with the same seed and draw (set.seed(20261009); sample of 10
# cover-weighted datasets), to be checked against the file on the machine
# where the models were built
used_file <- here("pipeline", "02_analysis", "02_hmsc", "data_clean",
                  c("understory_reduced_datasets_oct09.csv", "understory_reduced_datasets_used.csv"))
used_file <- used_file[file.exists(used_file)][1]
if (!is.na(used_file)) {
  used <- read_csv(used_file, show_col_types = FALSE)
  rep_col <- intersect(c("reduced_rep", "rep"), names(used))[1]
  stopifnot("No reduced_rep column in the datasets file" = !is.na(rep_col))
  reps_hmsc <- sort(used[[rep_col]])
  cat("HMSC datasets read from", basename(used_file), "\n")
} else {
  set.seed(20261009)
  reps_hmsc <- sort(sample(reps_all, 10))
  warning("understory_reduced_datasets_used.csv not found: HMSC datasets reproduced from the seed (",
          paste(reps_hmsc, collapse = ", "), "); check against the file on the modelling machine")
}
stopifnot("HMSC datasets not found among the reduced datasets" = all(reps_hmsc %in% reps_all))
cat("Reduced datasets:", length(reps_all), "| used in HMSC:", paste(reps_hmsc, collapse = ", "), "\n")

species <- read_csv(here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                         "understory_species_pool_strict_full_func_clim_cov.csv"),
                    show_col_types = FALSE) %>% pull(species) %>% unique()

plots <- inp$plot_info %>% filter(!is.na(elevation)) %>%
  mutate(elev_sc = as.numeric(scale(elevation))) %>% select(plot_id, elev_sc)
present <- inp$understory_long %>%
  filter(occurrence == 1, species %in% !!species, period == "present") %>% distinct(plot_id, species)

design <- bind_rows(plots %>% mutate(time = "past"), plots %>% mutate(time = "present")) %>%
  mutate(time = factor(time, levels = c("past", "present")))
X_past    <- model.matrix(~ time * elev_sc, mutate(plots, time = factor("past", levels = c("past", "present"))))
X_present <- model.matrix(~ time * elev_sc, mutate(plots, time = factor("present", levels = c("past", "present"))))

# Change in occurrence probability for one species and one reduced past dataset (as in 06)
fit_change <- function(past_records, sp) {
  y_past <- plots$plot_id %in% past_records$plot_id[past_records$species == sp]
  y_pres <- plots$plot_id %in% present$plot_id[present$species == sp]
  d <- mutate(design, y = as.integer(c(y_past, y_pres)))
  fit <- suppressWarnings(glm(y ~ time * elev_sc, family = binomial(link = "probit"), data = d))
  b <- coef(fit)
  if (anyNA(b)) return(NA_real_)
  mean(pnorm(X_present %*% b) - pnorm(X_past %*% b))
}

# -----------------------------------------------------------------------------
# 2. GLMs on each of the 100 reduced datasets
# -----------------------------------------------------------------------------
cat("Fitting", length(species), "species x", length(reps_all), "datasets...\n")
fits <- map_dfr(reps_all, function(r) {
  pr <- filter(rarefied, rep == r)
  tibble(rep = r, species = species, change = map_dbl(species, ~ fit_change(pr, .x)))
})

summarise_set <- function(d, set) d %>%
  group_by(species) %>%
  summarise(mean = mean(change, na.rm = TRUE),
            q025 = quantile(change, 0.025, na.rm = TRUE),
            q975 = quantile(change, 0.975, na.rm = TRUE), .groups = "drop") %>%
  mutate(rel = mean - mean(mean, na.rm = TRUE),
         direction = case_when(q025 > 0 ~ "increase", q975 < 0 ~ "decrease", TRUE ~ "overlaps zero"),
         set = set)
res <- bind_rows(summarise_set(fits, "all_100"),
                 summarise_set(filter(fits, rep %in% reps_hmsc), "hmsc_10"))
write_csv(res, file.path(out_dir, "glm_reduced_10_vs_100.csv"))

# -----------------------------------------------------------------------------
# 3. Comparison
# -----------------------------------------------------------------------------
w <- res %>% select(species, set, mean, rel, direction, q025, q975) %>%
  pivot_wider(names_from = set, values_from = c(mean, rel, direction, q025, q975))
summ <- tibble(
  n_species        = nrow(w),
  r_mean           = cor(w$mean_hmsc_10, w$mean_all_100, use = "complete.obs"),
  r_relative       = cor(w$rel_hmsc_10, w$rel_all_100, use = "complete.obs"),
  mean_abs_diff    = mean(abs(w$mean_hmsc_10 - w$mean_all_100), na.rm = TRUE),
  max_abs_diff     = max(abs(w$mean_hmsc_10 - w$mean_all_100), na.rm = TRUE),
  mean_change_100  = mean(w$mean_all_100, na.rm = TRUE),
  mean_change_10   = mean(w$mean_hmsc_10, na.rm = TRUE),
  same_direction   = sum(w$direction_hmsc_10 == w$direction_all_100, na.rm = TRUE),
  ratio_range_width = median((w$q975_hmsc_10 - w$q025_hmsc_10) / (w$q975_all_100 - w$q025_all_100), na.rm = TRUE)
) %>% mutate(across(where(is.double), ~ signif(.x, 3)))
write_csv(summ, file.path(out_dir, "glm_reduced_10_vs_100_summary.csv"))
cat("\nGLM species changes, 10 HMSC datasets vs all 100:\n"); print(as.data.frame(summ))
cat("\nDirection of change (rows: 100 datasets, columns: 10 datasets):\n")
print(table(w$direction_all_100, w$direction_hmsc_10))

# -----------------------------------------------------------------------------
# 4. Figure: mean change per species, 10 vs 100 datasets (2.5-97.5% ranges)
# -----------------------------------------------------------------------------
lim <- range(c(w$q025_all_100, w$q975_all_100, w$q025_hmsc_10, w$q975_hmsc_10), na.rm = TRUE)
fig <- ggplot(w, aes(mean_all_100, mean_hmsc_10)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  # 2.5-97.5% range across reduced datasets: round-ended lines, no end caps
  geom_segment(aes(x = q025_all_100, xend = q975_all_100, yend = mean_hmsc_10),
               colour = "#E07B39", alpha = 0.45, linewidth = 0.9, lineend = "round") +
  geom_segment(aes(xend = mean_all_100, y = q025_hmsc_10, yend = q975_hmsc_10),
               colour = "#E07B39", alpha = 0.45, linewidth = 0.9, lineend = "round") +
  geom_point(size = 2.6, colour = "#E07B39") +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 4.2,
           label = sprintf("r = %.2f\nmean |difference| = %.3f", summ$r_mean, summ$mean_abs_diff)) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Change in occurrence probability\n(GLM, mean over 100 reduced datasets)",
       y = "Change in occurrence probability\n(GLM, mean over 10 reduced datasets)") +
  theme_mora()
fig_file <- file.path(fig_dir, "glm_reduced_10_vs_100.png")
ggsave(fig_file, fig, width = 6.5, height = 6.5, dpi = 300, bg = "white")
cat("\nFigure written to", fig_dir, "\n")

# Open the figure (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)
