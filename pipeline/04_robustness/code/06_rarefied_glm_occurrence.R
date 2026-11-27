# =============================================================================
# 06_rarefied_glm_occurrence.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Species-level change in understory occurrence from single-species GLMs,
# before and after rarefaction of the richer survey (05; the past survey), for the species of the
# understory HMSC models. Following Wilson et al. (2024, Ecology 105: e4242),
# models are fitted to the observed data and to each rarefied data set, and
# results are averaged across rarefactions with their 2.5-97.5% range, for
# each weighting of 05 (cover-weighted and uniform removal).
#
# Model per species: probit GLM, presence ~ time * elevation (elevation centred
# and scaled across plots), the fixed-effect structure of the main HMSC model.
# Change between surveys: predicted occurrence probability in the present
# minus the past at each plot's elevation, averaged over all plots (as in
# 04_hmsc_results.R). Species with quasi-complete separation (very large
# coefficients or fitting warnings) are flagged.
#
# Also reported: change relative to the mean change across species; the
# correlation between observed and rarefied changes and with the HMSC
# full-model change (04); and the relationship of the adjustment introduced by
# rarefaction (rarefied - observed) with mean past cover and with the area
# coefficient of 03, which indicates whether the adjustment reflects the
# sampled-area effect or the weighting assumption.
#
# Inputs : pipeline/04_robustness/tmp_files/robustness_inputs.rds (00)
#          pipeline/04_robustness/outputs/rarefaction/rarefied_past_occurrence.rds (05)
#          pipeline/01_data_prep/02_species_pool/outputs/understory_species_pool_strict_full_func_clim_cov.csv
#          optional: pipeline/02_analysis/02_hmsc/outputs/<hmsc_file>
# Outputs: pipeline/04_robustness/outputs/rarefaction/
#            glm_occurrence_change_rarefied.csv (species x weighting)
#            glm_occurrence_change_rarefied_summary.csv
#            glm_occurrence_change_rarefied_robust_sign.csv (same sign under both weightings)
#          pipeline/04_robustness/figures/rarefaction/glm_vs_hmsc_observed.png (replaced by the pooled HMSC comparison, Figure S5.3)
#            (GLM vs HMSC, relative change)
#          pipeline/04_robustness/figures/rarefaction/
#            glm_occurrence_change_rarefied.png      (observed and both rarefactions)
#            rarefaction_adjustment_covariates.png   (adjustment vs past cover and area coefficient)
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)

# HMSC full-model species changes for comparison (04); set to NULL to skip
hmsc_file <- "understory_strict_medium_m0_time_nsp_rel_res_delta_species_full.csv"

tmp_dir <- here("pipeline", "04_robustness", "tmp_files")
out_dir <- here("pipeline", "04_robustness", "outputs", "rarefaction")
fig_dir <- here("pipeline", "04_robustness", "figures", "rarefaction")
sup_dir <- here("pipeline", "04_robustness", "figures")   # supplementary and diagnostic figures
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)
walk(c(out_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(),
          axis.line  = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 11),
          axis.title.y = element_text(colour = "black", size = 11),
          legend.text  = element_text(colour = "black"),
          plot.title   = element_text(colour = "black", size = 11, face = "bold"),
          plot.margin  = margin(12, 10, 10, 10))
}
col_obs <- "grey55"; col_rare <- "#E07B39"

inp      <- readRDS(file.path(tmp_dir, "robustness_inputs.rds"))
rarefied <- readRDS(file.path(out_dir, "rarefied_past_occurrence.rds"))
hmsc_species <- read_csv(here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                              "understory_species_pool_strict_full_func_clim_cov.csv"),
                         show_col_types = FALSE) %>% pull(species) %>% unique()

plots <- inp$plot_info %>% filter(!is.na(elevation)) %>%
  mutate(elev_sc = as.numeric(scale(elevation))) %>% select(plot_id, elev_sc)
occ <- inp$understory_long %>% filter(occurrence == 1, species %in% hmsc_species) %>%
  distinct(period, plot_id, species)
present <- occ %>% filter(period == "present") %>% select(plot_id, species)
past    <- occ %>% filter(period == "past") %>% select(plot_id, species)

# Design: every plot in both surveys
design <- bind_rows(plots %>% mutate(time = "past"), plots %>% mutate(time = "present")) %>%
  mutate(time = factor(time, levels = c("past", "present")))
X_past    <- model.matrix(~ time * elev_sc, plots %>% mutate(time = factor("past", levels = c("past", "present"))))
X_present <- model.matrix(~ time * elev_sc, plots %>% mutate(time = factor("present", levels = c("past", "present"))))

# Survey rarefied in 05 (the richer one); its records are replaced by the
# rarefied records, the other survey is used as observed
rarefied_survey <- read_csv(file.path(out_dir, "rarefaction_richness.csv"), show_col_types = FALSE)$rarefied_survey[1]
cat("Rarefied survey:", rarefied_survey, "\n")

fit_change <- function(records, sp) {
  past_records    <- if (rarefied_survey == "past") records else past
  present_records <- if (rarefied_survey == "present") records else present
  y_past <- plots$plot_id %in% past_records$plot_id[past_records$species == sp]
  y_pres <- plots$plot_id %in% present_records$plot_id[present_records$species == sp]
  d <- design %>% mutate(y = as.integer(c(y_past, y_pres)))
  warn <- FALSE
  fit <- withCallingHandlers(glm(y ~ time * elev_sc, family = binomial(link = "probit"), data = d),
                             warning = function(w) { warn <<- TRUE; invokeRestart("muffleWarning") })
  b <- coef(fit)
  if (anyNA(b)) return(tibble(change = NA_real_, separation = TRUE))
  tibble(change = mean(pnorm(X_present %*% b) - pnorm(X_past %*% b)),
         separation = warn || any(abs(b) > 8))
}

# -----------------------------------------------------------------------------
# 1. Observed data
# -----------------------------------------------------------------------------
observed_records <- if (rarefied_survey == "past") past else present
observed <- map_dfr(hmsc_species, ~ fit_change(observed_records, .x) %>% mutate(species = .x, .before = 1))

# -----------------------------------------------------------------------------
# 2. Rarefied data sets (each weighting of 05)
# -----------------------------------------------------------------------------
if (!"weighting" %in% names(rarefied)) rarefied$weighting <- "cover"
weightings <- unique(rarefied$weighting)
reps <- sort(unique(rarefied$rep))
cat("Fitting", length(hmsc_species), "species x", length(reps), "rarefied data sets x",
    length(weightings), "weightings...\n")
rare_fits <- map_dfr(weightings, function(wt) map_dfr(reps, function(r) {
  pr <- rarefied %>% filter(weighting == wt, rep == r)
  map_dfr(hmsc_species, ~ fit_change(pr, .x) %>% mutate(weighting = wt, rep = r, species = .x, .before = 1))
}))

rare_summary <- rare_fits %>%
  group_by(weighting, species) %>%
  summarise(q025 = quantile(change, 0.025, na.rm = TRUE),
            q975 = quantile(change, 0.975, na.rm = TRUE),
            change = mean(change, na.rm = TRUE),
            prop_separation = mean(separation),
            .groups = "drop")

# -----------------------------------------------------------------------------
# 3. Species covariates used to interpret the adjustment (rarefied - observed)
#    - mean cover of the species where present in the past survey (the
#      quantity that drives removal under cover weighting)
#    - area coefficient from 03 (sensitivity of past occurrence to plot area)
# -----------------------------------------------------------------------------
past_cover <- inp$understory_long %>%
  filter(period == "past", occurrence == 1, species %in% hmsc_species, !is.na(cover), cover > 0) %>%
  group_by(species) %>% summarise(mean_past_cover = mean(cover), .groups = "drop")
area_path <- here("pipeline", "04_robustness", "outputs", "occurrence", "area_sensitivity_species.csv")
area <- if (file.exists(area_path)) read_csv(area_path, show_col_types = FALSE) %>% select(species, area_coef) else
  tibble(species = character(), area_coef = numeric())

result <- rare_summary %>%
  left_join(observed %>% rename(change_observed = change, separation_observed = separation), by = "species") %>%
  rename(change_rarefied = change, change_rarefied_q025 = q025, change_rarefied_q975 = q975,
         prop_separation_rarefied = prop_separation) %>%
  group_by(weighting) %>%
  mutate(rel_observed = change_observed - mean(change_observed, na.rm = TRUE),
         rel_rarefied = change_rarefied - mean(change_rarefied, na.rm = TRUE),
         adjustment   = change_rarefied - change_observed) %>%
  ungroup() %>%
  left_join(past_cover, by = "species") %>%
  left_join(area, by = "species")

# HMSC full-model change (04), occurrence response
hmsc_path <- if (!is.null(hmsc_file)) here("pipeline", "02_analysis", "02_hmsc", "outputs", hmsc_file) else ""
has_hmsc <- file.exists(hmsc_path)
if (has_hmsc) {
  hmsc <- read_csv(hmsc_path, show_col_types = FALSE) %>%
    filter(response == "occurrence") %>%
    transmute(species = str_replace_all(species, "_", " "), change_hmsc = mean) %>%
    filter(species %in% hmsc_species) %>%
    mutate(rel_hmsc = change_hmsc - mean(change_hmsc))
  result <- result %>% left_join(hmsc, by = "species")
}
result <- result %>% left_join(inp$life_form, by = "species") %>% arrange(weighting, change_rarefied)
write_csv(result, file.path(out_dir, "glm_occurrence_change_rarefied.csv"))

# -----------------------------------------------------------------------------
# 4. Summary
# -----------------------------------------------------------------------------
cor_c <- function(x, y) round(cor(x, y, use = "complete.obs"), 3)
summ <- result %>% group_by(weighting) %>%
  summarise(mean_observed = mean(change_observed, na.rm = TRUE),
            mean_rarefied = mean(change_rarefied, na.rm = TRUE),
            n_above0 = sum(change_rarefied_q025 > 0, na.rm = TRUE),
            n_below0 = sum(change_rarefied_q975 < 0, na.rm = TRUE),
            n_overlap0 = sum(change_rarefied_q025 <= 0 & change_rarefied_q975 >= 0, na.rm = TRUE),
            r_obs_rarefied = cor_c(change_observed, change_rarefied),
            r_adjust_logcover = cor_c(adjustment, log10(mean_past_cover)),
            r_adjust_area = cor_c(adjustment, area_coef),
            r_hmsc_observed = if (has_hmsc) cor_c(change_hmsc, change_observed) else NA_real_,
            r_hmsc_rarefied = if (has_hmsc) cor_c(change_hmsc, change_rarefied) else NA_real_,
            .groups = "drop")
cat("\nSummary by weighting:\n"); print(as.data.frame(summ))
r_cover_uniform <- result %>% select(weighting, species, change_rarefied) %>%
  pivot_wider(names_from = weighting, values_from = change_rarefied)
if (all(c("cover", "uniform") %in% names(r_cover_uniform)))
  cat("Correlation cover-weighted vs uniform rarefied changes:",
      cor_c(r_cover_uniform$cover, r_cover_uniform$uniform), "\n")
cat("Species with separation in the observed fit:",
    sum(distinct(result, species, separation_observed)$separation_observed, na.rm = TRUE), "\n")
write_csv(summ, file.path(out_dir, "glm_occurrence_change_rarefied_summary.csv"))

# Species whose rarefied change has the same sign under both weightings, with
# the 2.5-97.5% range across rarefactions excluding zero in both
robust_sign <- result %>%
  group_by(species) %>%
  summarise(robust_increase = all(change_rarefied_q025 > 0),
            robust_decrease = all(change_rarefied_q975 < 0), .groups = "drop") %>%
  mutate(robust_direction = case_when(robust_increase ~ "increase", robust_decrease ~ "decrease",
                                      TRUE ~ "not robust"))
write_csv(robust_sign, file.path(out_dir, "glm_occurrence_change_rarefied_robust_sign.csv"))
cat("\nRarefied change with consistent sign under both weightings:\n")
print(count(robust_sign, robust_direction))
cat("Increase:", paste(robust_sign$species[robust_sign$robust_direction == "increase"], collapse = ", "), "\n")
cat("Decrease:", paste(robust_sign$species[robust_sign$robust_direction == "decrease"], collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 5. Figures
# -----------------------------------------------------------------------------
col_w <- c(Observed = col_obs, cover = col_rare, uniform = "#1A759F")
lab_w <- c(Observed = "Observed", cover = "Reduced, cover-weighted", uniform = "Reduced, uniform")

# (a) Species change: observed and both rarefactions
ord <- result %>% filter(weighting == weightings[1]) %>% arrange(change_rarefied) %>% pull(species)
fig_a <- result %>%
  mutate(species = factor(species, levels = ord)) %>%
  ggplot(aes(y = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(data = ~ distinct(.x, species, change_observed), aes(x = change_observed, colour = "Observed"), size = 1.8) +
  geom_linerange(aes(xmin = change_rarefied_q025, xmax = change_rarefied_q975, colour = weighting),
                 position = position_dodge(width = 0.6)) +
  geom_point(aes(x = change_rarefied, colour = weighting), size = 1.8, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = col_w, labels = lab_w, breaks = c("Observed", weightings), name = NULL) +
  labs(x = "Change in occurrence probability (single-species GLM, present - past)", y = NULL) +
  theme_mora() +
  theme(axis.text.y = element_text(face = "italic", colour = "black", size = 8), legend.position = "top")
ggsave(file.path(fig_dir, "glm_occurrence_change_rarefied.png"), fig_a,
       width = 7.5, height = 2 + 0.18 * length(ord), dpi = 300, bg = "white", limitsize = FALSE)

# (b) GLM vs HMSC: change relative to the mean across species
if (has_hmsc) {
  cmp <- bind_rows(
    result %>% distinct(species, rel_hmsc, rel_observed) %>% transmute(species, rel_hmsc, rel_glm = rel_observed, set = "Observed"),
    result %>% transmute(species, rel_hmsc, rel_glm = rel_rarefied, set = weighting)) %>%
    mutate(set = factor(set, levels = c("Observed", weightings), labels = lab_w[c("Observed", weightings)]))
  r_lab <- cmp %>% group_by(set) %>% summarise(r = cor(rel_hmsc, rel_glm, use = "complete.obs"), .groups = "drop") %>%
    mutate(label = paste0("r = ", sprintf("%.2f", r)))
  lim <- range(c(cmp$rel_hmsc, cmp$rel_glm), na.rm = TRUE)
  fig_b <- ggplot(cmp, aes(rel_hmsc, rel_glm, colour = set)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(size = 1.8, alpha = 0.85) +
    geom_text(data = r_lab, aes(x = lim[1], y = lim[2], label = label), hjust = 0, vjust = 1,
              colour = "black", size = 3.5, inherit.aes = FALSE) +
    facet_wrap(~ set, nrow = 1) +
    scale_colour_manual(values = setNames(unname(col_w[c("Observed", weightings)]), levels(cmp$set)), guide = "none") +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = "HMSC full-model change (relative to mean across species)",
         y = "GLM change (relative to mean across species)") +
    theme_mora() + theme(strip.text = element_text(colour = "black", face = "bold", size = 10))
  ggsave(file.path(fig_dir, "glm_vs_hmsc_observed.png"), fig_b,
         width = 3.3 * (1 + length(weightings)), height = 4, dpi = 300, bg = "white")
}

# (c) Adjustment introduced by rarefaction vs mean past cover and area coefficient
adj <- result %>%
  transmute(weighting, species, adjustment, mean_past_cover = log10(mean_past_cover), area_coef) %>%
  pivot_longer(c(mean_past_cover, area_coef), names_to = "covariate", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(covariate = recode(covariate, mean_past_cover = "Mean past cover (log10 %)",
                            area_coef = "Area coefficient (probit per doubling)"))
fig_c <- ggplot(adj, aes(value, adjustment, colour = weighting)) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
  geom_point(size = 1.6, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, formula = y ~ x) +
  facet_wrap(~ covariate, scales = "free_x") +
  scale_colour_manual(values = col_w[weightings], labels = lab_w[weightings], name = NULL) +
  labs(x = NULL, y = "Reduced - observed change") +
  theme_mora() + theme(legend.position = "top", strip.text = element_text(colour = "black", face = "bold", size = 10))
ggsave(file.path(fig_dir, "rarefaction_adjustment_covariates.png"), fig_c,
       width = 8, height = 4.2, dpi = 300, bg = "white")

cat("\nFigures written to", fig_dir, "\n")
