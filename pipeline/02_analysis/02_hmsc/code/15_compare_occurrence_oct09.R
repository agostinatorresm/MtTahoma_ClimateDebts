# =============================================================================
# 15_compare_occurrence_oct09.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (after 04_hmsc_results.R for the
#           observed-data and pooled reduced-data occurrence models)
# Author: Agostina Torres
#
# Understory occurrence: comparison of the HMSC model fitted to the observed
# data with the pooled posterior of the models fitted to 10 past datasets
# reduced to the present sampling effort (main estimate), and with single-
# species GLMs over the 100 reductions (04_robustness/06).
#
#   1. Species change in occurrence (full model, all random effects; change in
#      occurrence probability averaged over plots), relative to the mean
#      across species: observed vs pooled HMSC, and pooled HMSC vs GLM on the
#      cover-weighted reductions. Relative values are compared because the
#      change common to all species is removed by the reduction by
#      construction.
#   2. Agreement in direction: species with supported change (posterior
#      support >= 0.90) in the pooled model, and their agreement with the
#      GLM species of consistent sign (Table SR1).
#   3. Trait effects on the time response (Gamma, occurrence): observed vs
#      pooled, posterior mean and 90% interval.
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/understory_strict_long_<model>_rel_res_delta_species_full.csv
#          pipeline/02_analysis/02_hmsc/outputs/understory_strict_long_<model>_rel_res_gamma_time.csv
#            (<model> = m1_time_elev_nsp_nos, m1_time_elev_nsp_nos_redpool)
#          pipeline/04_robustness/outputs/rarefaction/glm_occurrence_change_rarefied.csv
#          pipeline/04_robustness/outputs/rarefaction/glm_occurrence_change_rarefied_robust_sign.csv
# Outputs: pipeline/02_analysis/02_hmsc/outputs/run_oct09_compare/*.csv
#          pipeline/02_analysis/02_hmsc/figures/run_oct09_compare/*.png
#          pipeline/04_robustness/figures/figSR4_hmsc_occurrence_observed_vs_reduced.png
#          pipeline/04_robustness/figures/figSR6_gamma_occurrence_observed_vs_reduced.png
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)
library(ggpubr)

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
res_dir    <- file.path(module_dir, "outputs")
out_dir    <- file.path(res_dir, "run_oct09_compare")
fig_dir    <- file.path(module_dir, "figures", "run_oct09_compare")
walk(c(out_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)
glm_dir    <- here("pipeline", "04_robustness", "outputs", "rarefaction")
sup_dir <- here("pipeline", "04_robustness", "figures")   # supplementary and diagnostic figures
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)

models <- c(observed = "m1_time_elev_nsp_nos", reduced_pooled = "m1_time_elev_nsp_nos_redpool")
res_path <- function(model, name) {
  file.path(res_dir, paste0("understory_strict_long_", model, "_rel_res_", name))
}

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
          strip.text   = element_text(colour = "black", face = "bold", size = 10),
          plot.title   = element_text(colour = "black", size = 11, face = "bold"),
          plot.margin  = margin(12, 10, 10, 10))
}
col_model <- c(observed = "grey55", reduced_pooled = "#E07B39")
lab_model <- c(observed = "Observed data", reduced_pooled = "Reduced past data (10 fits pooled)")
trait_labels <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity",
                  SLA = "SLA", seed_mass = "Seed mass", woodinesswoody = "Woody (vs non-woody)")
sp_key <- function(x) gsub("[ _]+", "_", trimws(x))

# -----------------------------------------------------------------------------
# 1. Species change in occurrence
# -----------------------------------------------------------------------------
delta <- imap_dfr(models, function(m, nm) {
  read_csv(res_path(m, "delta_species_full.csv"), show_col_types = FALSE) %>%
    filter(response == "occurrence") %>%
    mutate(model = nm)
}) %>%
  mutate(species = sp_key(species)) %>%
  group_by(model) %>%
  mutate(rel = mean - mean(mean)) %>%
  ungroup()

cat("Mean change in occurrence across species:\n")
print(delta %>% group_by(model) %>% summarise(mean_change = round(mean(mean), 3), .groups = "drop"))

glm <- read_csv(file.path(glm_dir, "glm_occurrence_change_rarefied.csv"), show_col_types = FALSE) %>%
  mutate(species = sp_key(species))
glm_cover <- glm %>% filter(weighting == "cover") %>%
  select(species, glm_change = change_rarefied, glm_q025 = change_rarefied_q025,
         glm_q975 = change_rarefied_q975, glm_observed = change_observed) %>%
  mutate(glm_rel = glm_change - mean(glm_change, na.rm = TRUE),
         glm_obs_rel = glm_observed - mean(glm_observed, na.rm = TRUE))
robust <- read_csv(file.path(glm_dir, "glm_occurrence_change_rarefied_robust_sign.csv"), show_col_types = FALSE)
robust_sp <- sp_key(robust$species[robust$robust_direction != "not robust"])

wide <- delta %>%
  select(species, model, mean, rel, q05, q95, support, direction) %>%
  pivot_wider(names_from = model, values_from = c(mean, rel, q05, q95, support, direction)) %>%
  left_join(glm_cover, by = "species") %>%
  left_join(robust %>% transmute(species = sp_key(species), glm_robust_direction = robust_direction),
            by = "species") %>%
  mutate(glm_robust = species %in% robust_sp)

cors <- tibble(
  comparison = c("HMSC observed vs HMSC reduced (pooled)",
                 "HMSC reduced (pooled) vs GLM reduced (cover, 100)",
                 "HMSC observed vs GLM observed"),
  r = c(cor(wide$rel_observed, wide$rel_reduced_pooled, use = "complete.obs"),
        cor(wide$rel_reduced_pooled, wide$glm_rel, use = "complete.obs"),
        cor(wide$rel_observed, wide$glm_obs_rel, use = "complete.obs"))
) %>% mutate(r = round(r, 3))
print(cors)

# -----------------------------------------------------------------------------
# 2. Direction of change
# -----------------------------------------------------------------------------
direction_tbl <- wide %>%
  transmute(species,
            pooled_direction   = direction_reduced_pooled,
            observed_direction = direction_observed,
            pooled_support     = round(support_reduced_pooled, 3),
            glm_robust,
            glm_sign = glm_robust_direction) %>%
  arrange(pooled_direction, desc(pooled_support))
cat("\nSupported changes (posterior support >= 0.90):\n")
print(count(direction_tbl, pooled_direction, observed_direction))
agree <- direction_tbl %>%
  filter(glm_robust) %>%
  mutate(agrees = case_when(pooled_direction == "not supported" ~ "pooled: not supported",
                            as.character(pooled_direction) == glm_sign ~ "same direction",
                            TRUE ~ "opposite direction"))
cat("\nGLM species of consistent sign (Table SR1) in the pooled HMSC model:\n")
print(count(agree, agrees))

write_csv(wide, file.path(out_dir, "occurrence_change_observed_vs_reduced.csv"))
write_csv(direction_tbl, file.path(out_dir, "occurrence_direction_pooled.csv"))
write_csv(cors, file.path(out_dir, "occurrence_correlations.csv"))

# -----------------------------------------------------------------------------
# 3. Gamma (occurrence): observed vs pooled
# -----------------------------------------------------------------------------
gamma <- imap_dfr(models, function(m, nm) {
  read_csv(res_path(m, "gamma_time.csv"), show_col_types = FALSE) %>%
    filter(response == "occurrence") %>%
    mutate(model = nm)
}) %>%
  mutate(trait_lab = factor(coalesce(trait_labels[trait], trait), levels = rev(unique(coalesce(trait_labels[trait], trait)))))
print(gamma %>% select(model, trait, mean, q05, q95, support) %>% mutate(across(where(is.numeric), ~ round(.x, 3))))
write_csv(gamma, file.path(out_dir, "gamma_occurrence_observed_vs_reduced.csv"))

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
pa <- ggplot(wide, aes(rel_observed, rel_reduced_pooled)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(shape = glm_robust), size = 2, colour = "grey20") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "Other", `TRUE` = "Consistent sign (GLM)"), name = NULL) +
  labs(x = "HMSC, observed data", y = "HMSC, reduced past data (pooled)",
       title = paste0("Change in occurrence relative to mean (r = ", cors$r[1], ")")) +
  theme_mora() + theme(legend.position = "bottom")

pb <- ggplot(wide, aes(glm_rel, rel_reduced_pooled)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(shape = glm_robust), size = 2, colour = "grey20") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  labs(x = "GLM, 100 cover-weighted reductions (mean)", y = "HMSC, reduced past data (pooled)",
       title = paste0("HMSC vs GLM on reduced data (r = ", cors$r[2], ")")) +
  theme_mora()

pc <- ggplot(gamma, aes(mean, trait_lab, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = q05, xmax = q95), width = 0, linewidth = 0.8,
                position = position_dodge(width = 0.5), orientation = "y") +
  geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = col_model, labels = lab_model, name = NULL) +
  labs(x = "Gamma (trait effect on the occurrence response to time)", y = NULL,
       title = "Trait effects, occurrence") +
  theme_mora() + theme(legend.position = "bottom")

fig <- ggarrange(pa, pb, pc, ncol = 3, labels = c("A", "B", "C"), widths = c(1, 1, 1.1))
ggsave(file.path(fig_dir, "occurrence_observed_vs_reduced.png"), fig,
       width = 16, height = 5.5, dpi = 300, bg = "white")
# Diagnostic figures: species change (panels A-B) and Gamma (panel C); final version in Figure S5.3 (05f)
ggsave(file.path(sup_dir, "figSR4_hmsc_occurrence_observed_vs_reduced.png"),
       ggarrange(pa, pb, ncol = 2, labels = c("A", "B")),
       width = 11, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(sup_dir, "figSR6_gamma_occurrence_observed_vs_reduced.png"), pc,
       width = 6.5, height = 5, dpi = 300, bg = "white")

ord <- wide %>% arrange(mean_reduced_pooled) %>% pull(species)
pd <- delta %>%
  mutate(species = factor(str_replace_all(species, "_", " "), levels = str_replace_all(ord, "_", " "))) %>%
  ggplot(aes(mean, species, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = q05, xmax = q95), width = 0, alpha = 0.7,
                position = position_dodge(width = 0.6), orientation = "y") +
  geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = col_model, labels = lab_model, name = NULL) +
  labs(x = "Change in occurrence probability (present - past; mean and 90% interval)", y = NULL) +
  theme_mora() +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "italic", size = 8))
ggsave(file.path(fig_dir, "occurrence_species_observed_vs_reduced.png"), pd,
       width = 8, height = 2 + 0.17 * length(ord), dpi = 300, bg = "white", limitsize = FALSE)

cat("\nFigures written to", fig_dir, "\n")
