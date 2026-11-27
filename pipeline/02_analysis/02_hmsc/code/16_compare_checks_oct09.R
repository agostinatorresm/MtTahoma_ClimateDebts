# =============================================================================
# 16_compare_checks_oct09.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine
#           after 04_hmsc_results.R for "rel" and "relall")
# Author: Agostina Torres
#
# Robustness checks of the understory models (run of 2026-10-09):
#   (a) Relative cover paired by plot (main) vs all presences with cover
#       (census-only occurrences missing at both times): species change in
#       log10 relative cover (full model) and trait effects on the time
#       response (Gamma), from the outputs of 04_hmsc_results.R.
#   (b) Models without design-sensitive species vs the corresponding main
#       models fitted to the same data:
#         occurrence without the 10 sparse species (medium) vs the
#           occurrence model on the observed data (long)
#         paired relative cover without the 18 patchy species (medium) vs the
#           paired relative-cover model (long)
#       compared through Gamma for time and, for the species in both models,
#       Beta for time (change at mean elevation, link scale).
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/understory_strict_long_m1_time_elev_nsp_nos_redpool_<rel|relall>_res_*.csv (04)
#          pipeline/02_analysis/02_hmsc/outputs/run_oct09/understory_strict_<run_mode>_<model_id>.rds
# Outputs: pipeline/02_analysis/02_hmsc/outputs/run_oct09_compare/checks_*.csv
#          pipeline/02_analysis/02_hmsc/figures/run_oct09_compare/checks_*.png
#          pipeline/04_robustness/figures/figSR8_hmsc_checks.png
# =============================================================================

library(here)
library(tidyverse)
library(Hmsc)
library(hrbrthemes)
library(ggpubr)

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
res_dir    <- file.path(module_dir, "outputs")
fit_dir    <- file.path(res_dir, "run_oct09")
out_dir    <- file.path(res_dir, "run_oct09_compare")
fig_dir    <- file.path(module_dir, "figures", "run_oct09_compare")
walk(c(out_dir, fig_dir), dir.create, recursive = TRUE, showWarnings = FALSE)
sig_level  <- 0.90
sup_dir <- here("pipeline", "04_robustness", "figures")   # supplementary and diagnostic figures
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)

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
          legend.text  = element_text(colour = "black"),
          strip.text   = element_text(colour = "black", face = "bold", size = 10),
          plot.title   = element_text(colour = "black", size = 11, face = "bold"),
          plot.margin  = margin(12, 10, 10, 10))
}
trait_labels <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity",
                  SLA = "SLA", seed_mass = "Seed mass", woodinesswoody = "Woody (vs non-woody)")
trait_lab <- function(x) factor(coalesce(trait_labels[x], x), levels = rev(trait_labels))

summ <- function(s) {
  tibble(mean = colMeans(s), q05 = apply(s, 2, quantile, 0.05), q95 = apply(s, 2, quantile, 0.95),
         support = colMeans(s > 0)) %>%
    mutate(direction = case_when(support >= sig_level ~ "increase",
                                 support <= 1 - sig_level ~ "decrease", TRUE ~ "not supported"))
}

# -----------------------------------------------------------------------------
# (a) Paired vs all presences (relative cover)
# -----------------------------------------------------------------------------
res_path <- function(resp, name) {
  file.path(res_dir, paste0("understory_strict_long_m1_time_elev_nsp_nos_redpool_", resp, "_res_", name))
}
cover_variants <- c(paired = "rel", all_presences = "relall")

delta_cov <- imap_dfr(cover_variants, function(resp, nm) {
  read_csv(res_path(resp, "delta_species_full.csv"), show_col_types = FALSE) %>%
    filter(response == "log10_abundance") %>% mutate(variant = nm)
})
gamma_cov <- imap_dfr(cover_variants, function(resp, nm) {
  read_csv(res_path(resp, "gamma_time.csv"), show_col_types = FALSE) %>%
    filter(response == "abundance") %>% mutate(variant = nm)
})

delta_cov_w <- delta_cov %>%
  select(species, variant, mean, direction) %>%
  pivot_wider(names_from = variant, values_from = c(mean, direction)) %>%
  mutate(rel_paired = mean_paired - mean(mean_paired, na.rm = TRUE),
         rel_all    = mean_all_presences - mean(mean_all_presences, na.rm = TRUE))
r_cov <- round(cor(delta_cov_w$rel_paired, delta_cov_w$rel_all, use = "complete.obs"), 3)
cat("(a) Relative cover, paired vs all presences\n")
print(delta_cov %>% group_by(variant) %>%
        summarise(n_species = n(), mean_change = round(mean(mean), 3),
                  n_increase = sum(direction == "increase"), n_decrease = sum(direction == "decrease"),
                  .groups = "drop"))
cat("Correlation of species change (relative to the mean):", r_cov, "\n")
print(gamma_cov %>% select(variant, trait, mean, q05, q95, support) %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

# -----------------------------------------------------------------------------
# (b) Models without design-sensitive species
# -----------------------------------------------------------------------------
read_fit <- function(run_mode, id) readRDS(file.path(fit_dir, paste0("understory_strict_", run_mode, "_", id, ".rds")))
post_time <- function(m) {
  post <- poolMcmcChains(m$postList)
  i_t  <- which(m$covNames == "timepresent")
  B <- t(sapply(post, function(p) p$Beta[i_t, ]));  colnames(B) <- m$spNames
  G <- t(sapply(post, function(p) p$Gamma[i_t, ])); colnames(G) <- m$trNames
  list(beta = B, gamma = G[, colnames(G) != "(Intercept)", drop = FALSE])
}

checks <- list(
  occurrence_no_area_sensitive = list(main = c("long", "m1_time_elev_nsp_nos_pa"),
                                      check = c("medium", "m1_time_elev_nsp_nos_noarea_pa")),
  relcover_no_patchy           = list(main = c("long", "m1_time_elev_nsp_nos_rel"),
                                      check = c("medium", "m1_time_elev_nsp_nos_nonpatchy_rel"))
)

check_res <- imap(checks, function(ck, nm) {
  pm <- post_time(read_fit(ck$main[1],  ck$main[2]));  gc()
  pc <- post_time(read_fit(ck$check[1], ck$check[2])); gc()
  gamma <- bind_rows(
    summ(pm$gamma) %>% mutate(trait = colnames(pm$gamma), model = "main"),
    summ(pc$gamma) %>% mutate(trait = colnames(pc$gamma), model = "check")
  ) %>% mutate(check = nm)
  shared <- intersect(colnames(pm$beta), colnames(pc$beta))
  beta <- tibble(species = shared,
                 main  = colMeans(pm$beta[, shared, drop = FALSE]),
                 check = colMeans(pc$beta[, shared, drop = FALSE])) %>%
    mutate(check_name = nm)
  list(gamma = gamma, beta = beta,
       r_beta = round(cor(beta$main, beta$check), 3),
       n_removed = length(setdiff(colnames(pm$beta), shared)))
})

gamma_chk <- bind_rows(map(check_res, "gamma"))
beta_chk  <- bind_rows(map(check_res, "beta"))
cat("\n(b) Models without design-sensitive species\n")
iwalk(check_res, ~ cat(.y, ": species removed", .x$n_removed, "| r of Beta(time) for shared species", .x$r_beta, "\n"))
print(gamma_chk %>% select(check, model, trait, mean, q05, q95, support) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf)

write_csv(delta_cov_w, file.path(out_dir, "checks_relcover_paired_vs_all.csv"))
write_csv(gamma_cov,   file.path(out_dir, "checks_gamma_relcover_paired_vs_all.csv"))
write_csv(gamma_chk,   file.path(out_dir, "checks_gamma_design_sensitive.csv"))
write_csv(beta_chk,    file.path(out_dir, "checks_beta_design_sensitive.csv"))

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
col_var <- c(paired = "#9A3F12", all_presences = "grey55", main = "black", check = "#1A759F")
gamma_panel <- function(d, var, title, labels) {
  ggplot(d, aes(mean, trait_lab(trait), colour = .data[[var]])) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = q05, xmax = q95), width = 0, linewidth = 0.8,
                  position = position_dodge(width = 0.5), orientation = "y") +
    geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = col_var, labels = labels, name = NULL) +
    labs(x = "Gamma (trait effect on the response to time)", y = NULL, title = title) +
    theme_mora() + theme(legend.position = "bottom")
}
pa <- ggplot(delta_cov_w, aes(rel_all, rel_paired)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2, colour = "grey20") +
  labs(x = "All presences with cover", y = "Paired by plot (main)",
       title = paste0("Change in relative cover, relative to mean (r = ", r_cov, ")")) +
  theme_mora()
pb <- gamma_panel(gamma_cov, "variant", "Relative cover: paired vs all presences",
                  c(paired = "Paired by plot (main)", all_presences = "All presences with cover"))
pc <- gamma_panel(filter(gamma_chk, check == "occurrence_no_area_sensitive"), "model",
                  "Occurrence without sparse species",
                  c(main = "All species (observed data)", check = "Without 10 sparse species"))
pd <- gamma_panel(filter(gamma_chk, check == "relcover_no_patchy"), "model",
                  "Relative cover without patchy species",
                  c(main = "All species (paired)", check = "Without 18 patchy species"))
fig <- ggarrange(pa, pb, pc, pd, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(file.path(fig_dir, "checks_cover_and_design_sensitive.png"), fig,
       width = 13, height = 10, dpi = 300, bg = "white")
ggsave(file.path(sup_dir, "figSR8_hmsc_checks.png"), fig, width = 13, height = 10, dpi = 300, bg = "white")
cat("\nFigure written to", fig_dir, "\n")
