# =============================================================================
# 05_rarefaction_occurrence.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Rarefaction of understory occurrence to equalise sampling effort between
# surveys, following Wilson et al. (2024, Ecology 105: e4242; Appendix S1,
# Section S1.3), who resurveyed historical plots with the same contemporary
# design (50-m2 census plus four 1-m2 cover subplots).
#
# Procedure (as in Wilson et al. 2024, applied to the survey with higher
# richness):
#   1. Species: those of the understory HMSC models (frequent pool with
#      complete traits and climatic affinity). Richness per plot is computed
#      over these species in each survey.
#   2. The difference in mean richness per plot between surveys (d, rounded)
#      is taken as a proxy for the difference in sampling effort and removed
#      from every plot of the richer survey. Plots with d species or fewer are
#      left unchanged.
#   3. Species retained in a plot are drawn without replacement with
#      probability proportional to their cover in that plot, so that species
#      of low cover are the most likely to be removed (weighting "cover", as
#      in Wilson et al. 2024). Missing cover is set to 0.1%. As a comparison,
#      species are also removed with equal probability (weighting "uniform"),
#      which tests how much the species-level results depend on the
#      assumption that low-cover species are the ones missed.
#   4. Repeated n_rare times for each weighting.
# In this study the past survey is the richer one: past occurrence was
# recorded over the whole plot (250-1000 m2), present occurrence over the
# 50-m2 census plus four 1-m2 subplots.
#
# The rarefied data sets are used for occurrence only: removing low-cover
# records would bias abundance conditional on presence.
#
# Per species: prevalence (proportion of plots) in the past, observed and
# rarefied (mean and 2.5-97.5% range over rarefactions), in the present, and
# the change between surveys before and after rarefaction.
#
# Inputs : pipeline/04_robustness/tmp_files/robustness_inputs.rds (00)
#          pipeline/01_data_prep/02_species_pool/outputs/understory_species_pool_strict_full_func_clim_cov.csv
# Outputs: pipeline/04_robustness/outputs/rarefaction/
#            rarefaction_richness.csv           (mean richness per survey, d)
#            rarefaction_species_change.csv     (one row per species x weighting)
#            rarefied_past_occurrence.rds       (rarefied records: weighting, rep, plot_id, species)
#          pipeline/04_robustness/figures/figS3_2_rarefaction_prevalence.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)

n_rare      <- 100
proxy_cover <- 0.1
weightings  <- c("cover", "uniform")
set.seed(20261008)

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
col_obs <- "grey60"; col_rare <- "#E07B39"

inp <- readRDS(file.path(tmp_dir, "robustness_inputs.rds"))
hmsc_species <- read_csv(here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                              "understory_species_pool_strict_full_func_clim_cov.csv"),
                         show_col_types = FALSE) %>%
  pull(species) %>% unique()
cat("HMSC understory species:", length(hmsc_species), "\n")

occ <- inp$understory_long %>%
  filter(occurrence == 1, species %in% hmsc_species) %>%
  distinct(period, plot_id, species, .keep_all = TRUE) %>%
  select(period, plot_id, species, cover)
plots   <- unique(inp$understory_long$plot_id)
n_plots <- length(plots)

# -----------------------------------------------------------------------------
# 1. Richness per plot and number of species to remove
# -----------------------------------------------------------------------------
richness <- occ %>% count(period, plot_id, name = "richness") %>%
  complete(period, plot_id = plots, fill = list(richness = 0))
rich_summary <- richness %>% group_by(period) %>%
  summarise(mean_richness = mean(richness), sd_richness = sd(richness), .groups = "drop")
richer <- rich_summary$period[which.max(rich_summary$mean_richness)]
d <- round(abs(diff(rich_summary$mean_richness)))
rich_summary <- rich_summary %>% mutate(rarefied_survey = richer, n_removed_per_plot = d)
write_csv(rich_summary, file.path(out_dir, "rarefaction_richness.csv"))
cat("Mean richness per plot (HMSC species):\n"); print(rich_summary)
cat("Rarefying the", richer, "survey:", d, "species removed per plot\n")

# -----------------------------------------------------------------------------
# 2. Rarefaction
# -----------------------------------------------------------------------------
rich_split <- occ %>% filter(period == richer) %>%
  mutate(weight = if_else(is.na(cover) | cover <= 0, proxy_cover, cover)) %>%
  split(.$plot_id)

rarefy_once <- function(weighting) {
  map_dfr(rich_split, function(x) {
    n <- nrow(x)
    if (n <= d) return(x[, c("plot_id", "species")])
    w <- if (weighting == "cover") x$weight else rep(1, n)
    keep <- sample(seq_len(n), n - d, prob = w)
    x[keep, c("plot_id", "species")]
  })
}
rarefied <- map_dfr(weightings, function(wt)
  map_dfr(seq_len(n_rare), ~ rarefy_once(wt) %>% mutate(weighting = wt, rep = .x, .before = 1)))
saveRDS(rarefied, file.path(out_dir, "rarefied_past_occurrence.rds"))

# -----------------------------------------------------------------------------
# 3. Change in prevalence per species, observed and rarefied
# -----------------------------------------------------------------------------
prev_obs <- occ %>% count(period, species) %>% mutate(prev = n / n_plots) %>%
  select(-n) %>% pivot_wider(names_from = period, values_from = prev, values_fill = 0) %>%
  complete(species = hmsc_species, fill = list(past = 0, present = 0))

prev_rare <- rarefied %>% count(weighting, rep, species) %>%
  complete(nesting(weighting, rep), species = hmsc_species, fill = list(n = 0)) %>%
  mutate(prev = n / n_plots) %>%
  group_by(weighting, species) %>%
  summarise(rarefied_mean = mean(prev), rarefied_q025 = quantile(prev, 0.025),
            rarefied_q975 = quantile(prev, 0.975), .groups = "drop")

species_change <- prev_rare %>%
  left_join(prev_obs, by = "species") %>%
  mutate(change_observed = present - past,
         change_rarefied = if (richer == "past") present - rarefied_mean else rarefied_mean - past,
         change_rarefied_q025 = if (richer == "past") present - rarefied_q975 else rarefied_q025 - past,
         change_rarefied_q975 = if (richer == "past") present - rarefied_q025 else rarefied_q975 - past) %>%
  left_join(inp$life_form, by = "species") %>%
  arrange(weighting, change_rarefied)
write_csv(species_change, file.path(out_dir, "rarefaction_species_change.csv"))

cat("\nChange in prevalence across species, by weighting:\n")
print(species_change %>% group_by(weighting) %>%
        summarise(mean_observed = mean(change_observed), mean_rarefied = mean(change_rarefied),
                  n_increase = sum(change_rarefied_q025 > 0), n_decrease = sum(change_rarefied_q975 < 0),
                  r_observed_rarefied = cor(change_observed, change_rarefied), .groups = "drop"))
cat("(2.5-97.5% range over rarefactions; reflects rarefaction variability only)\n")

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
ord <- species_change %>% filter(weighting == "cover") %>% arrange(change_rarefied) %>% pull(species)
fig <- species_change %>%
  mutate(species = factor(species, levels = ord)) %>%
  ggplot(aes(y = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(data = ~ distinct(.x, species, change_observed), aes(x = change_observed, colour = "Observed"), size = 1.8) +
  geom_linerange(aes(xmin = change_rarefied_q025, xmax = change_rarefied_q975, colour = weighting),
                 position = position_dodge(width = 0.6)) +
  geom_point(aes(x = change_rarefied, colour = weighting), size = 1.8, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(Observed = col_obs, cover = col_rare, uniform = "#1A759F"),
                      labels = c(Observed = "Observed", cover = "Reduced, cover-weighted", uniform = "Reduced, uniform"),
                      breaks = c("Observed", "cover", "uniform"), name = NULL) +
  labs(x = "Change in prevalence (present - past)", y = NULL) +
  theme_mora() +
  theme(axis.text.y = element_text(face = "italic", colour = "black", size = 8), legend.position = "top")
ggsave(file.path(sup_dir, "figS3_2_rarefaction_prevalence.png"), fig,
       width = 7.5, height = 2 + 0.18 * length(ord), dpi = 300, bg = "white", limitsize = FALSE)
cat("\nFigure written to", fig_dir, "\n")
