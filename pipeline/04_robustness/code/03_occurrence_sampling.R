# =============================================================================
# 03_occurrence_sampling.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Sensitivity of understory occurrence to the difference in sampled area
# between surveys: past occurrence was recorded over the whole plot (250, 500
# or 1000 m2), present occurrence over a 50-m2 census plus four 1-m2
# subplots. Two diagnostics:
#
#   1. Richness per plot in each survey (all species), and its relationship
#      with past plot area.
#
#   2. Area sensitivity per species, from the variation in plot size in the
#      past survey: probit GLM of past occurrence on log2(plot area) and
#      elevation (linear and quadratic). The area coefficient is the change in
#      probit occurrence per doubling of sampled area; species with a positive,
#      FDR-significant coefficient are flagged as sparse: they are the
#      species most likely to be missed when the sampled area decreases.
#
# Species-level results are reported for the evaluated set (00: more than
# five plots in each survey). The rarefaction of past occurrence is in 05-06.
#
# Inputs : pipeline/04_robustness/tmp_files/robustness_inputs.rds (00)
# Outputs: pipeline/04_robustness/outputs/occurrence/
#            richness_per_plot.csv, area_sensitivity_species.csv
#          pipeline/04_robustness/figures/occurrence/
#            richness.png            (richness per survey and by past plot area)
#            area_sensitivity.png    (area coefficient per species)
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

min_plots_area_model <- 10     # species occurring in at least this many past plots

tmp_dir <- here("pipeline", "04_robustness", "tmp_files")
out_dir <- here("pipeline", "04_robustness", "outputs", "occurrence")
fig_dir <- here("pipeline", "04_robustness", "figures", "occurrence")
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
col_past <- "#F4A259"; col_present <- "#B5501C"; col_flag <- c("Flagged" = "#9A3F12", "Not flagged" = "grey60")

inp <- readRDS(file.path(tmp_dir, "robustness_inputs.rds"))
ul  <- inp$understory_long %>% filter(occurrence == 1)
plots <- unique(inp$understory_long$plot_id)
n_plots <- length(plots)

# -----------------------------------------------------------------------------
# 1. Richness per plot
# -----------------------------------------------------------------------------
richness <- ul %>% count(period, plot_id, name = "richness") %>%
  complete(period, plot_id = plots, fill = list(richness = 0)) %>%
  left_join(inp$plot_info, by = "plot_id")
write_csv(richness, file.path(out_dir, "richness_per_plot.csv"))

mean_rich <- richness %>% group_by(period) %>% summarise(mean = mean(richness), sd = sd(richness), .groups = "drop")
cat("Mean species richness per plot:\n"); print(mean_rich)
cat("\nPast richness by plot area:\n")
print(richness %>% filter(period == "past") %>% group_by(plot_area) %>%
        summarise(n_plots = n(), mean_richness = mean(richness), .groups = "drop"))

# -----------------------------------------------------------------------------
# 2. Area sensitivity in the past survey
# -----------------------------------------------------------------------------
past_plots <- inp$plot_info %>% filter(!is.na(plot_area), !is.na(elevation)) %>%
  mutate(log2_area = log2(plot_area), elev_sc = as.numeric(scale(elevation)))
cat("\nPast plot areas:\n"); print(count(past_plots, plot_area))

past_occ <- ul %>% filter(period == "past") %>% distinct(plot_id, species) %>% mutate(y = 1L)
sp_area <- past_occ %>% count(species) %>% filter(n >= min_plots_area_model, species %in% inp$species_set) %>% pull(species)

area_fit <- map_dfr(sp_area, function(sp) {
  dd <- past_plots %>% left_join(past_occ %>% filter(species == sp), by = "plot_id") %>%
    mutate(y = coalesce(y, 0L))
  fit <- tryCatch(suppressWarnings(glm(y ~ log2_area + elev_sc + I(elev_sc^2),
                                       family = binomial(link = "probit"), data = dd)),
                  error = function(e) NULL)
  if (is.null(fit) || !"log2_area" %in% rownames(summary(fit)$coefficients)) return(NULL)
  cf <- summary(fit)$coefficients["log2_area", ]
  tibble(species = sp, n_past_plots = sum(dd$y), area_coef = cf[["Estimate"]],
         area_se = cf[["Std. Error"]], area_p = cf[["Pr(>|z|)"]])
}) %>%
  mutate(area_p_adj = p.adjust(area_p, method = "BH"),
         area_sensitive = area_coef > 0 & area_p_adj < 0.05) %>%
  left_join(inp$life_form, by = "species")
write_csv(area_fit, file.path(out_dir, "area_sensitivity_species.csv"))

cat("\nArea sensitivity (past survey), species with >=", min_plots_area_model, "past plots:", nrow(area_fit),
    "| sparse:", sum(area_fit$area_sensitive), "\n")
cat("Mean area coefficient (probit per doubling of area):", round(mean(area_fit$area_coef), 3), "\n")

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
p_a <- ggplot(richness, aes(period, richness, colour = period)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  scale_colour_manual(values = c(past = col_past, present = col_present), guide = "none") +
  scale_x_discrete(labels = c(past = "Past", present = "Present")) +
  labs(x = NULL, y = "Species richness per plot", title = "Richness per plot") +
  theme_mora()
p_b <- richness %>% filter(period == "past") %>%
  ggplot(aes(factor(plot_area), richness)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = col_past) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8, colour = col_past) +
  labs(x = "Past plot area (m²)", y = "Species richness per plot", title = "Past richness and plot area") +
  theme_mora()
p_d <- area_fit %>%
  mutate(flag = if_else(area_sensitive, "Flagged", "Not flagged")) %>%
  ggplot(aes(area_coef, reorder(species, area_coef), colour = flag)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = area_coef - 1.96 * area_se, xmax = area_coef + 1.96 * area_se), width = 0) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = col_flag, name = NULL) +
  labs(x = "Effect of doubling plot area on past occurrence (probit)", y = NULL,
       title = "Area sensitivity (past survey)") +
  theme_mora() + theme(axis.text.y = element_text(face = "italic", colour = "black", size = 8), legend.position = "bottom")

fig1 <- (p_a + p_b) + plot_annotation(tag_levels = "A")
ggsave(file.path(fig_dir, "richness.png"), fig1, width = 11, height = 5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "area_sensitivity.png"), p_d,
       width = 6.5, height = 2 + 0.16 * nrow(area_fit), dpi = 300, bg = "white", limitsize = FALSE)
cat("\nFigure written to", fig_dir, "\n")

# Panels for the combined diagnosis figures (08_diagnosis_figures.R)
panel_dir <- file.path(here("pipeline", "04_robustness", "tmp_files"), "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(richness = p_a, richness_area = p_b), file.path(panel_dir, "panels_03.rds"))
