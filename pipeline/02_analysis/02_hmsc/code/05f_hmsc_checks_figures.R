# =============================================================================
# 05f_hmsc_checks_figures.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the local computer, after
#           04_hmsc_results.R, 16_compare_checks_oct09.R,
#           04_robustness/code/11_glm_reduced_10_vs_100.R and
#           28_pool_compare_sparse12.R)
# Author: Agostina Torres
#
# Figures of the checks of the understory HMSC models (Supporting Information
# S5), one per model, with species changes and trait effects:

#   Figure S5.3, occurrence model after reduction of past data
#     (A) change in occurrence probability per species (species on the y
#         axis, posterior mean and 90% interval): model fitted to the observed
#         data vs pooled model fitted to 10 reduced past datasets (16)
#     (B) trait effects on the occurrence response to time (Gamma), same models
#     (C) species change from single-species GLMs on reduced past data: mean
#         over the 10 reductions used in HMSC vs mean over all 100
#     (D) species change without the 12 sparse species (10 species with
#         sparse past records plus Veratrum viride and Epilobium
#         angustifolium): change in occurrence probability at mean elevation
#         for the species in both pooled reduced models, relative to the mean
#         across species, with 90% intervals (28_pool_compare_sparse12.R)
#   Figure S5.4, relative-cover model after pairing by plot
#     (A) change in relative cover per species (species on the y axis):
#         paired by plot (main) vs all presences with cover (16)
#     (B) Gamma, same models
#     (C) species change without the patchy species (Beta for time, as in S5.3D)
#
# In the scatter panels (S5.3D, S5.4C) species changes are shown relative to
# the mean across species, so that differences among species are compared.
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/
#            understory_strict_long_m1_time_elev_nsp_nos_rel_res_{delta_species_plotRE,gamma_time}.csv (04; observed data)
#            understory_strict_long_m1_time_elev_nsp_nos_redpool_rel_res_{delta_species_plotRE,gamma_time}.csv (04; pooled)
#            (searched in outputs/ and its subfolders; occurrence rows used)
#          pipeline/02_analysis/02_hmsc/outputs/run_oct09_compare/ (16)
#            checks_relcover_paired_vs_all.csv, checks_gamma_relcover_paired_vs_all.csv,
#            checks_beta_design_sensitive.csv
#          pipeline/04_robustness/outputs/rarefaction/glm_reduced_10_vs_100.csv (11)
#          pipeline/02_analysis/02_hmsc/outputs/run_sparse12/sparse12_species_change_comparison.csv (28)
# Outputs: pipeline/02_analysis/02_hmsc/figures/m1_final/ms1__figS5_3_occurrence_checks_m1_final.png
#          pipeline/02_analysis/02_hmsc/figures/m1_final/ms1__figS5_4_cover_checks_m1_final.png
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)
library(patchwork)

tag        <- "m1_final"
fig_prefix <- "ms1__"
sig_level  <- 0.90
col_occ    <- "#E07B39"        # occurrence (as in the species-change figures)
col_abu    <- "#9A3F12"        # abundance
col_cmp    <- "black"          # comparison models (observed data; all presences with cover)
col_abu_raw <- "grey70"        # relative cover, paired by plot, not centred
# the final estimates (reduced occurrence; centred relative cover) in orange
trait_labels <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity",
                  SLA = "SLA", seed_mass = "Seed mass", woodinesswoody = "Woody (vs non-woody)")

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  theme_ipsum(...) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line.x = element_line(colour = "black", linewidth = 0.6),
          axis.line.y = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black", size = 12),
          axis.title.x = element_text(colour = "black", size = 13, hjust = 1),
          axis.title.y = element_text(colour = "black", size = 13, hjust = 1),
          legend.text = element_text(colour = "black", size = 12),
          plot.title = element_text(colour = "black", size = 13, face = "bold"),
          plot.margin = margin(12, 10, 10, 10))
}

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
res_dir    <- file.path(module_dir, "outputs")
cmp_dir    <- file.path(res_dir, "run_oct09_compare")
fig_dir    <- file.path(module_dir, "figures", tag)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

need <- function(path) {
  if (!file.exists(path)) stop("File not found (run the script that writes it, or copy it from the computing machine): ", path)
  path
}
# 04 result tables of a model pair (occurrence _pa + relative cover _rel; one
# table with a response column, named <model>_rel_res_<name>), searched in
# outputs/ and its subfolders; NA if not found
res_file <- function(model, name, abu = "rel") {
  pattern <- paste0("^understory_strict_long_", model, "_", abu, "_res_", name, "$")
  hits <- list.files(res_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) { warning("No file matching ", pattern, " in ", res_dir); return(NA_character_) }
  if (length(hits) > 1) message("Several matches for ", pattern, "; using ", hits[1])
  hits[1]
}
read_res <- function(model, name, abu = "rel") {
  f <- res_file(model, name, abu)
  if (is.na(f)) return(NULL)
  read_csv(f, show_col_types = FALSE)
}
occ_only <- function(d) if ("response" %in% names(d)) filter(d, response == "occurrence") else d

# -----------------------------------------------------------------------------
# 1. Inputs
# -----------------------------------------------------------------------------
models_occ <- c(observed = "m1_time_elev_nsp_nos", reduced = "m1_time_elev_nsp_nos_redpool")

delta_occ <- imap(models_occ, ~ read_res(.x, "delta_species_plotRE.csv"))
gamma_occ <- imap(models_occ, ~ read_res(.x, "gamma_time.csv"))
# Panels A and B need the 04 tables of both models; without the observed-data
# model they are drawn as empty panels
have_obs <- !any(map_lgl(c(delta_occ, gamma_occ), is.null))
if (have_obs) {
  delta_occ <- imap_dfr(delta_occ, ~ occ_only(.x) %>% mutate(model = .y)) %>%
    group_by(model) %>% mutate(rel = mean - mean(mean)) %>% ungroup() %>%
    select(species, model, rel) %>% pivot_wider(names_from = model, values_from = rel)
  gamma_occ <- imap_dfr(gamma_occ, ~ occ_only(.x) %>% mutate(model = .y))
}

# Species changes with 90% intervals and support, observed vs pooled reduced (16)
occ16 <- read_csv(need(file.path(cmp_dir, "occurrence_change_observed_vs_reduced.csv")), show_col_types = FALSE)
col_or_na <- function(d, name) if (name %in% names(d)) d[[name]] else NA_real_
long_change <- function(d, variants, prefix_mean) imap_dfr(variants, ~ tibble(
  species = d$species, model = .y,
  mean    = d[[paste0(prefix_mean, .x)]],
  q05     = col_or_na(d, paste0("q05_", .x)),
  q95     = col_or_na(d, paste0("q95_", .x)),
  support = col_or_na(d, paste0("support_", .x))))
occ_long <- long_change(occ16, c(observed = "observed", reduced = "reduced_pooled"), "mean_")

glm <- read_csv(need(here("pipeline", "04_robustness", "outputs", "rarefaction", "glm_reduced_10_vs_100.csv")),
                show_col_types = FALSE) %>%
  select(species, set, mean, q025, q975) %>%
  pivot_wider(names_from = set, values_from = c(mean, q025, q975))

beta_chk  <- read_csv(need(file.path(cmp_dir, "checks_beta_design_sensitive.csv")), show_col_types = FALSE) %>%
  group_by(check_name) %>% mutate(main = main - mean(main), check = check - mean(check)) %>% ungroup()
cover_pa  <- read_csv(need(file.path(cmp_dir, "checks_relcover_paired_vs_all.csv")), show_col_types = FALSE)
# Change in relative cover per species with 90% intervals: abundance rows of
# the 04 tables of the main relative-cover model (paired by plot, _rel) and of
# the model fitted to all presences with cover (_relall). Without these tables,
# the posterior means of 16 are used (no intervals).
# The species-change table holding the relative-cover rows is the first of the
# 04 delta tables with rows other than occurrence
read_cover_delta <- function(abu) {
  for (nm in c("delta_species_plotRE.csv", "delta_species_full.csv", "delta_species.csv")) {
    d <- read_res("m1_time_elev_nsp_nos_redpool", nm, abu)
    if (!is.null(d) && "response" %in% names(d) && any(d$response != "occurrence")) {
      cat(sprintf("Relative cover (%s): %s; responses: %s; columns: %s\n", abu, nm,
                  paste(unique(d$response), collapse = ", "), paste(names(d), collapse = ", ")))
      return(filter(d, response != "occurrence"))
    }
  }
  NULL
}
cov_tabs <- list(main = read_cover_delta("rel"), comparison = read_cover_delta("relall"))
if (!any(map_lgl(cov_tabs, is.null))) {
  cov_long <- imap_dfr(cov_tabs, ~ .x %>%
                         transmute(species, model = .y, mean,
                                   q05 = col_or_na(., "q05"), q95 = col_or_na(., "q95"),
                                   support = col_or_na(., "support")))
  # check against 16: correlation of the posterior means of the two variants
  chk <- cov_long %>% select(species, model, mean) %>% pivot_wider(names_from = model, values_from = mean) %>%
    inner_join(cover_pa %>% select(species, rel_paired, rel_all), by = "species")
  cat(sprintf("Relative cover, 04 vs 16: r(main, paired) = %.3f; r(comparison, all) = %.3f\n",
              cor(chk$main, chk$rel_paired, use = "complete.obs"), cor(chk$comparison, chk$rel_all, use = "complete.obs")))
} else {
  message("No 04 table with relative-cover rows found: posterior means of 16 used, without intervals")
  cov_long <- long_change(cover_pa, c(main = "paired", comparison = "all"), "rel_")
}
gamma_cov <- read_csv(need(file.path(cmp_dir, "checks_gamma_relcover_paired_vs_all.csv")), show_col_types = FALSE) %>%
  mutate(model = recode(variant, paired = "main", all_presences = "comparison"))

# -----------------------------------------------------------------------------
# 2. Panel functions
# -----------------------------------------------------------------------------
r_lab <- function(x, y) sprintf("r = %.2f", cor(x, y, use = "complete.obs"))

# Scatter of species changes (one point per species), 1:1 line
scatter_panel <- function(d, x, y, xlab, ylab, title, colour, ranges = NULL) {
  lim <- range(c(d[[x]], d[[y]], if (!is.null(ranges)) unlist(d[ranges])), na.rm = TRUE)
  p <- ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3)
  if (!is.null(ranges)) {
    p <- p +
      geom_segment(aes(x = .data[[ranges[1]]], xend = .data[[ranges[2]]], yend = .data[[y]]),
                   colour = colour, alpha = 0.4, linewidth = 0.8, lineend = "round") +
      geom_segment(aes(xend = .data[[x]], y = .data[[ranges[3]]], yend = .data[[ranges[4]]]),
                   colour = colour, alpha = 0.4, linewidth = 0.8, lineend = "round")
  }
  p + geom_point(size = 2.3, colour = colour) +
    annotate("text", x = lim[1], y = lim[2], label = r_lab(d[[x]], d[[y]]), hjust = 0, vjust = 1, size = 4.3) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = xlab, y = ylab, title = title) +
    theme_mora()
}

# Change per species of two models: species on the y axis (ordered by the
# change in the main model), change on the x axis, 90% intervals when
# available; faded when the posterior support is below sig_level
# Centred change of one model: species change minus the mean change across
# species, so that differences among species are compared independently of
# the change shared by all species. Intervals are shifted by the same mean
# (the uncertainty of the mean across species is not propagated).
add_centred <- function(long, m) {
  cen <- long %>% filter(model == m) %>%
    mutate(shift = mean(mean, na.rm = TRUE),
           across(c(mean, q05, q95), ~ .x - shift),
           model = paste0(m, "_centred")) %>%
    select(-shift)
  bind_rows(long, cen)
}

forest_panel <- function(long, cols, labels, xlab, title, order_by) {
  ord <- long %>% filter(model == order_by, !is.na(mean)) %>% arrange(mean) %>% pull(species)
  ord <- c(setdiff(unique(long$species), ord), ord)   # species missing in the main model at the bottom
  long <- long %>%
    mutate(species = factor(gsub("_", " ", species), levels = gsub("_", " ", ord)),
           model = factor(model, levels = names(labels)),
           supported = coalesce(support >= sig_level | support <= 1 - sig_level, TRUE))
  dodge <- position_dodge(width = 0.75)
  p <- ggplot(long, aes(mean, species, colour = model, alpha = supported)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50")
  if (any(!is.na(long$q05)))
    p <- p + geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.6, lineend = "round", position = dodge)
  p + geom_point(size = 1.7, position = dodge) +
    scale_colour_manual(values = cols, labels = labels, name = NULL) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    labs(x = xlab, y = NULL, title = title) +
    theme_mora() +
    guides(colour = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 8.5, face = "italic", colour = "black"))
}

# Trait effects (Gamma) of two models: points and 90% intervals, round-ended lines
gamma_panel <- function(d, cols, labels, title) {
  d <- d %>% mutate(trait = factor(coalesce(trait_labels[trait], trait), levels = rev(trait_labels)),
                    model = factor(model, levels = names(labels)),
                    supported = support >= sig_level | support <= 1 - sig_level)
  ggplot(d, aes(mean, trait, colour = model, alpha = supported)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.9, lineend = "round",
                   position = position_dodge(width = 0.55)) +
    geom_point(size = 2.6, position = position_dodge(width = 0.55)) +
    scale_colour_manual(values = cols, labels = labels, name = NULL) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
    labs(x = "Trait effect on the response to time (Γ)", y = NULL, title = title) +
    theme_mora() +
    theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())
}

# -----------------------------------------------------------------------------
# 3. Figure S5.3: occurrence model
# -----------------------------------------------------------------------------
empty_panel <- function(title) ggplot() + theme_void() + labs(title = title) +
  annotate("text", x = 0, y = 0, label = "04 results of the observed-data\noccurrence model not found", size = 4.5)
occ_labels <- c(observed = "Observed data", reduced = "Reduced past data (pooled)")
s3a <- forest_panel(occ_long, c(observed = col_cmp, reduced = col_occ), occ_labels,
                    "Change in occurrence probability", "Species change", order_by = "reduced")
s3b <- if (have_obs) gamma_panel(gamma_occ, c(observed = col_cmp, reduced = col_occ), occ_labels, "Trait effects") else
  empty_panel("Trait effects")
s3c <- scatter_panel(glm, "mean_all_100", "mean_hmsc_10",
                     "GLM, mean over 100 reductions", "GLM, mean over 10 reductions",
                     "Species change, 10 vs 100 reductions", col_occ,
                     ranges = c("q025_all_100", "q975_all_100", "q025_hmsc_10", "q975_hmsc_10"))
# Pooled reduced model with and without the 12 sparse species (28): change in
# occurrence probability for the species in both models, centred on the mean
# across species (intervals shifted by the same mean)
sparse <- read_csv(need(here("pipeline", "02_analysis", "02_hmsc", "outputs", "run_sparse12",
                             "sparse12_species_change_comparison.csv")), show_col_types = FALSE) %>%
  mutate(shift_main = mean(mean_main), shift_sparse = mean(mean_sparse12),
         across(c(mean_main, lower_90_main, upper_90_main), ~ .x - shift_main),
         across(c(mean_sparse12, lower_90_sparse12, upper_90_sparse12), ~ .x - shift_sparse))
s3d <- scatter_panel(sparse, "mean_main", "mean_sparse12",
                     "Main model (reduced past data)", "Without 12 sparse species",
                     "Species change without sparse species", col_occ,
                     ranges = c("lower_90_main", "upper_90_main", "lower_90_sparse12", "upper_90_sparse12"))

# A (species list) in the left column; B, C and D stacked on the right
fig3 <- (s3a | (s3b / s3c / s3d)) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A")
f3 <- file.path(fig_dir, paste0(fig_prefix, "figS5_3_occurrence_checks_", tag, ".png"))
ggsave(f3, fig3, width = 14, height = 17, dpi = 300, bg = "white")

# -----------------------------------------------------------------------------
# 4. Figure S5.4: relative-cover model
# -----------------------------------------------------------------------------
cov_labels <- c(main = "Paired by plot (main)", comparison = "All presences with cover")
s4a <- forest_panel(add_centred(cov_long, "main"),
                    c(main = col_abu_raw, main_centred = col_abu, comparison = col_cmp),
                    c(main = "Paired by plot (main)", main_centred = "Paired, centred",
                      comparison = "All presences with cover"),
                    "Change in relative cover", "Species change", order_by = "main")
s4b <- gamma_panel(gamma_cov, c(main = col_abu, comparison = col_cmp), cov_labels, "Trait effects")
patchy <- filter(beta_chk, check_name == "relcover_no_patchy")
s4c <- scatter_panel(patchy, "main", "check",
                     "All species", "Without patchy species",
                     "Species change without patchy species", col_abu)

fig4 <- (s4a | (s4b / s4c)) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A")
f4 <- file.path(fig_dir, paste0(fig_prefix, "figS5_4_cover_checks_", tag, ".png"))
ggsave(f4, fig4, width = 14, height = 13, dpi = 300, bg = "white")
cat("Figures written to", fig_dir, "\n")

# Open the figures (macOS: open; Windows: default viewer)
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
walk(c(f3, f4), open_fig)
