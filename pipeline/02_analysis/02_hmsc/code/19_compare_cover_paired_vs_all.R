# =============================================================================
# 19_compare_cover_paired_vs_all.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine,
#           or on the local computer after copying the tables listed under Inputs;
#           sections whose inputs are missing are skipped)
# Author: Agostina Torres
#
# Comparison of the two understory relative-cover models that differ only in
# the cover values used:
#   paired : cover kept only where a species has cover in both times of the
#            same plot (main model, m1_time_elev_nsp_nos_rel)
#   all    : every presence with cover; past values set to missing where the
#            present occurrence was recorded only in the 50-m2 subplot
#            (m1_time_elev_nsp_nos_relall)
# Both use the same structure (time x elevation, non-spatial plot level),
# traits, phylogeny and MCMC settings.
#
#   1. MCMC settings and run times (run log)
#   2. Convergence (PSRF, ESS) per parameter group
#   3. Species change in relative cover: both models in one panel
#   4. Trait effects on the relative-cover response to time (Gamma)
#   5. Variance partitioning and phylogenetic signal of relative cover
#
# Inputs : outputs/**/understory_strict_long_run_log.csv                      (02)
#          outputs/understory_strict_long_run_oct09_conv_parameters.csv       (03)
#          outputs/understory_strict_long_m1_time_elev_nsp_nos_redpool_{rel,relall}_res_delta_species_full.csv (04)
#          outputs/understory_strict_long_m1_time_elev_nsp_nos_redpool_{rel,relall}_res_gamma_time.csv    (04)
#          outputs/understory_strict_long_m1_time_elev_nsp_nos_redpool_{rel,relall}_res_{variance_partitioning,phylogenetic_signal}.csv (04)
#          01_data_prep/01_cleaning/outputs/understory_long.csv, plot_info.csv (species order)
# Outputs: outputs/cover_paired_vs_all_settings.csv
#          outputs/cover_paired_vs_all_convergence.csv
#          outputs/cover_paired_vs_all_species_change.csv
#          figures/m1_final/ms1__figS_cover_paired_vs_all_m1_final.png
# =============================================================================

library(here)
library(tidyverse)
library(patchwork)

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
res_dir    <- file.path(module_dir, "outputs")
fig_dir    <- file.path(module_dir, "figures", "m1_final")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

prefix     <- "understory_strict_long_m1_time_elev_nsp_nos_redpool_"
model_ids  <- c(paired = "m1_time_elev_nsp_nos_rel", all = "m1_time_elev_nsp_nos_relall")
model_labs <- c(paired = "Paired by plot (main)", all = "All presences with cover")
# Paired model in the colour used for abundance in the HMSC figures, all-presence model in black
col_model  <- c("Paired by plot (main)" = "#9A3F12", "All presences with cover" = "black")
sig_level  <- 0.90

has_ipsum <- requireNamespace("hrbrthemes", quietly = TRUE)
if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))
theme_mora <- function(...) {
  base <- if (has_ipsum) hrbrthemes::theme_ipsum(...) else theme_minimal(...)
  base +
    theme(panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line  = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text  = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black", size = 11),
          axis.title.y = element_text(colour = "black", size = 11),
          legend.text  = element_text(colour = "black"),
          plot.title   = element_text(colour = "black", size = 11, face = "bold"),
          plot.margin  = margin(12, 10, 10, 10))
}
fold_labels <- function(x) paste0("×", signif(10^x, 2))
find_file   <- function(pattern) list.files(res_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)

# -----------------------------------------------------------------------------
# 1. MCMC settings and run times
# -----------------------------------------------------------------------------
log_files <- find_file("^understory_strict_long_run_log\\.csv$")
if (length(log_files) > 0) {
  settings <- map_dfr(log_files, read_csv, show_col_types = FALSE) %>%
    filter(model_id %in% model_ids) %>%
    group_by(model_id) %>% slice_tail(n = 1) %>% ungroup() %>%
    mutate(model = names(model_ids)[match(model_id, model_ids)]) %>%
    select(model, model_id, samples, thin, transient, n_chains, elapsed_min, finished_at,
           any_of(c("machine", "r_version", "hmsc_version")))
  print(settings, width = Inf)
  write_csv(settings, file.path(res_dir, "cover_paired_vs_all_settings.csv"))
} else message("Run log not found: section 1 skipped")

# -----------------------------------------------------------------------------
# 2. Convergence
# -----------------------------------------------------------------------------
conv_file <- file.path(res_dir, "understory_strict_long_run_oct09_conv_parameters.csv")
if (file.exists(conv_file)) {
  conv <- read_csv(conv_file, show_col_types = FALSE) %>%
    filter(model_id %in% model_ids, !is.na(psrf)) %>%
    mutate(model = names(model_ids)[match(model_id, model_ids)]) %>%
    group_by(model, group) %>%
    summarise(n_parameters = n(), psrf_median = round(median(psrf), 3), psrf_max = round(max(psrf), 3),
              prop_psrf_above_1.1 = round(mean(psrf > 1.1), 3),
              ess_median = round(median(ess)), ess_min = round(min(ess)), .groups = "drop") %>%
    arrange(group, model)
  print(conv, n = Inf, width = Inf)
  write_csv(conv, file.path(res_dir, "cover_paired_vs_all_convergence.csv"))
} else message("Convergence table not found: section 2 skipped")

# -----------------------------------------------------------------------------
# 3. Species change in relative cover
# -----------------------------------------------------------------------------
read_model <- function(name) {
  imap_dfr(model_ids, function(id, m) {
    f <- file.path(res_dir, paste0(prefix, sub("^m1_time_elev_nsp_nos_", "", id), "_res_", name))
    if (!file.exists(f)) stop("Not found (run 04 for this model): ", basename(f))
    read_csv(f, show_col_types = FALSE) %>% mutate(model = model_labs[[m]])
  })
}
# Change predicted by the full model (fixed effects and plot-level random effect),
# averaged over plots, as in Supporting Information S5 (16_compare_checks_oct09.R)
delta <- read_model("delta_species_full.csv") %>% filter(response == "log10_abundance")
gamma <- read_model("gamma_time.csv")   %>% filter(response == "abundance")

# Species ordered by the mean elevation of the plots where they occur (as in 05)
cleaning_dir <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
plot_info <- read_csv(cleaning_dir("outputs", "plot_info.csv"), show_col_types = FALSE)
plot_key  <- intersect(c("plot_id", "plot", "Plot", "PlotID"), names(plot_info))[1]
elev_col  <- intersect(c("elevation", "elevation_m", "elev", "elev_m"), names(plot_info))[1]
plot_elev <- plot_info %>%
  transmute(plot_id = as.character(.data[[plot_key]]), elevation = as.numeric(.data[[elev_col]])) %>%
  filter(!is.na(elevation)) %>% distinct(plot_id, .keep_all = TRUE)
sp_order <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE) %>%
  filter(occurrence > 0) %>%
  mutate(plot_id = as.character(plot_id), species = str_replace_all(species, " ", "_")) %>%
  left_join(plot_elev, by = "plot_id") %>%
  group_by(species) %>% summarise(mean_elevation = mean(elevation, na.rm = TRUE), .groups = "drop") %>%
  filter(species %in% delta$species) %>%
  arrange(mean_elevation) %>% pull(species)

# Agreement between models (species in both)
wide <- delta %>% select(species, model, mean, direction) %>%
  pivot_wider(names_from = model, values_from = c(mean, direction))
m_p <- paste0("mean_", model_labs[["paired"]]); m_a <- paste0("mean_", model_labs[["all"]])
common <- wide %>% filter(!is.na(.data[[m_p]]), !is.na(.data[[m_a]]))
cat("\nSpecies: paired", sum(!is.na(wide[[m_p]])), "| all", sum(!is.na(wide[[m_a]])),
    "| in both", nrow(common), "\n")
cat("Correlation of species change (log10 fold, full model):", round(cor(common[[m_p]], common[[m_a]]), 3), "\n")
print(delta %>% count(model, direction) %>% pivot_wider(names_from = direction, values_from = n, values_fill = 0))
write_csv(wide, file.path(res_dir, "cover_paired_vs_all_species_change.csv"))

# Panels: the two models side by side within each species row (manual offset)
offset <- c("Paired by plot (main)" = 0.18, "All presences with cover" = -0.18)
supp_lab <- c(paste0("Supported (≥ ", sig_level, ")"), "Not supported")
prep <- function(d, levels) {
  d %>% mutate(row = match(.data$key, levels) + offset[model],
               supported = factor(direction != "not supported", levels = c(TRUE, FALSE), labels = supp_lab),
               model = factor(model, levels = names(col_model)))
}
range_layers <- function() list(
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.5),
  geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 1.2),
  geom_point(size = 1.6),
  scale_colour_manual(values = col_model, name = NULL, drop = FALSE),
  scale_alpha_manual(values = c(1, 0.3), name = NULL, drop = FALSE)
)

d_sp <- prep(delta %>% mutate(key = species), sp_order)
p_sp <- ggplot(d_sp, aes(y = row, x = mean, colour = model, alpha = supported)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  range_layers() +
  scale_x_continuous(labels = fold_labels) +
  scale_y_continuous(breaks = seq_along(sp_order), labels = str_replace_all(sp_order, "_", " "),
                     expand = expansion(add = 0.6)) +
  labs(x = "Change in relative cover when present", y = "Species (low to high mean elevation)",
       title = "A. Species change") +
  theme_mora() +
  theme(axis.text.y = element_text(face = "italic", colour = "black", size = 7), legend.position = "top")

trait_labels <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity",
                  SLA = "SLA", seed_mass = "Seed mass", woodinesswoody = "Woody (vs non-woody)")
gamma <- gamma %>% filter(trait %in% names(trait_labels))
tr_levels <- rev(names(trait_labels))
d_g <- prep(gamma %>% mutate(key = trait), tr_levels)
p_g <- ggplot(d_g, aes(y = row, x = mean, colour = model, alpha = supported)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  range_layers() +
  scale_y_continuous(breaks = seq_along(tr_levels), labels = trait_labels[tr_levels],
                     expand = expansion(add = 0.6)) +
  labs(x = "Trait effect on change in relative cover (Γ)", y = NULL, title = "B. Trait effects") +
  theme_mora() +
  theme(legend.position = "top")

# -----------------------------------------------------------------------------
# 5. Variance partitioning and phylogenetic signal (relative cover)
# -----------------------------------------------------------------------------
vp <- read_model("variance_partitioning.csv") %>%
  filter(response == "abundance") %>%
  group_by(model, component) %>%
  summarise(proportion = mean(proportion), .groups = "drop") %>%
  mutate(model = factor(model, levels = names(col_model)),
         component = recode(component, "Random: plot" = "Plot"))
print(vp %>% mutate(proportion = round(proportion, 3)) %>% pivot_wider(names_from = model, values_from = proportion))
p_vp <- ggplot(vp, aes(x = component, y = proportion, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = col_model, name = NULL, guide = "none") +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(x = NULL, y = "Explained variance", title = "C. Variance partitioning") +
  theme_mora()

rho <- read_model("phylogenetic_signal.csv") %>%
  filter(response == "abundance") %>%
  mutate(model = factor(model, levels = names(col_model)))
print(rho %>% select(model, mean, q05, q95))
p_rho <- ggplot(rho, aes(x = model, y = mean, ymin = q05, ymax = q95, colour = model)) +
  geom_pointrange(size = 0.6, linewidth = 0.8) +
  scale_colour_manual(values = col_model, guide = "none") +
  scale_x_discrete(labels = c("Paired", "All presences")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Phylogenetic signal (\u03c1)", title = "D. Phylogenetic signal") +
  theme_mora()

fig <- ((p_sp | (p_g / (p_vp | p_rho) + plot_layout(heights = c(1.4, 1))))) +
  plot_layout(widths = c(1.2, 1), guides = "collect") &
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "ms1__figS_cover_paired_vs_all_m1_final.png"), fig,
       width = 13, height = 11, dpi = 300, bg = "white")
cat("Figure written to", fig_dir, "\n")
