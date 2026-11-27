# =============================================================================
# 10_raw_cover_diagnosis.R
#
# Pipeline: pipeline/04_robustness/code/ (run on the local computer)
# Author: Agostina Torres
#
# Descriptive diagnosis of raw understory cover in the past and present, for
# the 58 HMSC species, to show how the smaller present sampling area affects
# the cover values available for comparison:
#   A. Distribution of raw cover (log10 %) by type of record:
#        paired       - species with cover at both times in the same plot
#                       (past and present values shown separately);
#        past only    - cover in the past, no cover in the present
#                       (absent, or recorded only in the 50-m2 census);
#        present only - cover in the present, absent in the past.
#      Present occurrences recorded only in the 50-m2 census have no cover and
#      are counted in the panel annotation.
#   B. Number of species with cover per plot in the past and present.
#   C. Change in raw cover (log10 present/past) of the paired records, for
#      patchy and not patchy species (01b). Both groups share any overall
#      difference between the past and present cover estimates; a different
#      change in patchy species indicates a bias linked to the small present
#      sampling area. The two groups are compared with a Wilcoxon rank-sum test
#      on the median change of each species, because records of the same
#      species are not independent; the test on all records is also printed.
# Raw cover is the percent cover recorded in the past (whole-plot visual
# estimate) and the mean across the sampled quadrats in the present.
#
# Inputs : pipeline/01_data_prep/01_cleaning/outputs/understory_long.csv
#          pipeline/02_analysis/02_hmsc/data_clean/understory_strict_hmsc_inputs.rds (species set)
#          pipeline/04_robustness/outputs/accumulation_curve/accumulation_curve_summary.csv
#            (patchy species, from 01b)
# Outputs: pipeline/04_robustness/outputs/raw_cover/raw_cover_by_record_type.csv
#          pipeline/04_robustness/outputs/raw_cover/species_with_cover_per_plot.csv
#          pipeline/04_robustness/outputs/raw_cover/paired_change_by_patchiness.csv
#          pipeline/04_robustness/tmp_files/panels/panels_10.rds (panels A-C of Figure S3.3,
#          assembled in 08_diagnosis_figures.R)
# =============================================================================

library(here)
library(tidyverse)
library(patchwork)

out_dir <- here("pipeline", "04_robustness", "outputs", "raw_cover")
panel_dir <- here("pipeline", "04_robustness", "tmp_files", "panels")
walk(c(out_dir, panel_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

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
# Past and present in the light and dark orange of the HMSC figures
col_time <- c(Past = "#F2A65A", Present = "#9A3F12")

# -----------------------------------------------------------------------------
# 1. Inputs: understory records of the HMSC species
# -----------------------------------------------------------------------------
veg <- read_csv(here("pipeline", "01_data_prep", "01_cleaning", "outputs", "understory_long.csv"),
                show_col_types = FALSE)
inp <- readRDS(here("pipeline", "02_analysis", "02_hmsc", "data_clean", "understory_strict_hmsc_inputs.rds"))
hmsc_species <- str_replace_all(rownames(inp$TrData), "_", " ")

pres_rec <- veg %>%
  filter(species %in% hmsc_species, occurrence == 1) %>%
  mutate(time = if_else(period == "past", "Past", "Present"),
         has_cover = !is.na(cover) & cover > 0)

# -----------------------------------------------------------------------------
# 2. Record types (plot x species)
# -----------------------------------------------------------------------------
by_pair <- pres_rec %>%
  group_by(plot_id, species) %>%
  summarise(cover_past    = cover[time == "Past" & has_cover][1],
            cover_present = cover[time == "Present" & has_cover][1],
            subplot_only  = any(time == "Present" & !has_cover), .groups = "drop")

records <- bind_rows(
  by_pair %>% filter(!is.na(cover_past), !is.na(cover_present)) %>%
    transmute(type = "Paired", time = "Past", cover = cover_past),
  by_pair %>% filter(!is.na(cover_past), !is.na(cover_present)) %>%
    transmute(type = "Paired", time = "Present", cover = cover_present),
  by_pair %>% filter(!is.na(cover_past), is.na(cover_present)) %>%
    transmute(type = "Past only", time = "Past", cover = cover_past),
  by_pair %>% filter(is.na(cover_past), !is.na(cover_present)) %>%
    transmute(type = "Present only", time = "Present", cover = cover_present)
) %>%
  mutate(type = factor(type, levels = c("Paired", "Past only", "Present only")))

summary_type <- records %>%
  group_by(type, time) %>%
  summarise(n = n(), median_cover = round(median(cover), 2),
            mean_log10_cover = round(mean(log10(cover)), 2), .groups = "drop")
n_subplot_only <- sum(by_pair$subplot_only & is.na(by_pair$cover_present))
cat("\nRaw cover by record type (HMSC species)\n"); print(summary_type)
cat("Present occurrences recorded only in the 50-m2 census (no cover):", n_subplot_only, "\n")
write_csv(summary_type, file.path(out_dir, "raw_cover_by_record_type.csv"))

# -----------------------------------------------------------------------------
# 3. Species with cover per plot
# -----------------------------------------------------------------------------
sp_per_plot <- pres_rec %>%
  filter(has_cover) %>%
  count(plot_id, time, name = "n_species") %>%
  complete(plot_id, time, fill = list(n_species = 0))
summary_plot <- sp_per_plot %>%
  group_by(time) %>%
  summarise(mean = round(mean(n_species), 1), median = median(n_species), .groups = "drop")
cat("\nSpecies with cover per plot\n"); print(summary_plot)
write_csv(sp_per_plot, file.path(out_dir, "species_with_cover_per_plot.csv"))

# -----------------------------------------------------------------------------
# 3b. Change in raw cover of paired records: patchy vs not patchy species
# -----------------------------------------------------------------------------
patchy_tab <- read_csv(here("pipeline", "04_robustness", "outputs", "accumulation_curve",
                            "accumulation_curve_summary.csv"), show_col_types = FALSE) %>%
  mutate(species = str_replace_all(species, "_", " "))
patchy_sp <- patchy_tab$species[patchy_tab$patchy]

paired_change <- by_pair %>%
  filter(!is.na(cover_past), !is.na(cover_present)) %>%
  mutate(change = log10(cover_present / cover_past),
         flag = factor(if_else(species %in% patchy_sp, "Patchy", "Not patchy"),
                       levels = c("Not patchy", "Patchy")))
summary_change <- paired_change %>%
  group_by(flag) %>%
  summarise(n_records = n(), n_species = n_distinct(species),
            median_change = round(median(change), 2),
            prop_increase = round(mean(change > 0), 2), .groups = "drop")
# Test on species medians (one value per species); the test on records is
# printed for reference only, since records of one species are not independent
species_change <- paired_change %>%
  group_by(species, flag) %>% summarise(median_change = median(change), .groups = "drop")
wt_sp  <- wilcox.test(median_change ~ flag, data = species_change, exact = FALSE)
wt_rec <- wilcox.test(change ~ flag, data = paired_change, exact = FALSE)
cat("\nChange in raw cover of paired records (log10 present/past)\n"); print(summary_change)
cat("Wilcoxon rank-sum test on species medians, patchy vs not patchy: W =", wt_sp$statistic,
    "| p =", signif(wt_sp$p.value, 3), "| n species =", table(species_change$flag), "\n")
cat("Same test on all records (not independent; reference only): p =", signif(wt_rec$p.value, 3), "\n")
write_csv(summary_change, file.path(out_dir, "paired_change_by_patchiness.csv"))

# -----------------------------------------------------------------------------
# 4. Panels for Figure S3.3
# -----------------------------------------------------------------------------
labs_n <- summary_type %>%
  mutate(label = paste0("n = ", n), y = max(records$cover) * 1.6)
cover_breaks <- c(0.1, 1, 10, 100)

p_a <- ggplot(records, aes(x = type, y = cover, fill = time)) +
  geom_violin(position = position_dodge(width = 0.8), width = 0.75, colour = NA, alpha = 0.85, scale = "width") +
  geom_boxplot(aes(group = interaction(type, time)),
               position = position_dodge(width = 0.8), width = 0.12, outlier.shape = NA,
               fill = "white", colour = "black", linewidth = 0.35) +
  geom_text(data = labs_n, aes(x = type, y = y, label = label, group = time),
            position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  scale_y_log10(breaks = cover_breaks, labels = cover_breaks) +
  scale_fill_manual(values = col_time, name = NULL) +
  labs(x = NULL, y = "Raw cover (%, log scale)", title = "Raw cover by type of record",
       caption = paste0("Present occurrences recorded only in the 50-m² census (no cover): ", n_subplot_only)) +
  theme_mora() +
  theme(legend.position = "top", plot.caption = element_text(hjust = 0, size = 9))

p_b <- ggplot(sp_per_plot, aes(x = time, y = n_species, fill = time)) +
  geom_violin(colour = NA, alpha = 0.85, scale = "width", width = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "black", linewidth = 0.35) +
  scale_fill_manual(values = col_time, guide = "none") +
  labs(x = NULL, y = "Species with cover per plot", title = "Species with cover per plot") +
  theme_mora()

labs_c <- summary_change %>%
  mutate(label = paste0("n = ", n_records, " (", n_species, " spp.)"),
         y = max(paired_change$change) * 1.15)
p_c <- ggplot(paired_change, aes(x = flag, y = change, fill = flag)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_violin(colour = NA, alpha = 0.85, scale = "width", width = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "black", linewidth = 0.35) +
  geom_text(data = labs_c, aes(x = flag, y = y, label = label), size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Not patchy" = "grey55", "Patchy" = "#9A3F12"), guide = "none") +
  labs(x = NULL, y = "Change in raw cover\n(log10 present / past)",
       title = "Paired records: change in cover",
       caption = paste0("Wilcoxon rank-sum test on species medians: p = ", signif(wt_sp$p.value, 2))) +
  theme_mora() +
  theme(plot.caption = element_text(hjust = 0, size = 9))

saveRDS(list(raw_cover = p_a, species_cover = p_b, patchy_change = p_c), file.path(panel_dir, "panels_10.rds"))
cat("\nPanels written to", panel_dir, "\n")