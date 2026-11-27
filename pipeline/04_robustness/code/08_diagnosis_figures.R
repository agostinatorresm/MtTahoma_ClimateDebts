# =============================================================================
# 08_diagnosis_figures.R
#
# Pipeline: pipeline/04_robustness/code/ (robustness module: past vs present
#           understory sampling design; run order 00 -> 09, or 00_run_all_robustness.R)
# Author: Agostina Torres
#
# Combined figures for the diagnosis of the past vs present understory
# sampling design (Supporting Information S3), assembled from the panels saved by the
# diagnostic scripts:
#   Figure S3.1 (occurrence): (A) richness per plot, past and present (03);
#     (B) past richness by plot area (03); (C) species accumulation over the
#     present sampling units (01a)
#   Figure S3.3 (cover): (A) raw cover by type of record (10); (B) species
#     with cover per plot (10); (C) change in raw cover of paired records for
#     patchy and not patchy species (10); (D) stability of relative cover with
#     the number of quadrats (01a); (E) bootstrap tests of the change in
#     relative cover (02); (F) distribution of mean quadrat occupancy with the
#     patchy threshold (01b)
#
# Inputs : pipeline/04_robustness/tmp_files/panels/panels_01a.rds, panels_01b.rds,
#          panels_02.rds, panels_03.rds, panels_10.rds
# Outputs: pipeline/04_robustness/figures/
#            figS3_1_occurrence_diagnosis.png
#            figS3_3_cover_diagnosis.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

panel_dir <- here("pipeline", "04_robustness", "tmp_files", "panels")
sup_dir <- here("pipeline", "04_robustness", "figures")   # Figures S3.1 and S3.3
dir.create(sup_dir, recursive = TRUE, showWarnings = FALSE)
if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))

read_panels <- function(id) {
  f <- file.path(panel_dir, paste0("panels_", id, ".rds"))
  if (!file.exists(f)) stop("Panels not found: ", f, " (run script ", id, " first)")
  readRDS(f)
}
p01a <- read_panels("01a"); p01b <- read_panels("01b"); p02 <- read_panels("02"); p03 <- read_panels("03")
p10 <- read_panels("10")

# Figure S3.1: occurrence
sr1 <- (p03$richness + p03$richness_area + p01a$accumulation) +
  plot_layout(nrow = 1, widths = c(1, 1, 1.3)) +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
ggsave(file.path(sup_dir, "figS3_1_occurrence_diagnosis.png"), sr1,
       width = 15, height = 5, dpi = 300, bg = "white")

# Figure S3.3: cover
# top row: raw cover (10); bottom row: precision of present cover (01a, 02, 01b)
sr2 <- (p10$raw_cover + p10$species_cover + p10$patchy_change) /
  (p01a$cover_stability + p02$bootstrap + p01b$occupancy) +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
ggsave(file.path(sup_dir, "figS3_3_cover_diagnosis.png"), sr2,
       width = 15, height = 10, dpi = 300, bg = "white")

cat("Figures written to", sup_dir, "\n")
