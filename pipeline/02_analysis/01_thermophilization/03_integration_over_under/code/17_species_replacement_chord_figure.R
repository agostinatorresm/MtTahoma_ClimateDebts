# =============================================================================
# 17_species_replacement_chord_figure.R
#
# Pipeline: pipeline/02_analysis/01_thermophilization/03_integration_over_under/code/
# Author: Agostina Torres
#
# Combines the overstory and understory species-replacement chord diagrams
# into one figure for the Supporting Information. The chord diagrams are drawn
# with circlize (base graphics) by
#   01_overstory/code/07_species_replacement_chord.R
#   02_understory/code/08_species_replacement_chord.R
# so this script only assembles the saved images; run those two scripts first.
#
# Panels (net flows, i.e. the dominant direction of replacement per pair):
#   A  overstory, species ordered and coloured by thermal affinity (MAT)
#   B  understory, same
#   C  overstory, species ordered and coloured by moisture affinity (MAP)
#   D  understory, same
# Understory flows are drawn for the species with the largest flows and the
# remaining species grouped by life form (UNDERSTORY_GROUPING); the version
# grouped by life form x thermal class can be used instead.
#
# Inputs : 01_overstory/figures/overstory_species_replacement_chord_{MAT,MAP}_net.png
#          02_understory/figures/understory_replacement_chord_<grouping>_{MAT,MAP}_net.png
# Output : 03_integration_over_under/figures/species_replacement_chord.png
# =============================================================================

library(here)
library(png)    # install.packages("png") if missing
library(grid)

thermo_dir <- function(...) here("pipeline", "02_analysis", "01_thermophilization", ...)
fig_dir <- thermo_dir("03_integration_over_under", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

UNDERSTORY_GROUPING <- "top_species"   # or "lifeform_thermal"

panels <- c(
  A = thermo_dir("01_overstory", "figures", "overstory_species_replacement_chord_MAT_net.png"),
  B = thermo_dir("02_understory", "figures", paste0("understory_replacement_chord_", UNDERSTORY_GROUPING, "_MAT_net.png")),
  C = thermo_dir("01_overstory", "figures", "overstory_species_replacement_chord_MAP_net.png"),
  D = thermo_dir("02_understory", "figures", paste0("understory_replacement_chord_", UNDERSTORY_GROUPING, "_MAP_net.png")))
missing <- panels[!file.exists(panels)]
if (length(missing)) stop("Chord diagrams not found (run the overstory and understory chord scripts first):\n",
                          paste(missing, collapse = "\n"))

# Column titles (layers) and panel letters
col_titles <- c("Overstory", "Understory")

fig_file <- file.path(fig_dir, "species_replacement_chord.png")
png(fig_file, width = 14, height = 14.6, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2, heights = unit(c(0.6, 1, 1), c("in", "null", "null")))))
for (j in 1:2) {
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = j))
  grid.text(col_titles[j], gp = gpar(fontsize = 18, fontface = "bold"))
  popViewport()
}
pos <- list(A = c(2, 1), B = c(2, 2), C = c(3, 1), D = c(3, 2))
for (p in names(panels)) {
  pushViewport(viewport(layout.pos.row = pos[[p]][1], layout.pos.col = pos[[p]][2]))
  grid.raster(readPNG(panels[[p]]), interpolate = TRUE)
  # panel letter at the top right, clear of the colour bar drawn at the top left
  grid.text(p, x = unit(0.97, "npc"), y = unit(0.98, "npc"), just = c("right", "top"),
            gp = gpar(fontsize = 22, fontface = "bold"))
  popViewport()
}
dev.off()
cat("Figure written to", fig_file, "\n")

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)