# ------------------------------------------------------------------------------
# 08_affinity_correlation_figure.R
# Module: pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: relationship between species' thermal (MAT) and moisture (MAP)
# affinities, for the SI (S2.2). CTI and CMI are interpreted as two separate
# dimensions of community change; this is valid only if the two affinities are
# not strongly correlated across species.
#
# Species: all species with climatic affinity, and the overstory and understory
# shared pools used for CTI and CMI. A species in both pools appears in both
# panels. Pearson and Spearman correlations are reported for each set.
#
# Output:
#   - outputs/affinity_MAT_MAP_correlation.csv
#   - figures/affinity_MAT_MAP_correlation.png
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(hrbrthemes)

# 0a. Module directories --------------------------------------------------------
this_module_dir  <- function(...) here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa", ...)
species_pool_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)
dir.create(this_module_dir("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("figures"), recursive = TRUE, showWarnings = FALSE)

# 1. Load ------------------------------------------------------------------------
affinity <- read_csv(this_module_dir("data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  select(species, MAT_mean, MAP_mean) %>%
  filter(!is.na(MAT_mean), !is.na(MAP_mean))

pool <- function(layer) {
  read_csv(species_pool_dir("outputs", paste0(layer, "_species_pool_broad_full_clim_cov.csv")),
           show_col_types = FALSE) %>% pull(species)
}

sets <- bind_rows(
  affinity %>% mutate(set = "All species"),
  affinity %>% filter(species %in% pool("overstory"))  %>% mutate(set = "Overstory"),
  affinity %>% filter(species %in% pool("understory")) %>% mutate(set = "Understory")
) %>%
  mutate(set = factor(set, levels = c("All species", "Overstory", "Understory")))

# 2. Correlations ----------------------------------------------------------------
cors <- sets %>%
  group_by(set) %>%
  summarise(n = n(),
            pearson_r  = cor(MAT_mean, MAP_mean),
            pearson_p  = cor.test(MAT_mean, MAP_mean)$p.value,
            spearman_rho = cor(MAT_mean, MAP_mean, method = "spearman"),
            spearman_p = suppressWarnings(cor.test(MAT_mean, MAP_mean, method = "spearman")$p.value),
            .groups = "drop")
print(cors)
write_csv(cors, this_module_dir("outputs", "affinity_MAT_MAP_correlation.csv"))

# Panel titles carry n and Pearson r
strip_lab <- cors %>%
  transmute(set, title = sprintf("%s (n = %d; r = %.2f, p = %.2f)", set, n, pearson_r, pearson_p)) %>%
  deframe()

# 3. Figure ------------------------------------------------------------------------
# Darker shade = overstory, lighter = understory, as in the other figures
set_cols <- c("All species" = "grey45", "Overstory" = "#3B3B3B", "Understory" = "#A6A6A6")

p <- ggplot(sets, aes(MAT_mean, MAP_mean)) +
  geom_point(aes(colour = set), size = 1.8, alpha = 0.7, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.5, fill = "grey85") +
  scale_colour_manual(values = set_cols) +
  facet_wrap(~ set, nrow = 1, labeller = as_labeller(strip_lab)) +
  labs(x = "Thermal affinity (MAT, °C)", y = "Moisture affinity (MAP, mm)") +
  theme_ipsum(base_size = 11, axis_title_size = 11, strip_text_size = 11) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.3),
        panel.spacing = unit(1.2, "lines"),
        plot.margin = margin(10, 10, 10, 10))

f_fig <- this_module_dir("figures", "affinity_MAT_MAP_correlation.png")
ggsave(f_fig, p, width = 10, height = 3.8, dpi = 300, bg = "white")

# 4. Open the figure ----------------------------------------------------------------
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(f_fig)
