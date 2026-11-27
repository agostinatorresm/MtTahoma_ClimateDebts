# ------------------------------------------------------------------------------
# 03_climate_pca_layout.R
# Author: Agostina Torres
# Combine:
#   A) Joint past/present climate PCA
#   B) PCA of plot-specific climatic changes
# Only the plots used in the analyses (324) are included.
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(hrbrthemes)
library(ggrepel)
library(patchwork)
library(grid)

this_module_dir <- function(...) {
  here(
    "pipeline",
    "01_data_prep",
    "03_climate",
    "02_explore_climate",
    ...
  )
}

dir.create(
  this_module_dir("figures"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 1. THEME
# ==============================================================================

climate_theme_large <- theme_ipsum(
  grid = FALSE,
  base_size = 18
) +
  theme(
    axis.line.x = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.line.y = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    
    axis.title.x = element_text(
      size = 22,
      colour = "black",
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 22,
      colour = "black",
      face = "bold",
      margin = margin(r = 10)
    ),
    
    axis.text = element_text(
      size = 17,
      colour = "black"
    ),
    
    legend.title = element_text(
      size = 17,
      colour = "black",
      face = "bold"
    ),
    legend.text = element_text(
      size = 16,
      colour = "black"
    ),
    
    plot.tag = element_text(
      size = 27,
      face = "bold",
      colour = "black",
      hjust = 0,
      margin = margin(b = 6)
    ),
    
    # Keep A/B outside plotting region
    plot.tag.position = "topleft",
    
    plot.margin = margin(
      t = 22,
      r = 15,
      b = 12,
      l = 12
    )
  )


# ==============================================================================
# 2. LOAD DATA
# ==============================================================================

climate_past_present_wide <- read_csv(
  this_module_dir(
    "data_clean",
    "climate_past_present_wide.csv"
  ),
  show_col_types = FALSE
)

# Only the plots used in the analyses (the 324 plots of the overstory
# climatic-affinity analyses, which contain the 316 understory plots)
analysis_plots <- read_csv(
  here("pipeline", "02_analysis", "01_thermophilization", "01_overstory", "outputs",
       "overstory_thermo_xero_rates.csv"),
  show_col_types = FALSE
) %>%
  distinct(plot_id = as.character(plot_id)) %>%
  pull(plot_id)

climate_past_present_wide <- climate_past_present_wide %>%
  mutate(plot_id = as.character(plot_id)) %>%
  filter(plot_id %in% analysis_plots)

cat("Plots used:", nrow(climate_past_present_wide), "of", length(analysis_plots), "analysis plots\n")


# ==============================================================================
# 3. PANEL A — PAST VS PRESENT CLIMATE
# ==============================================================================

climate_stems <- c(
  "MAT",
  "MAT_min",
  "MAP",
  "PAS",
  "temp_DJF",
  "temp_MAM",
  "temp_JJA",
  "temp_SON",
  "precip_DJF",
  "precip_MAM",
  "precip_JJA",
  "precip_SON"
)


# Past + present in long format
pp_long <- bind_rows(
  
  climate_past_present_wide %>%
    transmute(
      plot_id,
      window = "Past",
      across(
        all_of(paste0(climate_stems, "_past")),
        .names = "{str_remove(.col, '_past$')}"
      )
    ),
  
  climate_past_present_wide %>%
    transmute(
      plot_id,
      window = "Present",
      across(
        all_of(paste0(climate_stems, "_present")),
        .names = "{str_remove(.col, '_present$')}"
      )
    )
  
) %>%
  mutate(
    window = factor(
      window,
      levels = c("Past", "Present")
    )
  ) %>%
  drop_na()


# PCA
pca_pp <- prcomp(
  pp_long %>%
    select(all_of(climate_stems)),
  center = TRUE,
  scale. = TRUE
)


# Variance explained
var_pp <- 100 *
  pca_pp$sdev^2 /
  sum(pca_pp$sdev^2)


# Scores
scores_pp <- as_tibble(
  pca_pp$x
) %>%
  select(
    PC1,
    PC2
  ) %>%
  bind_cols(
    pp_long %>%
      select(
        plot_id,
        window
      )
  )


# Loadings
loadings_pp <- as_tibble(
  pca_pp$rotation,
  rownames = "variable"
)


# Scale arrows
arrow_scale_pp <- 0.65 *
  max(
    abs(scores_pp$PC1),
    abs(scores_pp$PC2)
  ) /
  max(
    abs(loadings_pp$PC1),
    abs(loadings_pp$PC2)
  )


loadings_pp_scaled <- loadings_pp %>%
  mutate(
    PC1_arrow = PC1 * arrow_scale_pp,
    PC2_arrow = PC2 * arrow_scale_pp,
    
    # Clean labels for plotting only
    label = recode(
      variable,
      
      "MAT"        = "MAT",
      "MAT_min"    = "MAT min",
      "MAP"        = "MAP",
      "PAS"        = "PAS",
      
      "temp_DJF"   = "T DJF",
      "temp_MAM"   = "T MAM",
      "temp_JJA"   = "T JJA",
      "temp_SON"   = "T SON",
      
      "precip_DJF" = "P DJF",
      "precip_MAM" = "P MAM",
      "precip_JJA" = "P JJA",
      "precip_SON" = "P SON"
    )
  )


# Centroids
centroids_pp <- scores_pp %>%
  group_by(window) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )


past_cent <- centroids_pp %>%
  filter(window == "Past")

present_cent <- centroids_pp %>%
  filter(window == "Present")


# Viridis scales
window_color_scale <- scale_color_viridis_d(
  option = "viridis",
  begin = 0.15,
  end = 0.75,
  name = NULL
)

window_fill_scale <- scale_fill_viridis_d(
  option = "viridis",
  begin = 0.15,
  end = 0.75,
  name = NULL
)


# Plot A
p_A <- ggplot(
  scores_pp,
  aes(
    x = PC1,
    y = PC2,
    colour = window,
    fill = window
  )
) +
  
  geom_point(
    alpha = 0.35,
    size = 2.2
  ) +
  
  stat_ellipse(
    geom = "polygon",
    alpha = 0.12,
    linewidth = 0.8,
    level = 0.68
  ) +
  
  # Loading arrows
  geom_segment(
    data = loadings_pp_scaled,
    aes(
      x = 0,
      y = 0,
      xend = PC1_arrow,
      yend = PC2_arrow
    ),
    inherit.aes = FALSE,
    arrow = arrow(
      length = unit(0.22, "cm")
    ),
    colour = "black",
    linewidth = 0.6
  ) +
  
  # Clean loading labels
  geom_text_repel(
    data = loadings_pp_scaled,
    aes(
      x = PC1_arrow,
      y = PC2_arrow,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5.2,
    colour = "black",
    max.overlaps = Inf,
    box.padding = 0.50,
    point.padding = 0.30,
    min.segment.length = 0,
    seed = 123
  ) +
  
  # Past -> present centroid displacement
  geom_segment(
    data = tibble(
      x = past_cent$PC1,
      y = past_cent$PC2,
      xend = present_cent$PC1,
      yend = present_cent$PC2
    ),
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    arrow = arrow(
      length = unit(0.28, "cm")
    ),
    linewidth = 1,
    colour = "black"
  ) +
  
  # Period centroids
  geom_point(
    data = centroids_pp,
    aes(
      x = PC1,
      y = PC2
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 4,
    stroke = 0.9,
    fill = "white",
    colour = "black"
  ) +
  
  window_color_scale +
  window_fill_scale +
  
  labs(
    x = paste0(
      "PC1 (",
      round(var_pp[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_pp[2], 1),
      "%)"
    )
  ) +
  
  climate_theme_large +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  
  guides(
    colour = guide_legend(
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE
    ),
    fill = guide_legend(
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE
    )
  )


# ==============================================================================
# 4. PANEL B — PCA OF CLIMATIC CHANGE
# ==============================================================================

delta_vars <- paste0(
  "delta_",
  c(
    "MAT",
    "MAT_min",
    "MAP",
    "PAS",
    "temp_DJF",
    "temp_MAM",
    "temp_JJA",
    "temp_SON",
    "precip_DJF",
    "precip_MAM",
    "precip_JJA",
    "precip_SON"
  )
)


# Prepare data
pca_input_delta <- climate_past_present_wide %>%
  select(
    plot_id,
    elevation,
    all_of(delta_vars)
  ) %>%
  drop_na()


pca_vars_delta <- c(
  delta_vars,
  "elevation"
)


# PCA
pca_fit_delta <- prcomp(
  pca_input_delta %>%
    select(all_of(pca_vars_delta)),
  center = TRUE,
  scale. = TRUE
)


# Variance explained
var_delta <- 100 *
  pca_fit_delta$sdev^2 /
  sum(pca_fit_delta$sdev^2)

# Variance explained by the first two axes of each PCA (reported in the Results)
cat("\nPanel A, past/present climate PCA -- PC1:", round(var_pp[1], 1), "% | PC2:", round(var_pp[2], 1),
    "% | PC1 + PC2:", round(sum(var_pp[1:2]), 1), "%\n")
cat("Panel B, PCA of climatic changes -- PC1:", round(var_delta[1], 1), "% | PC2:", round(var_delta[2], 1),
    "% | PC1 + PC2:", round(sum(var_delta[1:2]), 1), "%\n")


# Scores
scores_delta <- as_tibble(
  pca_fit_delta$x
) %>%
  select(
    PC1,
    PC2
  ) %>%
  bind_cols(
    pca_input_delta %>%
      select(
        plot_id,
        elevation
      )
  )


# Loadings
loadings_delta <- as_tibble(
  pca_fit_delta$rotation,
  rownames = "variable"
)


# Scale arrows
arrow_scale_delta <- 0.70 *
  max(
    abs(scores_delta$PC1),
    abs(scores_delta$PC2)
  ) /
  max(
    abs(loadings_delta$PC1),
    abs(loadings_delta$PC2)
  )


# Clean labels for B
loadings_delta_scaled <- loadings_delta %>%
  mutate(
    PC1_arrow = PC1 * arrow_scale_delta,
    PC2_arrow = PC2 * arrow_scale_delta,
    
    label = recode(
      variable,
      
      "delta_MAT"        = "ΔMAT",
      "delta_MAT_min"    = "ΔMAT min",
      "delta_MAP"        = "ΔMAP",
      "delta_PAS"        = "ΔPAS",
      
      "delta_temp_DJF"   = "ΔT DJF",
      "delta_temp_MAM"   = "ΔT MAM",
      "delta_temp_JJA"   = "ΔT JJA",
      "delta_temp_SON"   = "ΔT SON",
      
      "delta_precip_DJF" = "ΔP DJF",
      "delta_precip_MAM" = "ΔP MAM",
      "delta_precip_JJA" = "ΔP JJA",
      "delta_precip_SON" = "ΔP SON",
      
      "elevation"        = "Elevation"
    )
  )


# Plot B
p_B <- ggplot() +
  
  geom_point(
    data = scores_delta,
    aes(
      x = PC1,
      y = PC2,
      colour = elevation
    ),
    size = 2.2,
    alpha = 0.80
  ) +
  
  scale_color_viridis_c(
    option = "viridis",
    direction = -1,
    name = "Elevation (m)"
  ) +
  
  # Loading arrows
  geom_segment(
    data = loadings_delta_scaled,
    aes(
      x = 0,
      y = 0,
      xend = PC1_arrow,
      yend = PC2_arrow
    ),
    arrow = arrow(
      length = unit(0.22, "cm")
    ),
    colour = "black",
    linewidth = 0.6
  ) +
  
  # Clean loading labels
  geom_text_repel(
    data = loadings_delta_scaled,
    aes(
      x = PC1_arrow,
      y = PC2_arrow,
      label = label
    ),
    size = 5.2,
    colour = "black",
    max.overlaps = Inf,
    box.padding = 0.50,
    point.padding = 0.30,
    min.segment.length = 0,
    seed = 123
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(var_delta[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(var_delta[2], 1),
      "%)"
    )
  ) +
  
  climate_theme_large +
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.margin = margin(
      l = 10,
      r = 5
    )
  ) +
  
  guides(
    colour = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(6, "cm"),
      barwidth = unit(0.7, "cm"),
      ticks = TRUE,
      frame.colour = "black"
    )
  )


# ==============================================================================
# 5. COMBINE A + B
# ==============================================================================

p_climate_pca_combined <-
  p_A + p_B +
  
  plot_layout(
    ncol = 2,
    
    # Slightly more room for panel B because of elevation legend
    widths = c(1, 1.08)
  ) +
  
  plot_annotation(
    tag_levels = "A"
  )


# The figure is large; it is written to file and opened below instead of being
# printed in the plots pane


# ==============================================================================
# 6. SAVE HIGH-QUALITY PNG
# ==============================================================================

ggsave(
  filename = this_module_dir(
    "figures",
    "climate_pca_combined.png"
  ),
  plot = p_climate_pca_combined,
  
  width = 16,
  height = 8,
  units = "in",
  
  dpi = 600,
  scale = 1.25,
  
  bg = "white"
)


cat(
  "\nSaved:\n",
  this_module_dir(
    "figures",
    "climate_pca_combined.png"
  ),
  "\n"
)







climate_change_summary <- climate_past_present_wide %>%
  summarise(
    # Mean annual temperature
    MAT_past = mean(MAT_past, na.rm = TRUE),
    MAT_present = mean(MAT_present, na.rm = TRUE),
    delta_MAT = mean(delta_MAT, na.rm = TRUE),
    
    # Summer temperature — selected predictor
    T_JJA_past = mean(temp_JJA_past, na.rm = TRUE),
    T_JJA_present = mean(temp_JJA_present, na.rm = TRUE),
    delta_T_JJA = mean(delta_temp_JJA, na.rm = TRUE),
    
    # Annual precipitation — selected predictor
    MAP_past = mean(MAP_past, na.rm = TRUE),
    MAP_present = mean(MAP_present, na.rm = TRUE),
    delta_MAP = mean(delta_MAP, na.rm = TRUE)
  )

print(climate_change_summary)

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(this_module_dir("figures", "climate_pca_combined.png"))