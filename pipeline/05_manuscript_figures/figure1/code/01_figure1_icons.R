# =============================================================================
# 01_figure1_icons.R
# Author: Agostina Torres
#
# Step 1 of Figure 1 (step 2: 02_figure1_panels.R).
#
# Icons for the assembly filters of Figure 1, panel C, drawn in R so that they
# share the colours of the other panels. They are exported as transparent PNGs
# (600 dpi) and SVGs to be placed next to the filters in the PowerPoint layout:
#   - abiotic filter, temperature: thermometer in the CTI colour of Figure 2;
#   - abiotic filter, drought: empty water drop with a dashed outline in the
#     CMI colour, following the convention of the figure (dry = open, dashed);
#   - biotic filter, interactions: a canopy tree shading an understory plant
#     (competition for light);
#   - dispersal filter: three contrasting seeds, from large to small and from
#     gravity to wind dispersal (acorn, winged conifer seed, dandelion achene
#     with pappus), drawn together and also as separate files.
#
# Outputs: pipeline/05_manuscript_figures/figure1/panels/fig1_icon_*.png (and .svg)
# =============================================================================

library(tidyverse)
library(here)

out_dir <- here("pipeline", "05_manuscript_figures", "figure1", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

col_cti   <- "#F26B5B"   # as in Figure 2 and panel B
col_cmi   <- "#1A759F"
col_line  <- "grey25"
col_tree  <- "#5E7F55"
col_plant <- "#8DB27A"

save_icon <- function(p, name, w, h) {
  f <- file.path(out_dir, paste0("fig1_icon_", name, ".png"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(f, p, width = w, height = h, dpi = 600, bg = "transparent", device = ragg::agg_png)
  } else {
    ggsave(f, p, width = w, height = h, dpi = 600, bg = "transparent")
  }
  if (requireNamespace("svglite", quietly = TRUE)) ggsave(sub("\\.png$", ".svg", f), p, width = w, height = h, bg = "transparent")
  f
}
icon_theme <- theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

# -----------------------------------------------------------------------------
# 1. Thermometer (temperature)
# -----------------------------------------------------------------------------
circle <- function(x0, y0, r, n = 100) { t <- seq(0, 2 * pi, length.out = n); tibble(x = x0 + r * cos(t), y = y0 + r * sin(t)) }
# tube with a rounded top: rectangle plus half circle
tube <- function(x0, w, y0, y1, n = 50) {
  t <- seq(0, pi, length.out = n)
  bind_rows(tibble(x = x0 + w, y = y0), tibble(x = x0 + w, y = y1),
            tibble(x = x0 + w * cos(t), y = y1 + w * sin(t)),
            tibble(x = x0 - w, y = y1), tibble(x = x0 - w, y = y0))
}

p_thermo <- ggplot() +
  geom_polygon(data = tube(0, 0.16, 0, 1.55), aes(x, y), fill = "white", colour = col_line, linewidth = 0.9) +
  geom_polygon(data = circle(0, 0, 0.32), aes(x, y), fill = col_cti, colour = col_line, linewidth = 0.9) +
  geom_polygon(data = tube(0, 0.08, 0, 1.05), aes(x, y), fill = col_cti, colour = NA) +
  geom_polygon(data = circle(0, 0, 0.30), aes(x, y), fill = col_cti, colour = NA) +
  annotate("segment", x = 0.16, xend = 0.30, y = c(0.55, 0.85, 1.15, 1.45), yend = c(0.55, 0.85, 1.15, 1.45),
           colour = col_line, linewidth = 0.8, lineend = "round") +
  coord_fixed(xlim = c(-0.45, 0.45), ylim = c(-0.4, 1.8), expand = FALSE) +
  icon_theme
f1 <- save_icon(p_thermo, "thermometer", 0.45, 1.1)

# -----------------------------------------------------------------------------
# 2. Empty drop with dashed outline (drought)
# -----------------------------------------------------------------------------
t <- seq(0, 2 * pi, length.out = 200)
drop <- tibble(x = 0.9 * sin(t) * sin(t / 2)^1.6, y = 1.1 * cos(t))   # pointed tip at the top
p_drop <- ggplot(drop, aes(x, y)) +
  geom_polygon(fill = NA, colour = col_cmi, linewidth = 1.1, linetype = "22") +
  coord_fixed(xlim = c(-0.85, 0.85), ylim = c(-1.25, 1.25), expand = FALSE) +
  icon_theme
f2 <- save_icon(p_drop, "drought_drop", 0.6, 1.0)

# -----------------------------------------------------------------------------
# 3. Canopy tree shading an understory plant (biotic interactions)
# -----------------------------------------------------------------------------
tree <- bind_rows(
  tibble(x = c(-0.55, 0.55, 0), y = c(0.55, 0.55, 1.35), part = 1),
  tibble(x = c(-0.70, 0.70, 0), y = c(0.15, 0.15, 0.95), part = 2),
  tibble(x = c(-0.85, 0.85, 0), y = c(-0.25, -0.25, 0.55), part = 3))
trunk <- tibble(x = c(-0.09, 0.09, 0.09, -0.09), y = c(-0.6, -0.6, -0.25, -0.25))
shade <- tibble(x = c(0.35, 0.85, 1.75, 1.05), y = c(-0.15, -0.25, -0.6, -0.6))   # shadow cast on the plant
leaf  <- function(x0, y0, ang, len = 0.32, wid = 0.12) {
  s <- seq(0, 1, length.out = 40)
  base <- bind_rows(tibble(u = s, v = wid * sin(pi * s)), tibble(u = rev(s), v = -wid * sin(pi * rev(s))))
  tibble(x = x0 + len * (base$u * cos(ang) - base$v * sin(ang)),
         y = y0 + len * (base$u * sin(ang) + base$v * cos(ang)))
}
p_biotic <- ggplot() +
  geom_polygon(data = shade, aes(x, y), fill = "grey60", alpha = 0.5, colour = NA) +
  geom_polygon(data = trunk, aes(x, y), fill = "#7A5A3A", colour = NA) +
  geom_polygon(data = tree, aes(x, y, group = part), fill = col_tree, colour = NA) +
  annotate("segment", x = 1.35, xend = 1.35, y = -0.6, yend = -0.25, colour = col_plant, linewidth = 1, lineend = "round") +
  geom_polygon(data = leaf(1.35, -0.32, 2.4), aes(x, y), fill = col_plant) +
  geom_polygon(data = leaf(1.35, -0.38, 0.6), aes(x, y), fill = col_plant) +
  annotate("segment", x = -1, xend = 2, y = -0.6, yend = -0.6, colour = col_line, linewidth = 0.6) +
  coord_fixed(xlim = c(-1.05, 2.05), ylim = c(-0.7, 1.45), expand = FALSE) +
  icon_theme
f3 <- save_icon(p_biotic, "biotic_shading", 1.2, 0.85)

# -----------------------------------------------------------------------------
# 4. Seeds (dispersal): acorn, winged conifer seed, dandelion
# -----------------------------------------------------------------------------
ellipse <- function(x0, y0, a, b, ang = 0, n = 100) {
  t <- seq(0, 2 * pi, length.out = n)
  tibble(x = x0 + a * cos(t) * cos(ang) - b * sin(t) * sin(ang),
         y = y0 + a * cos(t) * sin(ang) + b * sin(t) * cos(ang))
}
col_nut  <- "#B07A3F"
col_cap  <- "#6E4B2A"
col_wing <- "#D9B97F"

# Acorn: nut with a pointed base, scaly cap and stalk
nut <- tibble(x = 0.40 * cos(seq(0, pi, length.out = 60)), y = 0.02) %>%
  bind_rows(tibble(u = seq(pi, 2 * pi, length.out = 80),
                   x = 0.40 * cos(u), y = 0.02 + 0.70 * sin(u) * (1 - 0.25 * (1 + sin(u)))) %>% select(x, y))
cap <- tibble(u = seq(0, pi, length.out = 60), x = 0.48 * cos(u), y = 0.02 + 0.28 * sin(u)) %>% select(x, y) %>%
  bind_rows(tibble(x = c(-0.48, 0.48), y = c(0.02, 0.02)))
scales_df <- crossing(i = 1:5, j = 1:2) %>%
  mutate(x = -0.36 + (i - 1) * 0.18 + (j - 1) * 0.09, y = 0.08 + (j - 1) * 0.10) %>% filter(abs(x) < 0.40)
acorn <- function(dx = 0) {
  list(geom_polygon(data = mutate(nut, x = x + dx), aes(x, y), fill = col_nut, colour = col_line, linewidth = 0.5),
       geom_polygon(data = mutate(cap, x = x + dx), aes(x, y), fill = col_cap, colour = col_line, linewidth = 0.5),
       geom_point(data = mutate(scales_df, x = x + dx), aes(x, y), shape = 4, size = 0.9, colour = "#4A321C", stroke = 0.5),
       annotate("segment", x = dx, xend = dx + 0.05, y = 0.30, yend = 0.42, colour = col_cap, linewidth = 1.2, lineend = "round"))
}

# Winged conifer seed: small dark seed at the base of an asymmetric wing (one
# nearly straight edge, one curved edge), tilted as when it spins in the air
rot <- function(d, a = -25 * pi / 180) mutate(d, x_r = x * cos(a) - y * sin(a), y = x * sin(a) + y * cos(a), x = x_r) %>% select(-x_r)
wing <- bind_rows(tibble(t = seq(0, 1, length.out = 60)) %>% mutate(x = 0.24 * sin(pi * t^0.65), y = 0.80 * t),
                  tibble(t = seq(1, 0, length.out = 30)) %>% mutate(x = -0.05 * sin(pi * t), y = 0.80 * t)) %>%
  select(x, y) %>% rot()
veins <- tibble(x = 0, y = 0.05, xend = c(0.03, 0.09, 0.15), yend = c(0.70, 0.62, 0.48)) %>%
  rot() %>% rename(x0 = x, y0 = y) %>%
  bind_cols(tibble(x = .$xend, y = .$yend) %>% rot() %>% rename(x1 = x, y1 = y)) %>%
  transmute(x = x0, y = y0, xend = x1, yend = y1)
pine <- function(dx = 0, dy = 0) {
  list(geom_polygon(data = mutate(wing, x = x + dx, y = y + dy), aes(x, y), fill = col_wing, alpha = 0.85, colour = col_line, linewidth = 0.4),
       geom_segment(data = mutate(veins, x = x + dx, xend = xend + dx, y = y + dy, yend = yend + dy),
                    aes(x, y, xend = xend, yend = yend), colour = "#A88A55", linewidth = 0.3),
       geom_polygon(data = ellipse(dx, dy - 0.04, 0.075, 0.12, ang = -25 * pi / 180), aes(x, y), fill = col_cap, colour = col_line, linewidth = 0.4))
}

# Dandelion: small achene, long beak and a pappus of radiating hairs
dandelion <- function(dx = 0, dy = 0) {
  ang <- seq(15, 165, length.out = 13) * pi / 180
  top <- c(dx, dy + 0.62)
  hairs <- tibble(x = top[1], y = top[2], xend = top[1] + 0.36 * cos(ang), yend = top[2] + 0.30 * sin(ang))
  list(geom_polygon(data = ellipse(dx, dy - 0.12, 0.035, 0.11), aes(x, y), fill = "grey35", colour = NA),
       annotate("segment", x = dx, xend = top[1], y = dy - 0.02, yend = top[2], colour = "grey45", linewidth = 0.4),
       geom_segment(data = hairs, aes(x, y, xend = xend, yend = yend), colour = "grey55", linewidth = 0.3),
       geom_point(data = hairs, aes(xend, yend), colour = "grey55", size = 0.25))
}

p_seeds <- ggplot() + acorn(0) + pine(0.95, -0.25) + dandelion(1.95, -0.35) +
  coord_fixed(xlim = c(-0.55, 2.4), ylim = c(-0.75, 0.6), expand = FALSE) + icon_theme
f4 <- save_icon(p_seeds, "seeds_dispersal", 1.45, 0.65)
f5 <- save_icon(ggplot() + acorn(0) + coord_fixed(xlim = c(-0.55, 0.55), ylim = c(-0.75, 0.5), expand = FALSE) + icon_theme,
                "seed_acorn", 0.5, 0.55)
f6 <- save_icon(ggplot() + pine(0, 0) + coord_fixed(xlim = c(-0.15, 0.45), ylim = c(-0.22, 0.78), expand = FALSE) + icon_theme,
                "seed_conifer", 0.33, 0.55)
f7 <- save_icon(ggplot() + dandelion(0, 0) + coord_fixed(xlim = c(-0.42, 0.42), ylim = c(-0.3, 0.98), expand = FALSE) + icon_theme,
                "seed_dandelion", 0.35, 0.55)

cat("Icons written to", out_dir, ":\n", paste(basename(c(f1, f2, f3, f4, f5, f6, f7)), collapse = "\n "), "\n")

# Open the icons (macOS: Preview; Windows: default viewer)
open_figs <- function(f) if (.Platform$OS.type == "windows") invisible(lapply(f, shell.exec)) else system2("open", shQuote(f))
open_figs(c(f1, f2, f3, f4, f5, f6, f7))