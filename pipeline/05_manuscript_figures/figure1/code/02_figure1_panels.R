# =============================================================================
# 02_figure1_panels.R
# Author: Agostina Torres
#
# Step 2 of Figure 1 (step 1: 01_figure1_icons.R, icons of the assembly
# filters in panel C).
#
# Data-like panels of the conceptual Figure 1, drawn in R so that they share
# the style of the result figures (theme_ipsum, Arial Narrow, plasma thermal
# palette, CTI/CMI colours of Figure 2, mountain profile of the species-change
# figures). The panels are exported as transparent PNGs (600 dpi) to be placed
# in the PowerPoint layout of Figure 1, where the text, the icons (step 1), the
# community diagrams and the sheet outlines are drawn.
#
# Panels (all conceptual, no data):
#   A  - species distributions along elevation (abundance curves coloured by
#        thermal affinity; filled = mesic, dashed = dry) next to the mountain
#        profile with the resurveyed plots; dry and mesic species alternate
#        irregularly along elevation, while thermal affinity decreases upslope
#   B  - expected change in the community temperature (dCTI) and moisture
#        (dCMI) indices along elevation; changes are larger at low elevation,
#        where warming was greater (same slope in both layers), and larger in
#        the understory; the layer labels state the demographic pace (slow in
#        the overstory, fast in the understory)
#   D  - expected trait effects on species' change over time: seed mass (-),
#        SLA (+), thermal affinity (+), moisture affinity (-); three ticks per
#        axis (-1, 0, 1)
#   Keys - thermal-affinity colour gradient, moisture-affinity gradient (from
#          open with a dashed outline, dry, to filled with a solid outline,
#          wet, as in panels A, B and D) and forest layers (tree and herb
#          icons), also combined in one row to place under the figure
#
# Outputs: pipeline/05_manuscript_figures/figure1/panels/fig1_*.png
# =============================================================================

library(tidyverse)
library(here)
library(hrbrthemes)
library(patchwork)

out_dir <- here("pipeline", "05_manuscript_figures", "figure1", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 0. Style shared with the result figures
# -----------------------------------------------------------------------------
font <- "Arial Narrow"
if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))

col_affinity <- viridisLite::plasma(256, end = 0.80, direction = 1)   # thermal affinity (as in 05b/05d)
mat_limits   <- c(0, 12.5)
col_cti <- "#F26B5B"                                                 # as in Figure 2
col_cmi <- "#1A759F"
emoji_tree <- "\U0001F332"
emoji_leaf <- "\U0001F33F"

theme_fig1 <- function(base = 14) {
  theme_ipsum(base_family = font, base_size = base, axis_title_size = base, plot_margin = margin(4, 4, 4, 4)) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.3),
          axis.text = element_text(colour = "black", size = base - 1),
          axis.title.x = element_text(colour = "black", hjust = 0.5),
          axis.title.y = element_text(colour = "black", hjust = 0.5),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA))
}

# Transparent PNG (600 dpi; ragg renders colour emoji when available) and SVG
save_panel <- function(p, name, w, h) {
  f <- file.path(out_dir, paste0("fig1_", name, ".png"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(f, p, width = w, height = h, dpi = 600, bg = "transparent", device = ragg::agg_png)
  } else {
    ggsave(f, p, width = w, height = h, dpi = 600, bg = "transparent")
  }
  # vector copy (scales without loss in PowerPoint: Insert > Pictures)
  if (requireNamespace("svglite", quietly = TRUE)) ggsave(sub("\\.png$", ".svg", f), p, width = w, height = h, bg = "transparent")
  f
}

# -----------------------------------------------------------------------------
# 1. Panel A: species distributions and mountain profile
# -----------------------------------------------------------------------------
# Mountain profile: same ridge as in the species-change figures (05b, 05d)
displace <- function(x, y, iters, sd, decay = 0.55) {
  for (i in seq_len(iters)) {
    n  <- length(x)
    mx <- (x[-n] + x[-1]) / 2
    my <- (y[-n] + y[-1]) / 2 + rnorm(n - 1, 0, sd)
    x  <- c(rbind(x[-n], mx), x[n])
    y  <- c(rbind(y[-n], my), y[n])
    sd <- sd * decay
  }
  tibble(x = x, y = y)
}
set.seed(7)
ridge <- displace(c(0, 0.16, 0.30, 0.42, 0.52, 0.60, 0.70, 0.80, 0.90, 1.00),
                  c(0.58, 0.66, 0.80, 0.90, 1.00, 0.88, 0.70, 0.48, 0.24, 0.00), iters = 6, sd = 0.045)
strokes <- list(list(x = c(0.52, 0.50, 0.48), y = c(0.99, 0.93, 0.88)),
                list(x = c(0.56, 0.58, 0.57), y = c(0.94, 0.89, 0.84)))
rocks <- imap_dfr(strokes, ~ displace(.x$x, .x$y, iters = 3, sd = 0.015) %>% mutate(g = .y))

# The profile from x = 0.30 onwards (as in the species-change figures) is placed
# to the right of the abundance curves
m_x0 <- 0.42; m_w <- 0.58
to_panel <- function(d) mutate(d, x = m_x0 + (x - 0.30) / 0.70 * m_w)
ridge_p <- to_panel(filter(ridge, x >= 0.30))
rocks_p <- to_panel(rocks)

# Resurveyed plots below the ridge
plots <- tibble(x = m_x0 + c(0.05, 0.12, 0.18, 0.24, 0.30, 0.36, 0.42, 0.48, 0.53, 0.27),
                drop = c(0.10, 0.18, 0.08, 0.22, 0.12, 0.15, 0.10, 0.12, 0.06, 0.40)) %>%
  mutate(y = approx(ridge_p$x, ridge_p$y, xout = x)$y - drop, y = pmax(y, 0.05))

# Abundance curves along elevation (x = abundance, y = elevation)
# Species with irregular positions, widths and heights; dry-affiliated
# (dashed) and mesic (filled) species follow no fixed order along the
# gradient, while thermal affinity (MAT) decreases with elevation.
sp <- tibble(sp  = 1:8,
             mu  = c(0.13, 0.23, 0.34, 0.45, 0.55, 0.65, 0.76, 0.87),
             sd  = c(0.055, 0.060, 0.075, 0.065, 0.070, 0.060, 0.065, 0.050),
             amp = c(0.30, 0.21, 0.27, 0.19, 0.29, 0.23, 0.27, 0.20),
             dry = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE)) %>%
  mutate(MAT = 12 - 11 * mu)                                   # warm species low, cold species high
x_base <- 0.02
curves <- sp %>%
  crossing(e = seq(0, 1, by = 0.004)) %>%
  mutate(g = exp(-0.5 * ((e - mu) / sd)^2)) %>%
  filter(g > 0.03) %>%
  mutate(x = x_base + amp * g)
mesic_poly <- curves %>% filter(!dry) %>% group_by(sp) %>%
  group_modify(~ bind_rows(tibble(x = x_base, e = min(.x$e), MAT = .x$MAT[1]),
                           select(.x, x, e, MAT),
                           tibble(x = x_base, e = max(.x$e), MAT = .x$MAT[1]))) %>%
  ungroup()

pA <- ggplot() +
  geom_polygon(data = mesic_poly, aes(x, e, group = sp, fill = MAT), alpha = 0.45, colour = NA) +
  geom_path(data = filter(curves, !dry), aes(x, e, group = sp, colour = MAT), linewidth = 0.6) +
  geom_path(data = filter(curves, dry), aes(x, e, group = sp, colour = MAT), linewidth = 0.7, linetype = "22") +
  scale_colour_gradientn(colours = col_affinity, limits = mat_limits, oob = scales::squish, guide = "none") +
  scale_fill_gradientn(colours = col_affinity, limits = mat_limits, oob = scales::squish, guide = "none") +
  geom_path(data = ridge_p, aes(x, y), colour = "grey50", alpha = 0.6, linewidth = 0.45,
            lineend = "round", linejoin = "round") +
  geom_path(data = rocks_p, aes(x, y, group = g), colour = "grey50", alpha = 0.6, linewidth = 0.3) +
  geom_point(data = plots, aes(x, y), size = 1.3, colour = "grey15") +
  annotate("text", x = m_x0 + 0.03, y = 0.06, label = "324 plots\n~1975 vs ~2015", hjust = 0, vjust = 0,
           family = font, size = 5.8, colour = "grey30", lineheight = 0.9) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = 1.04, linewidth = 0.4,
           arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("text", x = -0.05, y = 0.5, label = "Elevation", angle = 90, family = font, size = 5.3, colour = "grey30") +
  annotate("text", x = x_base, y = -0.05, label = "Occurrence / abundance", hjust = 0, vjust = 1,
           family = font, size = 5.3, colour = "grey30") +
  coord_cartesian(xlim = c(-0.08, 1), ylim = c(-0.17, 1.06), expand = FALSE, clip = "off") +
  theme_void(base_family = font) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA), plot.margin = margin(4, 4, 6, 8))

fA <- save_panel(pA, "A_distributions_mountain", 4.4, 2.9)

# -----------------------------------------------------------------------------
# 2. Panel B: expected change in CTI and CMI along elevation
# -----------------------------------------------------------------------------
# Changes are larger at low elevation and in the understory; overstory drawn
# opaque, understory transparent (as in the arrows of the conceptual figure)
b_data <- function(f_over, f_under) {
  tibble(elev = seq(650, 1800, length.out = 80)) %>%
    mutate(u = (elev - 650) / 1150, over = f_over(u), under = f_under(u)) %>%
    pivot_longer(c(over, under), names_to = "stratum", values_to = "delta")
}
plot_b <- function(d, col, ylab, ylim) {
  ends <- d %>% group_by(stratum) %>% filter(elev == max(elev)) %>% ungroup() %>%
    mutate(lab = if_else(stratum == "over", paste(emoji_tree, "slow"), paste(emoji_leaf, "fast")))
  ggplot(d, aes(elev, delta)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_line(data = filter(d, stratum == "under"), colour = col, alpha = 0.45, linewidth = 1, lineend = "round") +
    geom_line(data = filter(d, stratum == "over"), colour = col, linewidth = 1, lineend = "round") +
    geom_text(data = ends, aes(x = elev + 35, label = lab), hjust = 0, size = 4.9, family = font, colour = "grey25") +
    scale_x_continuous(breaks = c(750, 1250, 1750)) +
    scale_y_continuous(breaks = 0) +
    coord_cartesian(xlim = c(650, 1800), ylim = ylim, clip = "off") +
    labs(x = "Elevation (m)", y = ylab) +
    theme_fig1() +
    theme(plot.margin = margin(4, 60, 4, 4))
}
# same slope with elevation in both layers; larger change in the understory
d_cti <- b_data(function(u) 0.08 + 0.18 * (1 - u), function(u) 0.20 + 0.18 * (1 - u))
d_cmi <- b_data(function(u) -(8 + 18 * (1 - u)), function(u) -(20 + 18 * (1 - u)))
# Inset: replacement of species over time (cold- by warm-affiliated for CTI,
# wet- by dry-affiliated for CMI); dry = open dashed circle
circle_df <- function(x0, r = 0.28) { t <- seq(0, 2 * pi, length.out = 90); tibble(x = x0 + r * cos(t), y = r * sin(t)) }
repl_inset <- function(colA, fillA, ltyA, colB, fillB, ltyB, labA, labB, arrow_col) {
  ggplot() +
    geom_polygon(data = circle_df(0), aes(x, y), fill = fillA, colour = colA, linetype = ltyA, linewidth = 0.6) +
    geom_polygon(data = circle_df(2), aes(x, y), fill = fillB, colour = colB, linetype = ltyB, linewidth = 0.6) +
    annotate("segment", x = 0.42, xend = 1.55, y = 0, yend = 0, colour = arrow_col, linewidth = 0.6,
             arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
    annotate("text", x = 1, y = 0.24, label = "over time", family = font, size = 3.6, colour = "grey35", vjust = 0) +
    annotate("text", x = 0, y = -0.45, label = labA, family = font, size = 4.4, colour = colA, vjust = 1) +
    annotate("text", x = 2, y = -0.45, label = labB, family = font, size = 4.4, colour = colB, vjust = 1) +
    coord_fixed(xlim = c(-0.7, 2.7), ylim = c(-0.85, 0.55), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA))
}
col_cold <- col_affinity[round(1.5 / 12.5 * 255) + 1]
col_warm <- col_affinity[256]
inset_cti <- repl_inset(col_cold, col_cold, "solid", col_warm, col_warm, "solid", "cold", "warm", col_cti)
inset_cmi <- repl_inset("black", "black", "solid", "black", NA, "22", "wet", "dry", col_cmi)

# CTI: lines are high at low elevation, so the inset goes to the upper right;
# CMI: lines are low at low elevation, so the inset goes to the lower right
pB_cti <- plot_b(d_cti, col_cti, "\u0394CTI (\u00b0C/decade)", c(-0.04, 0.50)) +
  inset_element(inset_cti, left = 0.42, bottom = 0.6, right = 1.0, top = 1.02, align_to = "panel") &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))
pB_cmi <- plot_b(d_cmi, col_cmi, "\u0394CMI (mm/decade)", c(-50, 4)) +
  inset_element(inset_cmi, left = 0.42, bottom = 0.0, right = 1.0, top = 0.42, align_to = "panel") &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

fB1 <- save_panel(pB_cti, "B_delta_CTI", 3.6, 2.3)
fB2 <- save_panel(pB_cmi, "B_delta_CMI", 3.6, 2.3)

# -----------------------------------------------------------------------------
# 3. Panel D: expected trait effects on species' change over time
# -----------------------------------------------------------------------------
# species scattered around the expected relationship (fixed seed)
set.seed(11)
pts <- tibble(u = seq(0.07, 0.93, length.out = 14) + runif(14, -0.025, 0.025), jit = rnorm(14, 0, 0.11))
# Fixed aspect (one y unit drawn as mini_ratio x units) so that the moisture
# circles stay round; dot() draws one circle of radius r (in x units)
mini_ratio <- 0.65
dot <- function(x0, y0, r = 0.026) {
  t <- seq(0, 2 * pi, length.out = 40)
  tibble(x = x0 + r * cos(t), y = y0 + r / mini_ratio * sin(t))
}
# moisture steps as in the key (moist_steps, section 4), assigned by affinity
moist_points <- function(d) {
  step <- cut(d$u, breaks = seq(0, 1, length.out = 6), labels = FALSE, include.lowest = TRUE)
  tibble(x0 = d$u, y0 = d$v,
         fill_alpha = c(0, 0.25, 0.5, 0.75, 1)[step],
         outline = c("12", "22", "42", "82", "solid")[step])
}
mini <- function(slope, xlab, type = "grey") {
  d <- mutate(pts, v = 0.5 + slope * (u - 0.5) + jit,
              MAT = 1 + 11 * u, dry = u < 0.5, v = pmin(pmax(v, 0.04), 0.96))
  p <- ggplot(d, aes(u, v)) +
    geom_hline(yintercept = 0.5, linetype = "22", colour = "grey75", linewidth = 0.3) +
    annotate("segment", x = 0.08, xend = 0.92, y = 0.5 + slope * (0.08 - 0.5), yend = 0.5 + slope * (0.92 - 0.5),
             linewidth = 0.7, lineend = "round")
  p <- p + switch(type,
                  grey    = geom_point(colour = "grey60", size = 2.1),
                  thermal = geom_point(aes(colour = MAT), size = 2.1),
                  # moisture: the five steps of the key, from open with a dashed
                  # outline (dry-affiliated) to filled with a solid outline
                  # (wet-affiliated), drawn as small circles
                  moist   = pmap(moist_points(d), function(x0, y0, fill_alpha, outline)
                    geom_polygon(data = dot(x0, y0), aes(x, y), fill = scales::alpha("black", fill_alpha),
                                 colour = "black", linetype = outline, linewidth = 0.45)))
  if (type == "thermal") p <- p + scale_colour_gradientn(colours = col_affinity, limits = mat_limits, guide = "none")
  p +
    # three ticks per axis: standardized trait value and change over time
    # (-1, 0, 1); the axis text and ticks are set for x and y separately because
    # theme_ipsum() defines them per axis
    scale_x_continuous(breaks = c(0.1, 0.5, 0.9), labels = c("\u22121", "0", "1")) +
    scale_y_continuous(breaks = c(0.1, 0.5, 0.9), labels = c("\u22121", "0", "1")) +
    coord_fixed(ratio = mini_ratio, xlim = c(0, 1.02), ylim = c(0, 1), expand = FALSE, clip = "off") +
    labs(x = xlab, y = NULL) +
    theme_fig1(13) +
    theme(axis.text.x = element_text(colour = "grey30", size = 10, margin = margin(t = 1)),
          axis.text.y = element_text(colour = "grey30", size = 10, margin = margin(r = 1)),
          axis.ticks.x = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.y = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.length = unit(2, "pt"),
          axis.line.x = element_line(colour = "black", linewidth = 0.4,
                                     arrow = arrow(length = unit(1.3, "mm"), type = "closed")),
          axis.title.x = element_text(colour = "grey30", size = 13, hjust = 0.5, margin = margin(t = 2)),
          plot.margin = margin(3, 6, 3, 3))
}
# slope: sign of the expected effect on change over time
p_seed  <- mini(-0.7, "seed mass")
p_sla   <- mini( 0.7, "SLA")
p_therm <- mini( 0.7, "thermal affinity", "thermal")
p_moist <- mini(-0.7, "moisture affinity", "moist")
y_lab <- wrap_elements(full = grid::textGrob("Change over time", rot = 90,
                                             gp = grid::gpar(fontfamily = font, fontsize = 13.5, col = "grey30")))
pD <- y_lab + ((p_seed | p_sla) / (p_therm | p_moist)) + plot_layout(widths = c(0.07, 1)) &
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

fD <- save_panel(pD, "D_trait_expectations", 4.2, 2.9)

# -----------------------------------------------------------------------------
# 4. Keys
# -----------------------------------------------------------------------------
# Each key is saved on its own and the three together in one row, to place
# under the figure
key_theme <- theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA), plot.margin = margin(2, 6, 2, 6))
key_title <- function(label) annotate("text", x = 0, y = 0.9, label = label, hjust = 0, vjust = 0,
                                      family = font, size = 4.8, fontface = "bold")

# Thermal affinity: colour gradient
key_thermal <- ggplot(tibble(x = seq(0, 1, length.out = 256)), aes(x, 0, fill = x)) +
  geom_tile(height = 1) +
  scale_fill_gradientn(colours = col_affinity, guide = "none") +
  annotate("text", x = 0, y = -0.9, label = "cold", hjust = 0, vjust = 1, family = font, size = 4.4, colour = col_affinity[1]) +
  annotate("text", x = 1, y = -0.9, label = "warm", hjust = 1, vjust = 1, family = font, size = 4.4, colour = col_affinity[256]) +
  key_title("Thermal affinity") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-2.1, 2.3), expand = FALSE, clip = "off") +
  key_theme

# Moisture affinity: a gradient of circles from open with a dashed outline
# (dry) to filled with a solid outline (wet), linking the dashed and filled
# curves of panel A, the inset of panel B and the graded points of panel D
circle_key <- function(x0, r = 0.3) { t <- seq(0, 2 * pi, length.out = 90); tibble(x = x0 + r * cos(t), y = r * sin(t)) }
moist_steps <- tibble(x0 = c(0.35, 1.45, 2.55, 3.65, 4.75), fill_alpha = c(0, 0.25, 0.5, 0.75, 1),
                      outline = c("12", "22", "42", "82", "solid"))
key_moist <- ggplot() +
  pmap(moist_steps, function(x0, fill_alpha, outline)
    geom_polygon(data = circle_key(x0), aes(x, y), fill = scales::alpha("black", fill_alpha),
                 colour = "black", linetype = outline, linewidth = 0.6)) +
  annotate("text", x = 0, y = -0.6, label = "dry", hjust = 0, vjust = 1, family = font, size = 4.4) +
  annotate("text", x = 5.1, y = -0.6, label = "wet", hjust = 1, vjust = 1, family = font, size = 4.4) +
  annotate("text", x = 0, y = 0.7, label = "Moisture affinity", hjust = 0, vjust = 0,
           family = font, size = 4.8, fontface = "bold") +
  coord_fixed(xlim = c(0, 5.1), ylim = c(-1.3, 1.5), expand = FALSE, clip = "off") +
  key_theme

# Forest layers: the tree and herb icons used in panels B and D and in the
# result figures
key_layers <- ggplot() +
  annotate("text", x = 0, y = 0.75, label = "Forest layer", hjust = 0, vjust = 0,
           family = font, size = 4.8, fontface = "bold") +
  annotate("text", x = 0, y = 0, label = paste(emoji_tree, "overstory (trees)"), hjust = 0, family = font, size = 4.4) +
  annotate("text", x = 2.6, y = 0, label = paste(emoji_leaf, "understory"), hjust = 0, family = font, size = 4.4) +
  coord_cartesian(xlim = c(0, 4.6), ylim = c(-0.8, 1.6), expand = FALSE, clip = "off") +
  key_theme

fK1 <- save_panel(key_thermal, "key_thermal_affinity", 1.9, 0.8)
fK2 <- save_panel(key_moist, "key_moisture_affinity", 1.9, 0.8)
fK3 <- save_panel(key_layers, "key_forest_layers", 2.9, 0.8)
fK  <- save_panel(key_thermal + key_moist + key_layers + plot_layout(widths = c(1.9, 1.9, 2.9)) &
                    theme(plot.background = element_rect(fill = "transparent", colour = NA)),
                  "legend_bottom", 6.9, 0.8)

cat("Panels written to", out_dir, ":\n", paste(basename(c(fA, fB1, fB2, fD, fK1, fK2, fK3, fK)), collapse = "\n "), "\n")

# Open the panels (macOS: Preview; Windows: default viewer)
open_figs <- function(f) if (.Platform$OS.type == "windows") invisible(lapply(f, shell.exec)) else system2("open", shQuote(f))
open_figs(c(fA, fB1, fB2, fD, fK))
