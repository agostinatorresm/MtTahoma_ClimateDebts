# =============================================================================
# 05e_residual_associations_figure.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the local computer, after
#           26_residual_associations.R)
# Author: Agostina Torres
#
# Residual species associations at the plot level (correlations in the
# plot-level random effect of the final HMSC fits): overstory (top row) and
# understory (bottom row), occurrence (left column) and abundance (right
# column); row titles at the left, column titles once at the top. Lower triangle
# of each correlation matrix; cells are coloured by the posterior mean
# correlation when the probability of direction f >= 0.80 (posterior
# probability >= 0.80 that the correlation is positive, or negative), and left
# blank otherwise. Species names are coloured by the trait most strongly related
# to the residual associations in each layer (name_trait; key below each panel).
#
# Species are ordered by hierarchical clustering of the correlations (each
# panel its own clustering; average linkage, distance = 1 - correlation, as in
# 26_residual_associations.R). The dendrogram is drawn above each matrix and
# the groups obtained by cutting the tree into n_groups are outlined on the
# diagonal (G1, G2, ...). An inset in the empty corner of each panel with groups
# shows the trait profile of the groups (higher or lower than the mean of all
# species: upward or downward triangles). Orders by trait (species_orders.csv) can be added to
# orders_to_plot; a trait missing from a model falls back to the clustering.
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/associations/residual_associations.csv (26)
#          pipeline/02_analysis/02_hmsc/outputs/associations/species_orders.csv (26)
#          pipeline/02_analysis/02_hmsc/outputs/associations/species_traits.csv (26; inset)
#          pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/
#            data_clean/species_climatic_affinity.csv (name colours)
# Outputs: pipeline/02_analysis/02_hmsc/outputs/associations/cluster_groups.csv
#          pipeline/02_analysis/02_hmsc/figures/m1_final/
#            ms1__fig5_residual_associations_m1_final.png (Figure S6.5)
# =============================================================================

library(here)
library(tidyverse)
library(hrbrthemes)
library(patchwork)

tag        <- "m1_final"
fig_prefix <- "ms1__"
f_level    <- 0.80
# species orders to plot (trait orders, e.g. MAT_mean = "Ordered by thermal
# affinity", can be added to this vector)
orders_to_plot <- c(cluster = "Clustered by residual association")
fig_names <- c(cluster = "fig5_residual_associations",
               MAT_mean = "figS_residual_associations_MAT_mean",
               MAP_mean = "figS_residual_associations_MAP_mean",
               seed_mass = "figS_residual_associations_seed_mass")
# number of groups outlined in the clustering order, per model
n_groups <- c("overstory occurrence" = 2, "overstory abundance" = 2,
              "understory occurrence" = 5, "understory abundance" = 1)
# models cut at a fixed height instead (distance = 1 - correlation): used when
# most species have no associations, so cutting into k groups only separates
# single species; a height of 0.7 keeps groups whose species are joined at a
# mean correlation above 0.3
cut_heights <- c("understory abundance" = 0.7)
min_group_size <- 3     # groups with fewer species are not outlined
# trait used to colour the species names in each panel: the trait most strongly
# related to the residual associations (Mantel tests of association vs trait
# difference: moisture affinity in the overstory, SLA in the understory)
name_trait <- c("overstory occurrence" = "MAP_mean", "overstory abundance" = "MAP_mean",
                "understory occurrence" = "SLA", "understory abundance" = "SLA")
# Names use the end colours of the moisture scale of Figure 3 (dry - wet) and a
# brown scale for SLA (low - high), with a grey or mid-tone centre so that all
# names stay readable, and neither overlaps the purple-green scale of the
# correlations
name_palettes <- list(MAP_mean = c("#C98B2B", "#808080", "#0B3C5D"),                        # dry - wet
                      SLA      = c("#D9A441", "#9C6B30", "#4A2E14"),                        # low - high
                      MAT_mean = viridisLite::plasma(256, end = 0.80))
name_keys <- c(MAP_mean = "Names: dry \u2192 wet affinity", SLA = "Names: low \u2192 high SLA",
               MAT_mean = "Names: cold \u2192 warm affinity")
dendro_height <- 0.18   # dendrogram height as a proportion of the matrix side
# Correlations on a purple-green diverging scale (negative - none - positive),
# so that they are not read as cold and warm (blue and orange are used for the
# climatic affinities in the other figures)
col_assoc    <- c("#762A83", "#C2A5CF", "white", "#A6DBA0", "#1B7837")
col_affinity <- viridisLite::plasma(256, end = 0.80, direction = 1)   # as in 05b
mat_limits   <- c(0, 12.5)
strata_titles <- c(overstory = "\U0001F332 Overstory", understory = "\U0001F33F Understory")

if (.Platform$OS.type == "windows") windowsFonts(`Arial Narrow` = windowsFont("Arial Narrow"))
if (!interactive()) png(tempfile(fileext = ".png"))

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
in_dir  <- file.path(module_dir, "outputs", "associations")
fig_dir <- file.path(module_dir, "figures", tag)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

assoc  <- read_csv(file.path(in_dir, "residual_associations.csv"), show_col_types = FALSE)
orders <- read_csv(file.path(in_dir, "species_orders.csv"), show_col_types = FALSE)
affinity <- read_csv(here("pipeline", "01_data_prep", "05_species_distribution",
                          "02_integrate_gbif_chelsa", "data_clean", "species_climatic_affinity.csv"),
                     show_col_types = FALSE) %>%
  transmute(species = str_replace_all(species, " ", "_"), MAT_raw = MAT_mean)   # species_orders.csv has a MAT_mean rank column
sp_label <- function(x) str_replace_all(x, "_", " ")

# -----------------------------------------------------------------------------
# Clustering helpers: correlation matrix from the pair table, hierarchical
# clustering (as in 26_residual_associations.R: average linkage, distance =
# 1 - correlation), dendrogram segments and groups (tree cut into k groups)
# -----------------------------------------------------------------------------
cor_matrix <- function(a) {
  sp <- sort(unique(c(a$species_1, a$species_2)))
  R <- diag(1, length(sp)); dimnames(R) <- list(sp, sp)
  R[cbind(a$species_1, a$species_2)] <- a$mean
  R[cbind(a$species_2, a$species_1)] <- a$mean
  R
}
# Tree cut: at a fixed height for the models in cut_heights, otherwise into k groups
cut_tree <- function(hc, key) {
  if (key %in% names(cut_heights)) cutree(hc, h = cut_heights[[key]]) else cutree(hc, k = n_groups[[key]])
}
# Segments of the dendrogram, leaves at x = 1..n in hc$order, y = merge height
hc_segments <- function(hc) {
  n <- length(hc$order); pos_leaf <- match(seq_len(n), hc$order)
  xm <- numeric(n - 1); segs <- vector("list", n - 1)
  node_xy <- function(j) if (j < 0) c(pos_leaf[-j], 0) else c(xm[j], hc$height[j])
  for (i in seq_len(n - 1)) {
    a <- node_xy(hc$merge[i, 1]); b <- node_xy(hc$merge[i, 2]); h <- hc$height[i]
    xm[i] <- (a[1] + b[1]) / 2
    segs[[i]] <- tibble(x = c(a[1], b[1], a[1]), xend = c(a[1], b[1], b[1]),
                        y = c(a[2], b[2], h), yend = c(h, h, h))
  }
  bind_rows(segs)
}

# -----------------------------------------------------------------------------
# Panel: lower-triangle heat map of one model in a given species order. For the
# clustering order, the dendrogram is drawn above the matrix and the groups of
# the tree cut into n_groups are outlined on the diagonal (G1, G2, ...).
# -----------------------------------------------------------------------------
assoc_panel <- function(st, resp, order_name, show_legend = FALSE) {
  a <- filter(assoc, stratum == st, response == resp)
  o <- filter(orders, stratum == st, response == resp)
  used <- if (order_name %in% names(o) && !all(is.na(o[[order_name]]))) order_name else "cluster"
  if (used == "cluster") {
    R  <- cor_matrix(a)
    hc <- hclust(as.dist(1 - R), method = "average")
    o  <- tibble(species = rownames(R)[hc$order], pos = seq_along(hc$order))
    grp <- cut_tree(hc, paste(st, resp))[o$species]
    o$group <- match(grp, unique(grp))                      # groups numbered along the order
  } else {
    o <- o %>% arrange(.data[[used]]) %>% mutate(pos = row_number(), group = NA_integer_)
  }
  o <- left_join(o, affinity, by = "species")
  n <- nrow(o)
  # both triangles from the pair table, then keep the lower one in this order
  cells <- bind_rows(a, rename(a, species_1 = species_2, species_2 = species_1)) %>%
    mutate(x = o$pos[match(species_1, o$species)], y = o$pos[match(species_2, o$species)]) %>%
    filter(x < y, f >= f_level)
  title <- if (used != order_name) "(clustered: trait not in this model)" else NULL
  # species names coloured by the trait set in name_trait (standardized values);
  # thermal affinity when the trait table is not available
  nt <- name_trait[[paste(st, resp)]]
  if (!is.null(traits) && nt %in% names(traits)) {
    v <- traits %>% filter(stratum == st) %>% {setNames(.[[nt]], .$species)}
    v <- v[o$species]
    lim <- max(abs(v), na.rm = TRUE) * c(-1, 1)
    name_cols <- scales::col_numeric(name_palettes[[nt]], lim, na.color = "grey50")(v)
    key <- name_keys[[nt]]
  } else {
    name_cols <- scales::col_numeric(col_affinity, mat_limits)(pmin(pmax(o$MAT_raw, 0), 12.5))
    key <- name_keys[["MAT_mean"]]
  }
  
  p <- ggplot() +
    geom_tile(data = cells, aes(x = x, y = y, fill = mean), width = 0.95, height = 0.95) +
    annotate("rect", xmin = 0.5, xmax = n + 0.5, ymin = 0.5, ymax = n + 0.5,
             fill = NA, colour = "grey70", linewidth = 0.3)
  y_top <- n + 0.5
  if (used == "cluster") {
    # groups on the diagonal
    boxes <- o %>% group_by(group) %>%
      summarise(lo = min(pos) - 0.5, hi = max(pos) + 0.5, n_sp = n(), .groups = "drop") %>%
      filter(n_sp >= min_group_size, n_sp < n) %>%
      mutate(label = paste0("G", row_number()))
    # dendrogram above the matrix (height = dendro_height x n)
    seg <- hc_segments(hc)
    sc  <- dendro_height * n / max(hc$height)
    seg <- seg %>% mutate(y = n + 0.5 + y * sc, yend = n + 0.5 + yend * sc)
    y_top <- n + 0.5 + dendro_height * n * 1.02
    p <- p +
      geom_rect(data = boxes, aes(xmin = lo, xmax = hi, ymin = lo, ymax = hi),
                fill = NA, colour = "black", linewidth = 0.5) +
      # label inside the box, in its empty lower-right corner
      geom_text(data = boxes, aes(x = hi - 0.3, y = lo + 0.3, label = label),
                hjust = 1, vjust = 0, size = if (st == "overstory") 4.5 else 3.5, fontface = "bold") +
      geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
                   colour = "grey35", linewidth = 0.3)
  }
  p +
    scale_fill_gradientn(colours = col_assoc, limits = c(-1, 1), name = "Residual\ncorrelation",
                         guide = guide_colourbar(barheight = unit(4, "cm"))) +
    scale_x_continuous(breaks = o$pos, labels = sp_label(o$species), expand = c(0, 0)) +
    scale_y_continuous(breaks = o$pos, labels = sp_label(o$species), expand = c(0, 0)) +
    coord_fixed(xlim = c(0.5, n + 0.5), ylim = c(0.5, y_top), clip = "off") +
    labs(x = NULL, y = NULL, title = title) +
    theme_ipsum(base_size = 11) +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic",
                                     size = if (st == "overstory") 10 else 6.5, colour = name_cols),
          axis.text.y = element_text(face = "italic", size = if (st == "overstory") 10 else 6.5, colour = name_cols),
          plot.title = element_text(size = 10, face = "plain", colour = "grey40"),
          plot.caption = element_text(size = 10, colour = "grey30", hjust = 0, face = "plain"),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = if (show_legend) "right" else "none")
}

# Cluster groups per model (for the Results), written next to the inputs
groups <- map_dfr(names(n_groups), function(key) {
  st <- word(key, 1); resp <- word(key, 2)
  R  <- cor_matrix(filter(assoc, stratum == st, response == resp))
  hc <- hclust(as.dist(1 - R), method = "average")
  sp <- rownames(R)[hc$order]
  grp <- cut_tree(hc, key)[sp]
  g <- match(grp, unique(grp))
  size <- table(g)[as.character(g)]
  keep <- size >= min_group_size & size < length(sp)
  lab <- rep(NA_character_, length(sp))
  lab[keep] <- paste0("G", match(g[keep], unique(g[keep])))
  tibble(stratum = st, response = resp, species = sp, position = seq_along(sp), group = lab)
}) %>% left_join(affinity, by = "species")
write_csv(groups, file.path(in_dir, "cluster_groups.csv"))
# -----------------------------------------------------------------------------
# Trait profile of the groups (inset, drawn in the empty lower-right corner of
# each panel with groups): for each group and trait, an upward triangle when
# the group mean is above the mean of all species, downward when below, sized by
# the difference (standardized traits, SD units); a small dot when the
# difference is smaller than 0.25 SD
# -----------------------------------------------------------------------------
trait_file <- file.path(in_dir, "species_traits.csv")
traits <- if (file.exists(trait_file)) read_csv(trait_file, show_col_types = FALSE) else NULL
trait_labels <- c(MAT_mean = "Warm", MAP_mean = "Wet", SLA = "SLA", woodiness = "Woody",
                  seed_mass = "Seed\nmass", wood_density = "Wood\ndensity")
# The trait is given by the column label, so the symbols are drawn in a single
# colour; the triangle gives the direction and its size the difference
inset_col <- "grey20"

profile_inset <- function(st, resp) {
  g <- filter(groups, stratum == st, response == resp, !is.na(group))
  if (is.null(traits) || n_distinct(g$group) < 2) return(NULL)
  tr <- filter(traits, stratum == st) %>% select(-stratum)
  tr <- tr[, c("species", intersect(names(trait_labels), names(tr)))]
  tr <- tr[, colSums(is.na(tr)) == 0, drop = FALSE]
  d <- g %>% select(species, group) %>% left_join(tr, by = "species") %>%
    pivot_longer(-c(species, group), names_to = "trait", values_to = "z") %>%
    group_by(group, trait) %>% summarise(z = mean(z, na.rm = TRUE), .groups = "drop") %>%
    mutate(trait = factor(trait, levels = names(trait_labels)),
           shape = case_when(z >= 0.25 ~ "up", z <= -0.25 ~ "down", TRUE ~ "none"),
           size  = if_else(shape == "none", 1, 2 + 4 * pmin(abs(z), 1.2)))
  ggplot(d, aes(x = trait, y = factor(group, levels = rev(sort(unique(group)))))) +
    geom_point(aes(shape = shape, size = size), fill = inset_col, colour = inset_col, stroke = 0.3) +
    scale_shape_manual(values = c(up = 24, down = 25, none = 16), guide = "none") +
    scale_size_identity() +
    scale_x_discrete(labels = trait_labels, position = "top") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x.top = element_text(size = 8, colour = "black", lineheight = 0.85),
          axis.text.y = element_text(size = 9, face = "bold", colour = "black"),
          plot.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3),
          plot.margin = margin(3, 3, 3, 3))
}

# Colour key of the species names: a small colour bar with the two ends of the
# trait named, drawn in the empty half of each matrix, just above the
# trait-profile inset
name_key_labels <- list(MAP_mean = c("Dry", "Wet"), SLA = c("Low SLA", "High SLA"),
                        MAT_mean = c("Cold", "Warm"))
name_trait_used <- function(st, resp) {
  nt <- name_trait[[paste(st, resp)]]
  if (!is.null(traits) && nt %in% names(traits)) nt else "MAT_mean"
}
name_key_plot <- function(nt) {
  pal <- if (nt == "MAT_mean") col_affinity else name_palettes[[nt]]
  lab <- name_key_labels[[nt]]
  ggplot(tibble(x = seq(0, 1, length.out = 100)), aes(x = x, y = 1, fill = x)) +
    geom_tile(height = 1) +
    scale_fill_gradientn(colours = pal, guide = "none") +
    annotate("text", x = -0.04, y = 1, label = lab[1], hjust = 1, size = 3.6, colour = "black") +
    annotate("text", x = 1.04, y = 1, label = lab[2], hjust = 0, size = 3.6, colour = "black") +
    coord_cartesian(xlim = c(-0.55, 1.55), ylim = c(0.5, 1.5), expand = FALSE, clip = "off") +
    labs(title = "Species names") +
    theme_void() +
    theme(plot.title = element_text(size = 10, hjust = 0.5, colour = "black"))
}

# panel with the colour key of the names and the trait-profile inset in the
# lower-right corner (empty half of the matrix)
panel_with_inset <- function(st, resp, order_name, show_legend = FALSE) {
  # position of the key (panel units: left, right, bottom, top), clear of the
  # group boxes and their labels
  key_pos <- if (st == "overstory") c(0.45, 0.72, 0.26, 0.31) else c(0.70, 0.97, 0.355, 0.405)
  p <- assoc_panel(st, resp, order_name, show_legend) +
    inset_element(name_key_plot(name_trait_used(st, resp)),
                  left = key_pos[1], right = key_pos[2], bottom = key_pos[3], top = key_pos[4],
                  align_to = "panel")
  ins <- if (order_name == "cluster") profile_inset(st, resp) else NULL
  if (is.null(ins)) return(p)
  # overstory: lower inset, below the bottom edge of the largest group box
  top <- if (st == "overstory") 0.235 else 0.33
  p + inset_element(ins, left = 0.58, bottom = 0.02, right = 0.99, top = top, align_to = "panel")
}

# -----------------------------------------------------------------------------
# One figure per order
# -----------------------------------------------------------------------------
open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
for (on in names(orders_to_plot)) {
  # column titles (once, at the top) and row titles (vertical, at the left)
  txt <- function(label, rot = 0, size = 18) wrap_elements(grid::textGrob(label, rot = rot,
                                                                          gp = grid::gpar(fontsize = size, fontface = "bold")))
  fig <- wrap_plots(
    plot_spacer(), txt("Occurrence"), txt("Abundance"),
    txt(strata_titles[["overstory"]], rot = 90),
    panel_with_inset("overstory", "occurrence", on), panel_with_inset("overstory", "abundance", on, show_legend = TRUE),
    txt(strata_titles[["understory"]], rot = 90),
    panel_with_inset("understory", "occurrence", on), panel_with_inset("understory", "abundance", on),
    ncol = 3, widths = c(0.06, 1, 1), heights = c(0.04, 1, 1))
  # Figure S6.5 for the clustering order; other orders as supplementary versions
  fig_file <- file.path(fig_dir, paste0(fig_prefix, fig_names[[on]], "_", tag, ".png"))
  ggsave(fig_file, fig, device = "png", dpi = 300, width = 17, height = 17, bg = "white")
  cat("Figure written:", basename(fig_file), "\n")
}

cat("\nCluster groups (cluster_groups.csv):\n")
print(groups %>% group_by(stratum, response, group) %>%
        summarise(n = n(), mean_MAT = round(mean(MAT_raw, na.rm = TRUE), 1),
                  species = paste(sp_label(species), collapse = ", "), .groups = "drop"),
      n = Inf, width = Inf)

# Open the figures (macOS: open; Windows: default viewer)
walk(file.path(fig_dir, paste0(fig_prefix, fig_names[names(orders_to_plot)], "_", tag, ".png")), open_fig)
