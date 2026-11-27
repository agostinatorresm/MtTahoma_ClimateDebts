# =============================================================================
# 12_phylogeny_traits_figure.R
#
# Pipeline: pipeline/01_data_prep/04_traits/code/ (after 06_build_phylogeny.R
#           and 08_integration.R)
# Author: Agostina Torres
#
# Phylogeny of the species modelled with HMSC (restricted pool: 12 overstory
# and 58 understory species) with the traits used in the models shown next to
# each species, as a visual complement to the phylogenetic signal in Table
# S2.2. The same species and traits as in 11_phylogenetic_signal_traits.R are
# used:
#   overstory : thermal affinity, moisture affinity, seed mass, wood density
#   understory: thermal affinity, moisture affinity, seed mass, SLA, woodiness
#
# Trait values are shown on their natural scale (seed mass and SLA on log10,
# as in the HMSC models), not as z-scores, so that the colours can be read
# against the legends. Climatic affinities use the thermal and moisture scales
# of the other figures (Figure 3, Figure S2.1); each functional trait has its
# own scale, shared by the two layers.
#
# The tree is drawn with ape and ggplot2 (rectangular layout, ladderized), and
# each trait is a separate tile column aligned with the tips, so that every
# trait keeps its own colour scale without additional packages.
#
# Inputs : 04_traits/outputs/MORA_combined_phylo_tree.rda (06_build_phylogeny.R)
#          04_traits/outputs/<layer>_traits.csv (08_integration.R)
#          02_species_pool/outputs/<layer>_species_pool_strict_full_func_clim_cov.csv
#          05_species_distribution/02_integrate_gbif_chelsa/data_clean/species_climatic_affinity.csv
# Output : 04_traits/figures/phylogeny_traits_hmsc_species.png
# =============================================================================

library(tidyverse)
library(here)
library(ape)
library(patchwork)

traits_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)
pool_file  <- function(l) here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                               paste0(l, "_species_pool_strict_full_func_clim_cov.csv"))
dir.create(traits_dir("figures"), recursive = TRUE, showWarnings = FALSE)

# Palettes and limits shared with Figure 3 and Figure S2.1
COL_MAT    <- viridisLite::plasma(256, end = 0.80)
MAT_LIMITS <- c(0, 12.5)
COL_MAP    <- c("#C98B2B", "#E9D8A6", "#8EC3B0", "#1A759F", "#0B3C5D")
MAP_LIMITS <- c(800, 2700)
COL_WOODY  <- c(Woody = "#6B4E31", `Non-woody` = "#A8CF8E")

traits_by_layer <- list(overstory  = c("MAT_mean", "MAP_mean", "seed_mass", "wood_density"),
                        understory = c("MAT_mean", "MAP_mean", "seed_mass", "SLA", "woodiness"))
trait_labels <- c(MAT_mean = "Thermal\naffinity", MAP_mean = "Moisture\naffinity", seed_mass = "Seed\nmass",
                  SLA = "SLA", wood_density = "Wood\ndensity", woodiness = "Woodiness")

# -----------------------------------------------------------------------------
# 1. Tree, traits and affinities
# -----------------------------------------------------------------------------
env <- new.env(); load(traits_dir("outputs", "MORA_combined_phylo_tree.rda"), envir = env)
find_phylo <- function(o) {                       # the .rda may hold the tree directly or inside a list
  if (inherits(o, "phylo")) return(o)
  if (is.list(o)) for (el in o) { r <- find_phylo(el); if (!is.null(r)) return(r) }
  NULL
}
tree_all <- NULL
for (o in mget(ls(env), envir = env)) { tree_all <- find_phylo(o); if (!is.null(tree_all)) break }
if (is.null(tree_all)) stop("No phylo object found in MORA_combined_phylo_tree.rda")
tree_all$tip.label <- str_replace_all(tree_all$tip.label, " ", "_")
tree_all$node.label <- NULL

affinity <- read_csv(here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa",
                          "data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  mutate(species = str_replace_all(species, " ", "_"),
         species = recode(species, Xanthocyparis_nootkatensis = "Callitropsis_nootkatensis")) %>%
  select(species, MAT_mean, MAP_mean)

layer_data <- function(l) {
  traits <- read_csv(traits_dir("outputs", paste0(l, "_traits.csv")), show_col_types = FALSE) %>%
    mutate(species = str_replace_all(species, " ", "_")) %>%
    select(-any_of(c("MAT_mean", "MAP_mean"))) %>% distinct(species, .keep_all = TRUE)
  d <- read_csv(pool_file(l), show_col_types = FALSE) %>%
    transmute(species = str_replace_all(species, " ", "_")) %>%
    left_join(traits, by = "species") %>%
    left_join(affinity, by = "species")
  miss <- setdiff(traits_by_layer[[l]], names(d))
  if (length(miss)) stop(l, ": columns not found: ", paste(miss, collapse = ", "))
  not_in_tree <- setdiff(d$species, tree_all$tip.label)
  if (length(not_in_tree)) warning(l, ": species not in the tree, not shown: ", paste(not_in_tree, collapse = ", "))
  d %>% filter(species %in% tree_all$tip.label) %>%
    mutate(across(any_of(c("seed_mass", "SLA")), log10),
           across(any_of("woodiness"), ~ factor(if_else(tolower(as.character(.x)) %in% c("woody", "1", "true", "yes"),
                                                        "Woody", "Non-woody"), levels = c("Woody", "Non-woody"))))
}
dat <- map(set_names(names(traits_by_layer)), layer_data)
cat("Species shown:", paste(names(dat), map_int(dat, nrow), sep = " = ", collapse = "; "), "\n")

# Limits of the functional-trait scales, shared by the two layers
range_all <- function(v) range(unlist(map(dat, ~ .x[[v]][is.finite(.x[[v]])])), na.rm = TRUE)

# -----------------------------------------------------------------------------
# 2. Tree coordinates (rectangular layout)
# -----------------------------------------------------------------------------
tree_coords <- function(tree) {
  tree <- ladderize(tree, right = FALSE)
  n <- Ntip(tree)
  if (is.null(tree$edge.length)) tree <- compute.brlen(tree)
  node_x <- node.depth.edgelength(tree)                      # distance from the root
  tip_order <- tree$edge[tree$edge[, 2] <= n, 2]              # tips in drawing order
  node_y <- numeric(n + tree$Nnode); node_y[tip_order] <- seq_len(n)
  # internal node at the midpoint of its children, from the smallest clades
  # upwards (a node always has more descendant tips than its children)
  internal <- (n + 1):(n + tree$Nnode)
  for (node in internal[order(node.depth(tree)[internal])])
    node_y[node] <- mean(range(node_y[tree$edge[tree$edge[, 1] == node, 2]]))
  parent <- tree$edge[, 1]; child <- tree$edge[, 2]
  list(h = tibble(x = node_x[parent], xend = node_x[child], y = node_y[child], yend = node_y[child]),
       v = tibble(x = node_x[parent], xend = node_x[parent], y = node_y[child], yend = node_y[parent]),
       tips = tibble(species = tree$tip.label, y = node_y[seq_len(n)]))
}

# -----------------------------------------------------------------------------
# 3. Figure
# -----------------------------------------------------------------------------
theme_col <- theme_void() +
  theme(plot.title = element_text(size = 9, hjust = 0.5, vjust = 0, lineheight = 0.9),
        plot.margin = margin(0, 1, 0, 1))

tile_column <- function(d, var) {
  p <- ggplot(d, aes(x = 1, y = y, fill = .data[[var]])) +
    geom_tile(colour = "white", linewidth = 0.3, height = 0.95) +
    scale_y_continuous(limits = c(0.5, nrow(d) + 0.5), expand = c(0, 0)) +
    labs(title = trait_labels[[var]]) + theme_col
  legend_bar <- guide_colourbar(title.position = "top", barwidth = unit(2.2, "cm"), barheight = unit(0.3, "cm"))
  p + switch(var,
             MAT_mean     = scale_fill_gradientn(colours = COL_MAT, limits = MAT_LIMITS, oob = scales::squish, breaks = c(0, 5, 10),
                                                 name = "Thermal affinity (°C)", guide = legend_bar),
             MAP_mean     = scale_fill_gradientn(colours = COL_MAP, limits = MAP_LIMITS, oob = scales::squish, breaks = c(1000, 2000),
                                                 name = "Moisture affinity (mm)", guide = legend_bar),
             seed_mass    = scale_fill_viridis_c(option = "mako", direction = -1, limits = range_all("seed_mass"), n.breaks = 3,
                                                 na.value = "grey85", name = "Seed mass (log10 mg)", guide = legend_bar),
             SLA          = scale_fill_viridis_c(option = "viridis", limits = range_all("SLA"), n.breaks = 3,
                                                 na.value = "grey85", name = "SLA (log10 mm² mg⁻¹)", guide = legend_bar),
             wood_density = scale_fill_viridis_c(option = "cividis", direction = -1, limits = range_all("wood_density"), n.breaks = 3,
                                                 na.value = "grey85", name = "Wood density (g cm⁻³)", guide = legend_bar),
             woodiness    = scale_fill_manual(values = COL_WOODY, drop = FALSE, name = "Woodiness",
                                              guide = guide_legend(title.position = "top")))
}

layer_panel <- function(l, title, name_size) {
  d0 <- dat[[l]]
  tc <- tree_coords(keep.tip(tree_all, d0$species))
  d <- d0 %>% inner_join(tc$tips, by = "species") %>% mutate(label = str_replace_all(species, "_", " "))
  n <- nrow(d)
  p_tree <- ggplot() +
    geom_segment(data = tc$h, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4, lineend = "square") +
    geom_segment(data = tc$v, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4, lineend = "square") +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.01))) +
    labs(title = title) + theme_void() +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0), plot.margin = margin(0, 0, 0, 2))
  p_names <- ggplot(d, aes(x = 0, y = y, label = label)) +
    geom_text(hjust = 0, size = name_size, fontface = "italic") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(0, 2, 0, 2))
  cols <- map(traits_by_layer[[l]], ~ tile_column(d, .x))
  wrap_plots(c(list(p_tree, p_names), cols), nrow = 1,
             widths = c(2.2, 2.2, rep(0.45, length(cols))))
}

n_over  <- nrow(dat$overstory)
n_under <- nrow(dat$understory)
fig <- (layer_panel("overstory", "A. Overstory", 3.2) / layer_panel("understory", "B. Understory", 2.6)) +
  plot_layout(heights = c(n_over + 2, n_under + 2), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit(0.3, "cm"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7))

fig_file <- traits_dir("figures", "phylogeny_traits_hmsc_species.png")
ggsave(fig_file, fig, width = 8.5, height = 2.5 + 0.16 * (n_over + n_under), dpi = 300, bg = "white", limitsize = FALSE)
cat("Figure written to", fig_file, "\n")

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(fig_file)