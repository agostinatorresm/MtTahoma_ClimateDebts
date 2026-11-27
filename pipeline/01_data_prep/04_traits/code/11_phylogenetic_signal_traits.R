# =============================================================================
# 11_phylogenetic_signal_traits.R
#
# Pipeline: pipeline/01_data_prep/04_traits/code/ (after 06_build_phylogeny.R
# Author: Agostina Torres
#           and 10_plot_traits.R)
#
# Phylogenetic signal of the species traits used in the HMSC models, for the
# species modelled in each layer (frequent species pool with full trait and
# climatic-affinity coverage).
#   - Continuous traits: Blomberg's K and Pagel's lambda (phytools::phylosig),
#     tested by randomisation (K) and likelihood ratio (lambda). SLA and seed
#     mass are log10-transformed, as in the HMSC models.
#   - Woodiness (binary): Fritz & Purvis D (caper::phylo.d); D = 0 is the
#     expectation under Brownian motion and D = 1 a random distribution.
# Estimates with significant signal (p < 0.05) are shown in bold in Table S2.2.
#
# Inputs : 04_traits/outputs/MORA_combined_phylo_tree.rda (06_build_phylogeny.R)
#          pipeline/01_data_prep/02_species_pool/outputs/
#            <stratum>_species_pool_strict_full_func_clim_cov.csv (modelled species)
#          04_traits/outputs/<stratum>_traits.csv (SLA, seed mass, wood density, woodiness)
#          pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/
#            data_clean/species_climatic_affinity.csv (MAT_mean, MAP_mean)
# Outputs: 04_traits/outputs/trait_phylogenetic_signal.csv
#          04_traits/outputs/table_S2_phylogenetic_signal.csv / .docx
# =============================================================================

library(here)
library(ape)
library(phytools)
library(tidyverse)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
traits_dir <- here("pipeline", "01_data_prep", "04_traits")
tree_file  <- file.path(traits_dir, "outputs", "MORA_combined_phylo_tree.rda")   # combined tree, all species
pool_file <- function(s) here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                              paste0(s, "_species_pool_strict_full_func_clim_cov.csv"))
traits_by_stratum <- list(overstory  = c("MAT_mean", "MAP_mean", "seed_mass", "wood_density"),
                          understory = c("MAT_mean", "MAP_mean", "SLA", "seed_mass", "woodiness"))
log_traits <- c("SLA", "seed_mass")
nsim  <- 999
alpha <- 0.05
out_dir <- file.path(traits_dir, "outputs")

# Traits (04_traits) and climatic affinities, keyed on species with underscores
traits_tab <- function(s) read_csv(file.path(traits_dir, "outputs", paste0(s, "_traits.csv")), show_col_types = FALSE) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  select(-any_of(c("MAT_mean", "MAP_mean"))) %>% distinct(species, .keep_all = TRUE)
affinity <- read_csv(here("pipeline", "01_data_prep", "05_species_distribution", "02_integrate_gbif_chelsa",
                          "data_clean", "species_climatic_affinity.csv"), show_col_types = FALSE) %>%
  mutate(species = str_replace_all(species, " ", "_"),
         species = recode(species, Xanthocyparis_nootkatensis = "Callitropsis_nootkatensis")) %>%  # synonym used in the HMSC inputs
  select(species, MAT_mean, MAP_mean)

has_caper <- requireNamespace("caper", quietly = TRUE)
if (!has_caper) message("Package caper not installed: D statistic for woodiness skipped.")
set.seed(31)

# The .rda may hold the tree directly or inside a list (e.g. phylo.maker output, $scenario.3)
env <- new.env(); load(tree_file, envir = env)
objs <- mget(ls(env), envir = env)
find_phylo <- function(o) {
  if (inherits(o, "phylo")) return(o)
  if (is.list(o)) for (el in o) { r <- find_phylo(el); if (!is.null(r)) return(r) }
  NULL
}
tree_all <- NULL
for (o in objs) { tree_all <- find_phylo(o); if (!is.null(tree_all)) break }
if (is.null(tree_all)) stop("No phylo object found in ", tree_file, " (objects: ", paste(names(objs), collapse = ", "), ")")
tree_all$tip.label <- str_replace_all(tree_all$tip.label, " ", "_")
tree_all$node.label <- NULL                      # node labels can repeat tip names (not accepted by caper)

# -----------------------------------------------------------------------------
# 1. Phylogenetic signal per layer and trait
# -----------------------------------------------------------------------------
# phylosig() returns a list (with the test) or a bare number depending on the
# phytools version; read either form, and test by randomisation if no p-value
get_el <- function(res, el) {
  if (is.list(res)) return(as.numeric(res[[el]]))
  if (!is.null(names(res)) && el %in% names(res)) return(as.numeric(res[[el]]))
  if (el %in% c("K", "lambda")) return(as.numeric(res[1]))
  NA_real_
}

signal_stratum <- function(s) {
  pool <- read_csv(pool_file(s), show_col_types = FALSE) %>%
    transmute(species = str_replace_all(species, " ", "_")) %>%
    left_join(traits_tab(s), by = "species") %>%
    left_join(affinity, by = "species")
  miss_cols <- setdiff(traits_by_stratum[[s]], names(pool))
  if (length(miss_cols)) stop("Columns not found in ", basename(pool_file(s)), ": ", paste(miss_cols, collapse = ", "),
                              "\nAvailable: ", paste(names(pool), collapse = ", "))
  not_in_tree <- setdiff(pool$species, tree_all$tip.label)
  if (length(not_in_tree)) warning(s, ": species not in the tree, dropped: ", paste(not_in_tree, collapse = ", "))
  tree <- keep.tip(tree_all, intersect(pool$species, tree_all$tip.label))
  cat(s, ":", Ntip(tree), "species\n")
  
  map_dfr(traits_by_stratum[[s]], function(tr) {
    v <- setNames(pool[[tr]], pool$species)[tree$tip.label]
    if (tr == "woodiness") {
      if (!has_caper) return(tibble(stratum = s, trait = tr, type = "binary", n_species = length(v)))
      y <- as.integer(tolower(as.character(v)) %in% c("woody", "1", "true", "yes"))
      cd <- caper::comparative.data(tree, data.frame(species = names(v), y = y), names.col = "species")
      dres <- caper::phylo.d(cd, binvar = y, permut = nsim)
      return(tibble(stratum = s, trait = tr, type = "binary", n_species = length(v),
                    D = dres$DEstimate, p_D_random = dres$Pval1, p_D_brownian = dres$Pval0))
    }
    x <- as.numeric(v); names(x) <- names(v)
    if (tr %in% log_traits) x <- log10(x)
    x <- x[!is.na(x)]; tr_x <- keep.tip(tree, names(x))
    K <- phytools::phylosig(tr_x, x, method = "K", test = TRUE, nsim = nsim)
    L <- phytools::phylosig(tr_x, x, method = "lambda", test = TRUE)
    k <- get_el(K, "K"); pk <- get_el(K, "P"); l <- get_el(L, "lambda"); pl <- get_el(L, "P")
    perm_p <- function(stat, obs, n) {
      null <- replicate(n, get_el(phytools::phylosig(tr_x, setNames(sample(x), names(x)), method = stat), stat))
      (sum(null >= obs) + 1) / (n + 1)
    }
    if (is.na(pk)) pk <- perm_p("K", k, nsim)
    if (is.na(pl)) pl <- perm_p("lambda", l, 199)
    tibble(stratum = s, trait = tr, type = "continuous", n_species = length(x),
           K = k, p_K = pk, lambda = l, p_lambda = pl)
  })
}

sig <- map_dfr(names(traits_by_stratum), signal_stratum) %>%
  mutate(across(where(is.double), ~ signif(.x, 3)))
for (v in c("K", "p_K", "lambda", "p_lambda", "D", "p_D_random")) if (!v %in% names(sig)) sig[[v]] <- NA_real_
write_csv(sig, file.path(out_dir, "trait_phylogenetic_signal.csv"))
cat("\nPhylogenetic signal of traits:\n"); print(sig, n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 2. Table S2.2 (Word)
# -----------------------------------------------------------------------------
trait_labels <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity", seed_mass = "Seed mass",
                  wood_density = "Wood density", SLA = "SLA", woodiness = "Woodiness")
fmt_p <- function(p) case_when(is.na(p) ~ NA_character_, p < 0.001 ~ "< 0.001",
                               TRUE ~ formatC(p, format = "f", digits = 3))
fmt_est <- function(x, p) if_else(is.na(x), "–", paste0(formatC(x, format = "f", digits = 2), " (", fmt_p(p), ")"))

tab <- sig %>%
  filter(type == "continuous" | !is.na(D)) %>%
  mutate(Stratum = paste0(str_to_title(stratum), " (", n_species, " species)"),
         Trait = coalesce(trait_labels[trait], trait),
         K_cell = if_else(type == "binary", paste0("D = ", fmt_est(D, p_D_random)), fmt_est(K, p_K)),
         L_cell = if_else(type == "binary", "–", fmt_est(lambda, p_lambda)),
         K_bold = if_else(type == "binary", coalesce(p_D_random < alpha, FALSE), coalesce(p_K < alpha, FALSE)),
         L_bold = type != "binary" & coalesce(p_lambda < alpha, FALSE),
         order = match(trait, names(trait_labels))) %>%
  arrange(stratum, order)
# Layer label with the number of species of the layer (largest n among its traits)
tab <- tab %>% group_by(stratum) %>%
  mutate(Stratum = paste0(str_to_title(stratum), " (", max(n_species), " species)")) %>% ungroup()

out <- tab %>% transmute(Stratum, Trait, `Blomberg's K (p)` = K_cell, `Pagel's λ (p)` = L_cell)
write_csv(out, file.path(out_dir, "table_S2_phylogenetic_signal.csv"))

if (requireNamespace("flextable", quietly = TRUE) && requireNamespace("officer", quietly = TRUE)) {
  rule_thick <- officer::fp_border(color = "black", width = 1)
  rule_thin  <- officer::fp_border(color = "black", width = 0.5)
  stratum_ends <- which(out$Stratum != lead(out$Stratum, default = last(out$Stratum)))
  ft <- flextable::flextable(out) %>%
    flextable::merge_v(j = "Stratum") %>%
    flextable::border_remove() %>%
    flextable::hline_top(border = rule_thick, part = "header") %>%
    flextable::hline_bottom(border = rule_thin, part = "header") %>%
    flextable::hline(i = stratum_ends, border = rule_thin, part = "body") %>%
    flextable::hline_bottom(border = rule_thick, part = "body") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = "Stratum", part = "body") %>%
    flextable::bold(i = which(tab$K_bold), j = "Blomberg's K (p)", part = "body") %>%
    flextable::bold(i = which(tab$L_bold), j = "Pagel's λ (p)", part = "body") %>%
    flextable::align(align = "left", part = "all") %>%
    flextable::align(j = c("Blomberg's K (p)", "Pagel's λ (p)"), align = "center", part = "all") %>%
    flextable::valign(valign = "top", part = "body") %>%
    flextable::padding(padding.top = 2, padding.bottom = 2, part = "all") %>%
    flextable::width(j = names(out), width = c(1.8, 1.5, 1.5, 1.5)) %>%
    flextable::fix_border_issues()
  doc <- officer::read_docx() %>% flextable::body_add_flextable(ft)
  print(doc, target = file.path(out_dir, "table_S2_phylogenetic_signal.docx"))
  cat("Word table written to", out_dir, "\n")
} else {
  message("Install flextable and officer to write the Word table; the csv was written")
}
