# =============================================================================
# 31_phylogenetic_signal_traits.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run where the final fits are
#           stored; no refitting)
# Author: Agostina Torres
#
# Phylogenetic signal of the species traits used in the HMSC models, and the
# phylogenetic signal estimated by HMSC in the species responses.
#
# 1. Traits. For each layer, the trait matrix (m$Tr, i.e. the traits as they
#    enter the model, after the transformations in the trait formula) and the
#    phylogeny (m$phyloTree) are taken from the final occurrence fit.
#    Continuous traits: Blomberg's K and Pagel's lambda (phytools::phylosig),
#    with randomisation (K) and likelihood-ratio (lambda) tests. Binary traits
#    (woodiness): Fritz & Purvis D (caper::phylo.d), where D = 0 is the
#    expectation under Brownian motion and D = 1 a random distribution.
#
# 2. Responses. HMSC estimates rho, the proportion of the residual variation in
#    species responses (Beta) that follows the phylogeny after accounting for
#    traits (V x [rho C + (1 - rho) I]). Strong signal in the traits together
#    with a high rho means that the phylogenetic term and the trait effects can
#    share explanatory power; a low rho means that the phylogeny explains little
#    beyond the traits.
#
# 3. Table S2.2 (Word): phylogenetic signal of traits per layer; estimates
#    with significant signal (p < 0.05) in bold.
#
# Inputs : final fits (paths in section 0, as in 24_trait_effects_elevation.R)
# Outputs: outputs/phylo_signal/trait_phylogenetic_signal.csv
#          outputs/phylo_signal/hmsc_rho.csv
#          outputs/phylo_signal/table_S2_phylogenetic_signal.csv / .docx
# =============================================================================

library(here)
library(Hmsc)
library(ape)
library(phytools)
library(tidyverse)

# -----------------------------------------------------------------------------
# 0. Settings: final fits
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
fits <- tribble(
  ~stratum,     ~response,    ~path,
  "overstory",  "occurrence", file.path(module_dir, "outputs", "overstory_strict_medium_m1_time_elev_nf2_nsp_nos_pa.rds"),
  "overstory",  "abundance",  file.path(module_dir, "outputs", "overstory_strict_medium_m1_time_elev_nf2_nsp_nos_rel.rds"),
  "understory", "occurrence", file.path(module_dir, "outputs", "run_oct09", "understory_strict_long_m1_time_elev_nsp_nos_redpool_pa.rds"),
  "understory", "abundance",  file.path(module_dir, "outputs", "run_oct09", "understory_strict_long_m1_time_elev_nsp_nos_rel.rds")
)
nsim    <- 999
out_dir <- file.path(module_dir, "outputs", "phylo_signal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
missing <- fits$path[!file.exists(fits$path)]
if (length(missing) > 0) stop("Fits not found (edit section 0):\n", paste(missing, collapse = "\n"))
has_caper <- requireNamespace("caper", quietly = TRUE)
if (!has_caper) message("Package caper not installed: D statistic for binary traits skipped.")

set.seed(31)

# -----------------------------------------------------------------------------
# 1. Phylogenetic signal of traits (occurrence fit of each stratum)
# -----------------------------------------------------------------------------
# phylosig() returns a list (with the test) or a bare number depending on the
# phytools version and arguments; read either form
get_el <- function(res, el) {
  if (is.list(res)) return(as.numeric(res[[el]]))
  if (!is.null(names(res)) && el %in% names(res)) return(as.numeric(res[[el]]))
  if (el %in% c("K", "lambda")) return(as.numeric(res[1]))
  NA_real_
}

signal_traits <- function(m, stratum) {
  tree <- m$phyloTree
  Tr   <- m$Tr[, colnames(m$Tr) != "(Intercept)", drop = FALSE]
  rownames(Tr) <- m$spNames
  stopifnot("Tree tips and species differ" = setequal(tree$tip.label, rownames(Tr)))
  map_dfr(colnames(Tr), function(tr) {
    x <- setNames(Tr[, tr], rownames(Tr))[tree$tip.label]
    if (length(unique(x)) == 2) {                                   # binary trait
      if (!has_caper) return(tibble(stratum, trait = tr, type = "binary", n_species = length(x)))
      d  <- data.frame(species = names(x), y = as.integer(x == max(x)))
      tr_d <- tree; tr_d$node.label <- NULL                         # V.PhyloMaker node labels can repeat tip names
      cd <- caper::comparative.data(tr_d, d, names.col = "species")
      D  <- caper::phylo.d(cd, binvar = y, permut = nsim)
      tibble(stratum, trait = tr, type = "binary", n_species = length(x),
             D = D$DEstimate, p_D_random = D$Pval1, p_D_brownian = D$Pval0)
    } else {                                                        # continuous trait
      K <- phytools::phylosig(tree, x, method = "K", test = TRUE, nsim = nsim)
      L <- phytools::phylosig(tree, x, method = "lambda", test = TRUE)
      k <- get_el(K, "K"); pk <- get_el(K, "P"); l <- get_el(L, "lambda"); pl <- get_el(L, "P")
      # Fallback when no test is returned: randomisation of trait values among tips
      # (999 permutations for K, 199 for lambda, which needs a likelihood fit each)
      perm_p <- function(stat, obs, n) {
        null <- replicate(n, get_el(phytools::phylosig(tree, setNames(sample(x), names(x)), method = stat), stat))
        (sum(null >= obs) + 1) / (n + 1)
      }
      if (is.na(pk)) pk <- perm_p("K", k, nsim)
      if (is.na(pl)) pl <- perm_p("lambda", l, 199)
      tibble(stratum, trait = tr, type = "continuous", n_species = length(x),
             K = k, p_K = pk, lambda = l, p_lambda = pl)
    }
  })
}

trait_signal <- fits %>% filter(response == "occurrence") %>%
  pmap_dfr(function(stratum, response, path) signal_traits(readRDS(path), stratum)) %>%
  mutate(across(where(is.double), ~ signif(.x, 3)))
write_csv(trait_signal, file.path(out_dir, "trait_phylogenetic_signal.csv"))
cat("\nPhylogenetic signal of traits:\n"); print(trait_signal, n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# 2. Phylogenetic signal in species responses (rho) in each fit
# -----------------------------------------------------------------------------
# In this Hmsc version the posterior samples store the value of rho (on the grid
# m$rhopw[, 1]); older versions stored the grid index, which is converted here
rho_tab <- pmap_dfr(fits, function(stratum, response, path) {
  m   <- readRDS(path)
  rho <- unlist(lapply(poolMcmcChains(m$postList), `[[`, "rho"))
  if (all(rho >= 1 & rho == round(rho)) && max(rho) > 1) rho <- m$rhopw[rho, 1]
  tibble(stratum, response, n_species = m$ns, rho_mean = mean(rho),
         rho_q05 = unname(quantile(rho, 0.05)), rho_q95 = unname(quantile(rho, 0.95)),
         p_rho_above_0.5 = mean(rho > 0.5))
}) %>% mutate(across(where(is.double), ~ signif(.x, 3)))
write_csv(rho_tab, file.path(out_dir, "hmsc_rho.csv"))
cat("\nPhylogenetic signal in species responses (HMSC rho):\n"); print(rho_tab, width = Inf)

# -----------------------------------------------------------------------------
# 3. Table S2.2 (Word)
# -----------------------------------------------------------------------------
# Estimates with significant signal (p < 0.05) in bold. If woodiness was
# evaluated (Fritz & Purvis D), D and its p-value (test against a random
# distribution) are shown in the K column.
alpha <- 0.05

trait_labels <- c(MAT_mean = "Thermal affinity", MAP_mean = "Moisture affinity",
                  seed_mass = "Seed mass", wood_density = "Wood density",
                  SLA = "SLA", woodinesswoody = "Woodiness")

sig <- trait_signal
for (v in c("D", "p_D_random")) if (!v %in% names(sig)) sig[[v]] <- NA_real_

fmt_p <- function(p) case_when(is.na(p) ~ NA_character_, p < 0.001 ~ "< 0.001",
                               TRUE ~ formatC(p, format = "f", digits = 3))
fmt_est <- function(x, p) if_else(is.na(x), "–", paste0(formatC(x, format = "f", digits = 2), " (", fmt_p(p), ")"))

tab <- sig %>%
  filter(type == "continuous" | !is.na(D)) %>%                      # woodiness only if D was estimated
  mutate(Stratum = paste0(str_to_title(stratum), " (", n_species, " species)"),
         Trait = coalesce(trait_labels[trait], trait),
         K_cell = if_else(type == "binary", paste0("D = ", fmt_est(D, p_D_random)), fmt_est(K, p_K)),
         L_cell = if_else(type == "binary", "–", fmt_est(lambda, p_lambda)),
         K_bold = if_else(type == "binary", coalesce(p_D_random < alpha, FALSE), coalesce(p_K < alpha, FALSE)),
         L_bold = type != "binary" & coalesce(p_lambda < alpha, FALSE),
         order = match(trait, names(trait_labels))) %>%
  arrange(stratum, order)

out <- tab %>% transmute(Stratum, Trait, `Blomberg's K (p)` = K_cell, `Pagel's λ (p)` = L_cell)
write_csv(out, file.path(out_dir, "table_S2_phylogenetic_signal.csv"))
print(out, n = Inf)

# Word table, in the layout of the other tables: horizontal rules above and
# below the header, between layers and at the bottom; layer on the first row
# of each group; significant estimates in bold
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
