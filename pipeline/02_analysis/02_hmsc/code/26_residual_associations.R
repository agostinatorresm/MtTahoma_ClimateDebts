# =============================================================================
# 26_residual_associations.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run where the final fits are
#           stored; no refitting)
# Author: Agostina Torres
#
# Residual species associations at the plot level from the final HMSC fits:
# correlations between species in the plot-level random effect (latent
# factors), i.e. co-occurrence (occurrence models) or covariation in relative
# abundance (abundance models) not explained by time, elevation and their
# interaction. For each pair of species the posterior mean correlation and the
# posterior probability that it is positive (support) are computed with
# Hmsc::computeAssociations().
#
# Species orders used by the figure (05e_residual_associations_figure.R):
#   - cluster: hierarchical clustering (average linkage) of the posterior mean
#     correlation matrix, distance = 1 - correlation, so species with positive
#     residual associations are adjacent;
#   - each trait in the model's trait data (e.g. thermal and moisture
#     affinity, seed mass, SLA, wood density, woodiness), species ordered from
#     the lowest to the highest value.
#
# Inputs : final fits (paths in section 0, as in 23_species_change_centered.R)
# Outputs: pipeline/02_analysis/02_hmsc/outputs/associations/residual_associations.csv
#            (stratum, response, species pair, mean, support, f)
#          pipeline/02_analysis/02_hmsc/outputs/associations/species_orders.csv
#            (stratum, response, species, position in each order)
#          pipeline/02_analysis/02_hmsc/outputs/associations/species_traits.csv
#            (traits as used in the models, standardized)
# Figure : 05e_residual_associations_figure.R
# =============================================================================

library(here)
library(Hmsc)
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
random_level <- "plot"
out_dir <- file.path(module_dir, "outputs", "associations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
missing <- fits$path[!file.exists(fits$path)]
if (length(missing) > 0) stop("Fits not found (edit section 0):\n", paste(missing, collapse = "\n"))

# -----------------------------------------------------------------------------
# 1. Associations and species orders
# -----------------------------------------------------------------------------
assoc_all <- list(); orders_all <- list(); traits_all <- list()
for (i in seq_len(nrow(fits))) {
  st <- fits$stratum[i]; resp <- fits$response[i]
  m <- readRDS(fits$path[i])
  r <- which(m$rLNames == random_level)
  stopifnot("Random level not found in the model" = length(r) == 1)
  oc <- computeAssociations(m)[[r]]
  R <- oc$mean; S <- oc$support
  sp <- m$spNames
  dimnames(R) <- dimnames(S) <- list(sp, sp)
  
  # all pairs (upper triangle)
  idx <- which(upper.tri(R), arr.ind = TRUE)
  assoc_all[[paste(st, resp)]] <- tibble(
    stratum = st, response = resp,
    species_1 = sp[idx[, 1]], species_2 = sp[idx[, 2]],
    mean = R[idx], support = S[idx]) %>%
    mutate(f = pmax(support, 1 - support))
  
  # species orders: clustering of the mean correlations, and each trait
  hc <- hclust(as.dist(1 - R), method = "average")
  ord <- tibble(species = sp, cluster = match(sp, sp[hc$order]))
  tr <- m$TrData
  for (tn in names(tr)) {
    v <- tr[[tn]]
    if (is.numeric(v)) ord[[tn]] <- rank(v, ties.method = "first")
    else ord[[tn]] <- rank(as.numeric(factor(v)) + seq_along(v) / 1e6, ties.method = "first")
  }
  orders_all[[paste(st, resp)]] <- ord %>% mutate(stratum = st, response = resp, .before = 1)
  # traits as used in the model (standardized; factors as 0/1), for the inset of the figure
  traits_all[[paste(st, resp)]] <- as_tibble(tr) %>%
    mutate(across(where(Negate(is.numeric)), ~ as.numeric(factor(.x)) - 1)) %>%
    mutate(stratum = st, species = sp, .before = 1)
  
  n_pairs <- nrow(idx)
  cat(sprintf("%s %s: %d species, %d pairs; f >= 0.80: %d positive, %d negative; f >= 0.90: %d positive, %d negative\n",
              st, resp, length(sp), n_pairs,
              sum(S[idx] >= 0.80), sum(S[idx] <= 0.20), sum(S[idx] >= 0.90), sum(S[idx] <= 0.10)))
  rm(m, oc); gc()
}

assoc  <- bind_rows(assoc_all)
orders <- bind_rows(orders_all)
write_csv(assoc,  file.path(out_dir, "residual_associations.csv"))
write_csv(orders, file.path(out_dir, "species_orders.csv"))
write_csv(bind_rows(traits_all) %>% distinct(stratum, species, .keep_all = TRUE),
          file.path(out_dir, "species_traits.csv"))
cat("\nOrders available:", paste(setdiff(names(orders), c("stratum", "response", "species")), collapse = ", "), "\n")
cat("Tables written to", out_dir, "(figure: 05e_residual_associations_figure.R)\n")
