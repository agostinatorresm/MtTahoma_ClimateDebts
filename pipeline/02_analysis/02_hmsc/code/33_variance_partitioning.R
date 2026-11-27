# =============================================================================
# 33_variance_partitioning.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run where the final fits are
#           stored; no refitting)
# Author: Agostina Torres
#
# Variance partitioning of the final HMSC fits (Ovaskainen et al., 2017):
#   - variation in each species' occurrence or abundance attributed to
#     elevation (with the intercept), time (time and time x elevation) and the
#     plot-level random effect, averaged across species;
#   - proportion of the variation among species in their responses explained
#     by traits (R2T), for the response to time and overall.
# Together with the residual phylogenetic signal (rho; 31_phylogenetic_signal_traits.R)
# this summarises what the fixed effects, traits and phylogeny explain.
#
# Inputs : final fits (paths in section 0, as in 24_trait_effects_elevation.R)
# Outputs: outputs/variance_partitioning/variance_partitioning_species.csv
#          outputs/variance_partitioning/variance_partitioning_summary.csv
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
out_dir <- file.path(module_dir, "outputs", "variance_partitioning")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
missing <- fits$path[!file.exists(fits$path)]
if (length(missing) > 0) stop("Fits not found (edit section 0):\n", paste(missing, collapse = "\n"))

# Covariate groups: the intercept carries no variation among samples and is
# grouped with elevation; time includes its interaction with elevation
group_of <- function(cov) case_when(grepl("time", cov) ~ "Time", TRUE ~ "Elevation")

# -----------------------------------------------------------------------------
# 1. Variance partitioning per fit
# -----------------------------------------------------------------------------
res <- pmap(fits, function(stratum, response, path) {
  m <- readRDS(path)
  groupnames <- c("Elevation", "Time")
  group <- match(group_of(m$covNames), groupnames)
  vp <- computeVariancePartitioning(m, group = group, groupnames = groupnames)
  sp <- as_tibble(t(vp$vals), rownames = "species") %>%
    mutate(stratum = stratum, response = response, .before = 1)
  it <- which(m$covNames == "timepresent")
  r2t <- tibble(stratum = stratum, response = response,
                R2T_time = unname(vp$R2T$Beta[it]),                     # among-species variation in the time response explained by traits
                R2T_time_x_elev = unname(vp$R2T$Beta[which(m$covNames == "timepresent:elev_sc")]),
                R2T_overall = unname(vp$R2T$Y))                         # variation in species' predicted occurrence/abundance explained by traits
  rm(m); invisible(gc())
  list(sp = sp, r2t = r2t)
})
species_vp <- map_dfr(res, "sp")
write_csv(species_vp, file.path(out_dir, "variance_partitioning_species.csv"))

# Mean across species of each component, and R2T
summary_vp <- species_vp %>%
  pivot_longer(-c(stratum, response, species), names_to = "component", values_to = "prop") %>%
  group_by(stratum, response, component) %>%
  summarise(mean = mean(prop), min = min(prop), max = max(prop), .groups = "drop") %>%
  pivot_wider(names_from = component, values_from = c(mean, min, max)) %>%
  left_join(map_dfr(res, "r2t"), by = c("stratum", "response")) %>%
  mutate(across(where(is.double), ~ round(.x, 3)))
write_csv(summary_vp, file.path(out_dir, "variance_partitioning_summary.csv"))
cat("\nVariance partitioning (mean proportion across species; R2T = variation among species explained by traits):\n")
print(summary_vp, width = Inf)
