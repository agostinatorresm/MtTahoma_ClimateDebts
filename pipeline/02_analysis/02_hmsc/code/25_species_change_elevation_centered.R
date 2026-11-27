# =============================================================================
# 25_species_change_elevation_centered.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run where the final fits are
#           stored; no refitting)
# Author: Agostina Torres
#
# Understory change in relative abundance along the elevational gradient,
# relative to the mean change across the species present at each elevation
# ("centered" change), from the final HMSC abundance fit.
#
# Relative cover increases for most understory species because fewer species
# have cover in the present, so the raw change in abundance contains a common
# component that is not interpretable (see 23_species_change_centered.R).
# Centering is done within each posterior sample and at each elevation of the
# grid: species change minus the mean change across the species whose
# elevational range includes that elevation. The credible intervals and the
# posterior support of the centered change therefore carry the uncertainty of
# the common component.
#
# Change in log10 relative abundance at elevation e (fixed effects, as in
# 04_hmsc_results.R):
#     Beta[timepresent] + Beta[timepresent:elev_sc] * e
# The elevation grid and the species' elevational ranges are taken from the
# delta_elevation_species.csv table of 04, so the cells match the heat map.
#
# Inputs : final understory abundance fit (section 0)
#          pipeline/02_analysis/02_hmsc/outputs/
#            understory_strict_long_m1_time_elev_nsp_nos_redpool_rel_res_delta_elevation_species.csv (04)
# Outputs: pipeline/02_analysis/02_hmsc/outputs/centered/understory_abundance_change_elevation_centered.csv
# Figure : 05d_species_change_elevation_heatmap.R
# =============================================================================

library(here)
library(Hmsc)
library(tidyverse)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
fit_path   <- file.path(module_dir, "outputs", "run_oct09", "understory_strict_long_m1_time_elev_nsp_nos_rel.rds")
delta_file <- file.path(module_dir, "outputs", "understory_strict_long_m1_time_elev_nsp_nos_redpool_rel_res_delta_elevation_species.csv")
out_dir    <- file.path(module_dir, "outputs", "centered")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Fit not found (edit section 0)" = file.exists(fit_path),
          "delta_elevation_species.csv not found (edit section 0)" = file.exists(delta_file))

# -----------------------------------------------------------------------------
# 1. Elevation grid and species ranges (from 04)
# -----------------------------------------------------------------------------
delta <- read_csv(delta_file, show_col_types = FALSE) %>%
  filter(response %in% c("abundance", "log10_abundance"))
grid_tab <- delta %>% distinct(elev_sc, elevation) %>% arrange(elev_sc)
grid <- grid_tab$elev_sc
ranges <- delta %>% distinct(species, elev_min, elev_max, MAT_mean)

# -----------------------------------------------------------------------------
# 2. Change per posterior sample, species and elevation
# -----------------------------------------------------------------------------
m <- readRDS(fit_path)
it  <- which(m$covNames == "timepresent")
itx <- which(m$covNames == "timepresent:elev_sc")
stopifnot("Model has no time x elevation term" = length(it) == 1, length(itx) == 1)
other_time <- setdiff(grep("timepresent", m$covNames), c(it, itx))
stopifnot("Other time interactions in the model: the change is not Beta_t + Beta_tx * e" = length(other_time) == 0)

sp <- intersect(m$spNames, ranges$species)
stopifnot("Species names do not match the 04 table" = length(sp) == nrow(ranges))
rng <- ranges[match(sp, ranges$species), ]
# in_range[j, g]: grid point g within the elevational range of species j
in_range <- outer(rng$elev_min, grid, "<=") & outer(rng$elev_max, grid, ">=")

post <- poolMcmcChains(m$postList)
j_sp <- match(sp, m$spNames)    # Beta columns follow m$spNames (no dimnames in the posterior)
n_s <- length(post); n_j <- length(sp); n_g <- length(grid)
raw <- cen <- array(NA_real_, c(n_s, n_j, n_g))
for (s in seq_len(n_s)) {
  b0 <- post[[s]]$Beta[it,  j_sp]
  b1 <- post[[s]]$Beta[itx, j_sp]
  d  <- outer(b0, rep(1, n_g)) + outer(b1, grid)      # species x grid
  d[!in_range] <- NA
  raw[s, , ] <- d
  cen[s, , ] <- sweep(d, 2, colMeans(d, na.rm = TRUE))  # centered at each elevation
}
rm(m, post); gc()

# -----------------------------------------------------------------------------
# 3. Summaries
# -----------------------------------------------------------------------------
summ <- function(a, version) {
  map_dfr(seq_len(n_j), function(j) {
    x <- a[, j, ]
    keep <- in_range[j, ]
    tibble(species = sp[j], elev_sc = grid[keep],
           mean = colMeans(x[, keep, drop = FALSE]),
           q05 = apply(x[, keep, drop = FALSE], 2, quantile, 0.05),
           q25 = apply(x[, keep, drop = FALSE], 2, quantile, 0.25),
           q75 = apply(x[, keep, drop = FALSE], 2, quantile, 0.75),
           q95 = apply(x[, keep, drop = FALSE], 2, quantile, 0.95),
           support = colMeans(x[, keep, drop = FALSE] > 0))
  }) %>% mutate(version = version)
}
out <- bind_rows(summ(raw, "raw"), summ(cen, "centered")) %>%
  left_join(grid_tab, by = "elev_sc") %>%
  left_join(rng, by = "species") %>%
  mutate(stratum = "understory", response = "log10_abundance", .before = 1)

# Check: raw change reproduces the 04 table
chk <- out %>% filter(version == "raw") %>%
  inner_join(delta %>% select(species, elev_sc, mean_04 = mean), by = c("species", "elev_sc"))
cat("Raw change vs 04 table: r =", round(cor(chk$mean, chk$mean_04), 4),
    "; max abs difference =", signif(max(abs(chk$mean - chk$mean_04)), 3), "\n")

write_csv(out, file.path(out_dir, "understory_abundance_change_elevation_centered.csv"))
cat("\nCells with f >= 0.90:\n")
print(out %>% mutate(f = pmax(support, 1 - support)) %>%
        group_by(version) %>% summarise(prop_f90 = mean(f >= 0.90), prop_f80 = mean(f >= 0.80)))
cat("Table written to", out_dir, "\n")
