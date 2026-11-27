# =============================================================================
# 01_prepare_hmsc_inputs.R
# Author: Agostina Torres
#
# Builds and validates all inputs for the HMSC models of one vegetation layer
# (overstory or understory) and defines the unfitted model objects.
#
#   Y   : sample x species matrices
#         - Y_pa  : occurrence (1/0); for the understory this includes species
#                   recorded at the whole-plot scale without a cover value
#         - Y_list: one matrix per abundance response, log10 abundance
#                   conditional on presence (NA where absent or where the
#                   species was recorded without an abundance value), absolute
#                   or relative to the sample total (see responses_by_stratum)
#   Tr  : species x trait table (TrData), rows in the same order as Y columns
#   X   : sample-level covariates (XData), rows in the same order as Y rows
#   Phy : phylogeny pruned to the species in Y
#   studyDesign / random levels:
#         - sample : one unit per plot x survey period (residual species
#                    associations within a sampling occasion)
#         - plot   : permanent plot, spatially explicit (projected coordinates);
#                    links the historical and contemporary records of the same
#                    plot and accounts for spatial autocorrelation
#         The number of latent factors and whether the plot level is spatial
#         are set per random-level structure (rl_structures).
#
# Sample unit: sample_id = plot_id x period. Only plots surveyed in both periods
# and with coordinates, elevation and climate in both periods are retained.
#
# Fixed-effect structures (one model per structure and response type):
#   m0_time            : ~ time
#                        Change within plots averaged over the gradient; the
#                        baseline distribution of each species along the
#                        gradient is absorbed by the plot-level (spatial)
#                        random effect. Trait effects on the time response
#                        (Gamma) test whether change depends on traits.
#   m1_time_elev       : ~ time * elev_sc
#                        As m0, with the baseline and the change allowed to
#                        vary linearly with elevation.
#   m2_climchange      : ~ elev_sc + d_temp_JJA_sc + d_MAP_sc
#   m3_climchange_elev : ~ (d_temp_JJA_sc + d_MAP_sc) * elev_sc
#
# Climate-change predictors are the change experienced by a plot since the
# historical survey: 0 for historical samples and (present - past) for
# contemporary samples. Climate change is two-dimensional (temperature and
# precipitation axes, nearly orthogonal), so one variable per axis enters the
# models (summer temperature, temp_JJA, and annual precipitation, MAP). Predictors are scaled but NOT centred, so that
# 0 keeps its meaning of "no change since the historical survey". A single
# PCA-based axis is available as an alternative.
#
# Elevation is centred and scaled explicitly (used in interactions). Traits are
# log-transformed where required and scaled explicitly before entering TrData.
#
# Outputs (pipeline/02_analysis/02_hmsc/):
#   data_clean/<stratum>_<pool>_hmsc_inputs.rds            (Y, XData, TrData, tree, ...)
#   data_clean/<stratum>_<pool>_hmsc_models_unfitted.rds
#   outputs/<stratum>_<pool>_check_*.csv                   (input check tables)
# All file names start with the layer (overstory / understory) and pool
# (plus the species subset for sensitivity runs, e.g. understory_strict_nonpatchy).
# =============================================================================

library(here)
library(tidyverse)
library(Hmsc)
library(ape)
library(sf)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
stratum <- "overstory"   # "overstory" | "understory"
pool    <- "strict"      # "strict" | "broad" (pool file name pattern below)

# Traits entering TrFormula, per layer
traits_model <- list(
  overstory  = c("MAT_mean", "MAP_mean", "seed_mass", "wood_density"),
  understory = c("MAT_mean", "MAP_mean", "SLA", "seed_mass", "woodiness")
)
log_traits         <- c("SLA", "seed_mass")   # log10-transformed before scaling
categorical_traits <- c("woodiness")

# Climate-change predictors. Climate change between surveys is
# two-dimensional (Figure S1.4B): a temperature axis (delta MAT and related
# variables) and a precipitation axis (delta MAP and related variables), nearly
# orthogonal to each other. One representative variable per axis enters the
# models as a separate predictor.
#   clim_change_source = "variables": one predictor per variable in
#                        clim_change_vars (delta = present - past)
#   clim_change_source = "pca": a single predictor, the shift along PC1 of a
#                        PCA of clim_vars fitted on historical conditions
clim_change_source <- "variables"
clim_change_vars   <- c("temp_JJA", "MAP")   # summer temperature and annual precipitation

# Variables used by the PCA option (columns <variable>_past / _present in
# climate_past_present_wide.csv); spring/autumn omitted as largely redundant
# with the annual means.
clim_vars          <- c("MAT", "MAT_min", "MAP", "PAS",
                        "temp_DJF", "temp_JJA", "precip_DJF", "precip_JJA")
clim_temp_var      <- "MAT"   # used to orient PC1 so that positive = warmer

# Handling of empty samples after subsetting to the species pool: plots with
# no pool species in both periods are always dropped. With
# drop_plots_with_empty_sample = TRUE, plots in which either survey contains
# no pool species are also dropped (both samples), keeping only plots with at
# least one pool species in each survey.
drop_plots_with_empty_sample <- TRUE
drop_species_missing_traits <- TRUE   # FALSE = stop if any species lacks a trait
min_abu_obs                 <- 5      # minimum abundance observations per species (abundance models)
# Understory paired abundance (abu, rel): a species is kept when it has a cover
# value in both periods of at least one plot (2 values); species with fewer
# than min_abu_obs values are kept but flagged, and their individual
# responses are not interpreted (their estimates rest mainly on the trait
# and phylogenetic structure of the model)
min_abu_obs_paired          <- 2

# Models and response types to build (subset of the three structures defined
# in section 6: m0_time, m1_time_elev, m2_climchange, m3_climchange_elev)
models_to_build    <- c("m0_time", "m1_time_elev")
# Response types per layer. "pa" is occurrence (probit). Every other
# response is log10 abundance conditional on presence (normal), of the
# variable given in response_vars: absolute (names starting with "abu") or
# relative to the sample total (names starting with "rel"). Response names
# contain no underscore (they end the model id).
#   overstory  abu    : basal area as % of plot area
#              rel    : relative basal area (species basal area / sample total)
#              abucnt : stem density (stems per m2 of plot)
#              relcnt : relative density (species stems / sample total)
#   understory abu    : percent cover, paired by plot (main; see below)
#              rel    : relative cover, paired by plot (main)
#              abuall : percent cover, all presences with cover (check)
#              relall : relative cover, all presences with cover (check)
# All responses are fitted; relative basal area and relative cover are the
# main abundance responses (see "Relative abundance" below).
responses_by_stratum <- list(
  overstory  = c("pa", "abu", "rel", "abucnt", "relcnt"),
  understory = c("pa", "abu", "rel", "abuall", "relall")
)
response_vars <- list(
  overstory  = c(abu = "basal_area_perc", rel = "basal_area_perc",
                 abucnt = "density", relcnt = "count"),
  understory = c(abu = "cover", rel = "cover", abuall = "cover", relall = "cover")
)
is_relative <- function(r) startsWith(r, "rel")

spatial_method <- "Full"   # "Full" (~330 plots is feasible) | "NNGP"
projected_crs  <- 32610    # UTM zone 10N (metres)

# Random-level structures. Each structure sets the maximum number of latent
# factors per random level (nf_max), whether the plot level is spatially
# explicit, and whether a sample-level (plot x survey) random level is
# included. One unfitted model is built per fixed-effect structure x response
# x random-level structure. The "default" structure keeps the plain model id
# (e.g. m1_time_elev_pa); other structures add their name to the id
# (e.g. m1_time_elev_nf2_pa).
#   default     : nf_max = 5, spatial plot level, sample level
#   nf2         : nf_max = 2, spatial plot level, sample level (fewer latent
#                 factors; relevant when the number of species is small
#                 relative to 2 x nf_max)
#   nf2_nsp     : nf_max = 2, non-spatial plot level (plot identity only),
#                 sample level
#   nsp         : nf_max = 5, non-spatial plot level, sample level
#   nos         : nf_max = 5, spatial plot level, no sample level
#   nf2_nos     : nf_max = 2, spatial plot level, no sample level
#   nsp_nos     : nf_max = 5, non-spatial plot level, no sample level
#   nf2_nsp_nos : nf_max = 2, non-spatial plot level, no sample level
# Without the sample level, the change between surveys is carried entirely by
# the fixed effects (and by the plot-level effect only through its persistence
# across surveys), and residual species associations are estimated at the plot
# level only.
rl_structures <- list(
  default     = list(nf_max = 5, spatial = TRUE,  sample = TRUE),
  nf2         = list(nf_max = 2, spatial = TRUE,  sample = TRUE),
  nf2_nsp     = list(nf_max = 2, spatial = FALSE, sample = TRUE),
  nsp         = list(nf_max = 5, spatial = FALSE, sample = TRUE),
  nos         = list(nf_max = 5, spatial = TRUE,  sample = FALSE),
  nf2_nos     = list(nf_max = 2, spatial = TRUE,  sample = FALSE),
  nsp_nos     = list(nf_max = 5, spatial = FALSE, sample = FALSE),
  nf2_nsp_nos = list(nf_max = 2, spatial = FALSE, sample = FALSE)
)
structures_to_build <- names(rl_structures)

# Column names in the input tables (new name = name in file)
cols_veg       <- c(plot_id = "plot_id", period = "period", species = "species")

# Relative abundance. Dividing each species' value by the total of the sample
# removes differences between surveys that act on all species of a sample
# alike; on the log10 scale this is a sample-level offset shared by all
# species. Relative values express change in the share of the community,
# which is the quantity behind community thermophilization.
#   understory: percent cover is exposed to observer, phenological timing and
#     the change from whole-plot estimates in the historical survey to four
#     1-m2 subplots in the present one (absolute cover declined by ~60% in
#     most plots).
#   overstory: total basal area and stem number change with stand development
#     (growth, recruitment, self-thinning), which would load on the time
#     effect of every species. Relative basal area is dominated by large
#     trees; relative density weights every stem equally and responds earlier
#     to recruitment, but depends on the same minimum stem diameter being
#     tallied in both surveys.
# Relativisation does not remove species-specific biases (identification,
# detection of patchy species in small subplots, stem-size thresholds).
# Denominator: total of the response variable over the species recorded in
# both surveys somewhere in the layer (before restriction to the analysis
# pool, so that it does not depend on pool or modelling choices). Species
# recorded in one survey only are left out of the denominator: in the
# overstory, several of them (e.g. Acer circinatum, Taxus brevifolia, Alnus
# alnobetula) appear with many small stems in the present survey only, which
# points to a difference in which stems were tallied rather than to
# colonisation, and would lower the relative values of all other species in
# the present.
# Overstory basal area is expressed as a percentage of plot area
# (basal_area_spp_m2 / plot_area x 100) and stem density as stems per m2 of
# plot, so that both are comparable across plots of different size. With
# YScale = TRUE the unit does not affect the fitted models.

# Species excluded from all models of a layer, with the reason documented
# here and in the Supplementary Methods (S.1).
species_exclude <- list(
  overstory  = character(0),
  understory = character(0)
)

cols_plot_info <- c(plot_id = "plot_id", lat = "lat", long = "long",
                    elevation = "elevation")
cols_climate   <- c(plot_id = "plot_id")

# Sensitivity subset (understory; Supplementary Methods S3). Present-day
# understory data come from four 1-m2 subplots (cover) plus a 50-m2 census
# (presence only), whereas historical data cover the whole plot. Species
# detected in few subplots where they occur (mean subplot occupancy below
# `patchy_threshold`) are less reliably characterised by the present-day
# design. species_subset = "nonpatchy" removes them, so that the main
# conclusions can be compared between the full and the reduced species set.
# The subset name is added to the file tag (e.g. understory_strict_nonpatchy).
species_subset   <- "all"    # "all" | "nonpatchy"
patchy_threshold <- 0.4
cols_reliability <- c(species = "species", occupancy = "mean_subplot_occupancy")
period_codes   <- c(past = "past", present = "present")   # values used in `period`

# Understory abundance. Detection differs between the past (whole plot,
# 250-1000 m2) and the present (cover in four 1-m2 subplots, presence also in
# a 50-m2 census). Abundance is therefore built in two ways, both restricted
# to the modelled species (non-modelled species do not enter the totals):
# Relative cover is computed on the cover matrix of the modelled species
# (samples x modelled species): species cover / row total, i.e. the total
# cover of the modelled species with a cover value in that plot and period.
# The missing values below are applied afterwards and do not change the
# totals.
#   abu, rel (main): paired by plot. Only species with a cover value in both
#     the past and the present of the same plot keep their values; all other
#     plot x species pairs are missing (local losses and gains, which are part
#     of the occurrence model, and occurrences recorded only in the 50-m2
#     census).
#   abuall, relall (check): every presence with a cover value. Pairs whose
#     present occurrence was recorded only
#     in the 50-m2 census are set to missing in the past as well
#     (census_only_past_na), so that the past does not include low-cover
#     occurrences that the present design misses.
census_only_past_na <- TRUE

# Input files: NULL = located by file name anywhere in pipeline/ (tmp_files
# excluded); set an explicit here() path if several copies exist.
input_paths <- list(
  reliability  = NULL,   # composite_reliability_diagnostic.csv (sensitivity subset only)
  veg          = NULL,   # <stratum>_long.csv
  species_pool = here("pipeline", "01_data_prep", "02_species_pool", "outputs",
                      paste0(stratum, "_species_pool_", pool, "_full_func_clim_cov.csv")),
  plot_info    = NULL,   # plot_info.csv
  climate      = here("pipeline", "01_data_prep", "03_climate", "02_explore_climate",
                      "data_clean", "climate_past_present_wide.csv"),
  traits       = NULL,   # <stratum>_traits.csv
  affinity     = NULL,   # species_climatic_affinity.csv
  phylogeny    = here("pipeline", "01_data_prep", "04_traits", "outputs",
                      "MORA_combined_phylo_tree.rda")   # phy + C from 06_build_phylogeny.R
)

# -----------------------------------------------------------------------------
# Paths and helpers
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
if (!species_subset %in% c("all", "nonpatchy")) stop("Unknown species_subset: ", species_subset)
if (species_subset != "all" && stratum != "understory") {
  stop("species_subset = '", species_subset, "' applies to the understory only.")
}
# File tag: <stratum>_<pool>, plus the subset name for sensitivity subsets.
# Downstream scripts (01b, 02, 03) use the same tag through their `pool`
# setting, e.g. pool <- "strict_nonpatchy".
pool_label <- if (species_subset == "all") pool else paste(pool, species_subset, sep = "_")
run_tag    <- paste(stratum, pool_label, sep = "_")
clean_dir  <- file.path(module_dir, "data_clean")
out_dir    <- file.path(module_dir, "outputs")
walk(c(clean_dir, out_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

# File names always start with <stratum>_<pool>
check_file <- function(name) file.path(out_dir, paste0(run_tag, "_check_", name))
clean_file <- function(name) file.path(clean_dir, paste0(run_tag, "_", name))

find_input <- function(file_name, root = here("pipeline")) {
  pattern <- paste0("^", str_replace_all(file_name, "\\.", "\\\\."), "$")
  hits <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  hits <- hits[!str_detect(hits, "/tmp_files/")]
  if (length(hits) == 0) {
    # Suggest files sharing the extension and the first word of the name
    stem  <- str_extract(file_name, "^[A-Za-z]+")
    ext   <- str_extract(file_name, "\\.[A-Za-z]+$")
    near  <- list.files(root, pattern = paste0(stem, ".*", str_replace(ext, "\\.", "\\\\."), "$"),
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    near  <- near[!str_detect(near, "/tmp_files/")]
    stop("File not found in pipeline/: ", file_name,
         "\nSet its path explicitly in `input_paths`.",
         if (length(near) > 0) paste0("\nCandidates:\n", paste(head(near, 30), collapse = "\n")))
  }
  if (length(hits) > 1) {
    stop("Several copies of ", file_name, " found; set the path in `input_paths`:\n",
         paste(hits, collapse = "\n"))
  }
  hits
}

resolve_path <- function(key, file_name) {
  if (!is.null(input_paths[[key]])) input_paths[[key]] else find_input(file_name)
}

check_cols <- function(df, cols, df_name) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(df_name, " is missing columns: ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(names(df), collapse = ", "),
         "\nUpdate the column mapping in Settings.")
  }
}

# Species key used for all matrix names (Y columns, TrData rows, tree tips)
sp_key <- function(x) str_replace_all(str_squish(x), " ", "_")

# Synonyms between GBIF-based affinity names and project names
affinity_synonyms <- c(
  "Xanthocyparis nootkatensis" = "Callitropsis nootkatensis",
  "Brewerimitella breweri"     = "Mitella breweri",
  "Mitellastra caulescens"     = "Mitella caulescens",
  "Pectiantia pentandra"       = "Mitella pentandra",
  "Dodecatheon jeffreyanum"    = "Primula jeffreyi"
)

# -----------------------------------------------------------------------------
# 1. Load inputs
# -----------------------------------------------------------------------------
veg_raw <- read_csv(resolve_path("veg", paste0(stratum, "_long.csv")),
                    show_col_types = FALSE)
check_cols(veg_raw, cols_veg, "Vegetation table")
veg_raw <- veg_raw %>% rename(all_of(cols_veg))

# Species pool restricted to species with full functional-trait and
# climatic-affinity coverage (one file per stratum)
species_pool <- read_csv(input_paths$species_pool, show_col_types = FALSE)
check_cols(species_pool, "species", "Species pool")
if ("stratum" %in% names(species_pool)) {
  species_pool <- species_pool %>% filter(stratum == !!stratum)
}
pool_species <- species_pool %>% distinct(species) %>% pull(species) %>% sp_key()
message(length(pool_species), " species in the ", stratum, " ", pool, " pool.")

# Sensitivity subset: drop species flagged as patchy in the present-day design
if (species_subset == "nonpatchy") {
  reliability <- read_csv(resolve_path("reliability", "composite_reliability_diagnostic.csv"),
                          show_col_types = FALSE)
  check_cols(reliability, cols_reliability, "Reliability table")
  reliability <- reliability %>%
    rename(all_of(cols_reliability)) %>%
    mutate(species = sp_key(species)) %>%
    distinct(species, occupancy)

  no_flag <- setdiff(pool_species, reliability$species)
  if (length(no_flag) > 0) {
    stop("Pool species without a subplot-occupancy value: ", paste(no_flag, collapse = ", "))
  }
  patchy_species <- reliability %>%
    filter(species %in% pool_species, occupancy < patchy_threshold) %>%
    pull(species)
  write_csv(tibble(species = patchy_species), check_file("patchy_species_excluded.csv"))
  pool_species <- setdiff(pool_species, patchy_species)
  message(length(patchy_species), " patchy species excluded (mean subplot occupancy < ",
          patchy_threshold, "); ", length(pool_species), " species retained.")
}

plot_info_raw <- read_csv(resolve_path("plot_info", "plot_info.csv"), show_col_types = FALSE)
check_cols(plot_info_raw, cols_plot_info, "plot_info")
plot_info_raw <- plot_info_raw %>% rename(all_of(cols_plot_info))

# Climate table: one row per plot with past and present values per variable
# (e.g. tas_past / tas_present); reshaped to one row per plot x period.
climate_wide <- read_csv(input_paths$climate, show_col_types = FALSE)
check_cols(climate_wide, cols_climate[["plot_id"]], "Climate table")
climate_wide <- climate_wide %>% rename(plot_id = all_of(cols_climate[["plot_id"]]))

clim_read_vars <- union(clim_vars, clim_change_vars)
clim_pattern <- paste0("^(", paste(clim_read_vars, collapse = "|"), ")_(past|present)$")
clim_cols    <- str_subset(names(climate_wide), clim_pattern)
missing_clim <- setdiff(as.vector(outer(clim_read_vars, c("past", "present"), paste, sep = "_")),
                        clim_cols)
if (length(missing_clim) > 0) {
  stop("Climate table lacks columns: ", paste(missing_clim, collapse = ", "),
       "
Available columns: ", paste(names(climate_wide), collapse = ", "),
       "
Update `clim_vars` (expected <variable>_past / <variable>_present).")
}

climate_raw <- climate_wide %>%
  select(plot_id, all_of(clim_cols)) %>%
  pivot_longer(-plot_id, names_to = c(".value", "period"), names_pattern = clim_pattern)

traits_raw <- read_csv(resolve_path("traits", paste0(stratum, "_traits.csv")),
                       show_col_types = FALSE)

# Phylogeny: Newick (.tre/.nwk), .rds or .rda/.RData containing a "phylo" object
read_phylo <- function(path) {
  ext <- tolower(tools::file_ext(path))
  obj <- switch(ext,
    rds   = readRDS(path),
    rda   = ,
    rdata = { env <- new.env(); load(path, envir = env)
              Filter(function(x) inherits(x, "phylo"), as.list(env))[[1]] },
    read.tree(path))
  if (!inherits(obj, "phylo")) {
    # e.g. V.PhyloMaker2 output list with the tree in $scenario.3
    obj <- Filter(function(x) inherits(x, "phylo"), obj)[[1]]
  }
  if (is.null(obj) || !inherits(obj, "phylo")) stop("No phylo object found in ", path)
  obj
}

if (is.null(input_paths$phylogeny)) {
  stop("Set input_paths$phylogeny to the tree built by 06_build_phylogeny.R.\n",
       "Files with 'phylo' or 'tree' in their name:\n",
       paste(list.files(here("pipeline"), pattern = "phylo|tree", recursive = TRUE,
                        ignore.case = TRUE), collapse = "\n"))
}
tree_raw <- read_phylo(input_paths$phylogeny)

# Period coding
unexpected_periods <- setdiff(unique(veg_raw$period), period_codes)
if (length(unexpected_periods) > 0) {
  stop("Unexpected period values in vegetation table: ",
       paste(unexpected_periods, collapse = ", "), "\nUpdate `period_codes`.")
}
recode_period <- function(x) factor(names(period_codes)[match(x, period_codes)],
                                    levels = c("past", "present"))

veg <- veg_raw %>%
  mutate(period = recode_period(period), species = sp_key(species))

# Plot x period combinations surveyed, taken before removing absence rows so
# that a survey with no species recorded still counts as surveyed
veg_surveyed <- veg %>% distinct(plot_id, period)

# Only rows recording a presence (occurrence > 0, when the column exists)
if ("occurrence" %in% names(veg)) {
  veg <- veg %>% filter(occurrence > 0)
}

if (stratum == "overstory") {
  check_cols(veg, c("basal_area_spp_m2", "count", "plot_area"), "Overstory vegetation table")
  if (any(is.na(veg$plot_area) | veg$plot_area <= 0)) {
    stop("Missing or non-positive plot_area values in the overstory table.")
  }
  veg <- veg %>% mutate(basal_area_perc = basal_area_spp_m2 / plot_area * 100,
                        density         = count / plot_area)
}

abundance_responses <- setdiff(responses_by_stratum[[stratum]], "pa")
resp_var <- response_vars[[stratum]][abundance_responses]
check_cols(veg, unique(resp_var), "Vegetation table")

# Denominator of each relative response: total of its variable per sample over
# the species recorded with a positive value in both surveys
relative_totals <- map(set_names(keep(abundance_responses, is_relative)), function(r) {
  v <- veg %>% transmute(plot_id, period, species, value = .data[[resp_var[[r]]]]) %>%
    filter(!is.na(value), value > 0)
  both <- v %>% distinct(species, period) %>% count(species) %>% filter(n == 2) %>% pull(species)
  list(species_both = both,
       totals = v %>% filter(species %in% both) %>%
         group_by(plot_id, period) %>% summarise(total = sum(value), .groups = "drop"))
})

climate <- climate_raw %>%
  mutate(period = factor(period, levels = c("past", "present")))

# One row per sample x species is required
dup_rows <- veg %>% count(plot_id, period, species) %>% filter(n > 1)
if (nrow(dup_rows) > 0) {
  write_csv(dup_rows, check_file("veg_duplicated_sample_species.csv"))
  stop(nrow(dup_rows), " duplicated plot x period x species rows in the vegetation table ",
       "(see the check_veg_duplicated_sample_species file in outputs/); the table must be aggregated to one ",
       "row per sample and species.")
}

# -----------------------------------------------------------------------------
# 2. Plots and samples
# -----------------------------------------------------------------------------
# Plot-level coordinates and elevation: past and present GPS positions differ
# slightly for the same plot, so they are averaged per plot_id.
plot_meta <- plot_info_raw %>%
  group_by(plot_id) %>%
  summarise(lat = mean(lat, na.rm = TRUE), long = mean(long, na.rm = TRUE),
            elevation = mean(elevation, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(lat), !is.na(long), !is.na(elevation))

plots_veg_both <- veg_surveyed %>%
  count(plot_id) %>%
  filter(n == 2) %>%
  pull(plot_id)

plots_climate_both <- climate %>%
  filter(if_all(all_of(clim_read_vars), ~ !is.na(.x))) %>%
  distinct(plot_id, period) %>%
  count(plot_id) %>%
  filter(n == 2) %>%
  pull(plot_id)

plot_filter <- tibble(plot_id = unique(veg_surveyed$plot_id)) %>%
  mutate(veg_both_periods     = plot_id %in% plots_veg_both,
         has_coords_elevation = plot_id %in% plot_meta$plot_id,
         climate_both_periods = plot_id %in% plots_climate_both,
         kept = veg_both_periods & has_coords_elevation & climate_both_periods)

write_csv(plot_filter, check_file("plot_filter.csv"))
plot_filter %>% filter(!kept)

plots_kept <- plot_filter %>% filter(kept) %>% pull(plot_id)

samples <- expand_grid(plot_id = sort(plots_kept),
                       period  = factor(c("past", "present"), levels = c("past", "present"))) %>%
  mutate(sample_id = paste(plot_id, period, sep = "_"))

# -----------------------------------------------------------------------------
# 3. Y matrices
# -----------------------------------------------------------------------------
missing_pool_species <- setdiff(pool_species, unique(veg$species))
if (length(missing_pool_species) > 0) {
  stop("Species in the pool not found in the vegetation table: ",
       paste(missing_pool_species, collapse = ", "))
}

excluded_species <- intersect(species_exclude[[stratum]], pool_species)
if (length(excluded_species) > 0) {
  message("Species excluded by `species_exclude`: ", paste(excluded_species, collapse = ", "))
  pool_species <- setdiff(pool_species, excluded_species)
}

veg_pool <- veg %>%
  filter(species %in% pool_species) %>%
  inner_join(samples, by = c("plot_id", "period"))

# The vegetation table holds presences only. Absences are completed with 0 by
# crossing every retained sample with every pool species.
# Presence: any record of the species in the sample (whole-plot records
# without cover included). Abundance: only where a positive value exists;
# species present without an abundance value get NA (not 0), and absences get
# NA (abundance is modelled conditional on presence).
veg_grid <- samples %>%
  select(sample_id, plot_id, period) %>%
  expand_grid(species = sort(pool_species)) %>%
  left_join(veg_pool %>% select(plot_id, period, species, all_of(unique(resp_var))) %>%
              mutate(recorded = TRUE),
            by = c("plot_id", "period", "species")) %>%
  mutate(pa = as.integer(!is.na(recorded)))

if (stratum == "overstory") {
  for (r in abundance_responses) {
    v <- veg_grid[[resp_var[[r]]]]
    ok <- veg_grid$pa == 1 & !is.na(v) & v > 0
    if (is_relative(r)) {
      tot <- veg_grid %>% select(plot_id, period) %>%
        left_join(relative_totals[[r]]$totals, by = c("plot_id", "period")) %>% pull(total)
      ok <- ok & veg_grid$species %in% relative_totals[[r]]$species_both
      veg_grid[[paste0("obs_", r)]] <- if_else(ok, v / tot, NA_real_)
    } else {
      veg_grid[[paste0("obs_", r)]] <- if_else(ok, v, NA_real_)
    }
  }
} else {
  # Understory: see "Understory abundance" in the settings
  veg_grid <- veg_grid %>%
    mutate(has_cover = pa == 1 & !is.na(cover) & cover > 0) %>%
    group_by(plot_id, species) %>%
    mutate(paired = sum(has_cover) == 2) %>%
    group_by(sample_id) %>%
    mutate(total_all = sum(if_else(has_cover, cover, 0))) %>%
    ungroup() %>%
    mutate(obs_abu    = if_else(has_cover & paired, cover, NA_real_),
           obs_rel    = if_else(has_cover & paired, cover / total_all, NA_real_),
           obs_abuall = if_else(has_cover, cover, NA_real_),
           obs_relall = if_else(has_cover, cover / total_all, NA_real_))

  census_only <- veg_grid %>%
    filter(period == "present", pa == 1, !has_cover) %>%
    distinct(plot_id, species) %>%
    mutate(census_only = TRUE)
  veg_grid <- veg_grid %>%
    left_join(census_only, by = c("plot_id", "species")) %>%
    mutate(census_only = coalesce(census_only, FALSE))
  if (census_only_past_na) {
    past_dropped <- veg_grid %>%
      filter(period == "past", census_only, !is.na(obs_abuall)) %>%
      count(species, name = "n_past_values_set_missing")
    write_csv(past_dropped, check_file("census_only_past_missing.csv"))
    message(sum(past_dropped$n_past_values_set_missing),
            " past values set to missing in abuall/relall (present occurrence recorded only in the 50-m2 census)")
    veg_grid <- veg_grid %>%
      mutate(across(c(obs_abuall, obs_relall), ~ if_else(period == "past" & census_only, NA_real_, .x)))
  }
  paired_summary <- veg_grid %>%
    filter(period == "past") %>%
    group_by(species) %>%
    summarise(n_plots_present_both = sum(pa == 1 & paired),
              n_plots_cover_both   = sum(paired), .groups = "drop")
  write_csv(paired_summary, check_file("abundance_paired_plots.csv"))
  message("Understory abundance paired by plot: ", sum(!is.na(veg_grid$obs_abu)), " values (",
          sum(!is.na(veg_grid$obs_abuall)), " with all presences)")
  veg_grid <- select(veg_grid, -has_cover, -paired, -total_all, -census_only)
}

to_matrix <- function(df, value_col) {
  df %>%
    select(sample_id, species, all_of(value_col)) %>%
    pivot_wider(names_from = species, values_from = all_of(value_col)) %>%
    column_to_rownames("sample_id") %>%
    as.matrix()
}

Y_pa   <- to_matrix(veg_grid, "pa")
Y_list <- map(set_names(abundance_responses), ~ log10(to_matrix(veg_grid, paste0("obs_", .x))))

# Sample totals past vs present for each relative response (diagnostic for
# the choice between absolute and relative abundance)
total_abundance_change <- imap_dfr(relative_totals, function(rt, r) {
  samples %>%
    left_join(rt$totals, by = c("plot_id", "period")) %>%
    select(plot_id, period, total) %>%
    pivot_wider(names_from = period, values_from = total, names_prefix = "total_") %>%
    mutate(response = r, variable = resp_var[[r]], log10_ratio = log10(total_present / total_past),
           .before = 1)
})
write_csv(total_abundance_change, check_file("total_abundance_change.csv"))
total_abundance_change %>%
  group_by(response, variable) %>%
  summarise(n_plots            = sum(!is.na(log10_ratio)),
            median_log10_ratio = median(log10_ratio, na.rm = TRUE),
            prop_decline       = mean(log10_ratio < 0, na.rm = TRUE), .groups = "drop")

# Presence without an abundance value (e.g. whole-plot records)
abu_col <- if ("obs_abuall" %in% names(veg_grid)) "obs_abuall" else "obs_abu"
presence_without_abundance <- veg_grid %>%
  filter(pa == 1, is.na(.data[[abu_col]])) %>%
  count(species, period, name = "n_presence_without_abundance")
write_csv(presence_without_abundance, check_file("presence_without_abundance.csv"))
presence_without_abundance

# --- Empty samples and plots -------------------------------------------------
sample_richness <- tibble(sample_id = rownames(Y_pa), richness = rowSums(Y_pa)) %>%
  left_join(samples, by = "sample_id")

empty_samples <- sample_richness %>% filter(richness == 0)
empty_samples

empty_plots <- empty_samples %>%
  count(plot_id, name = "n_empty") %>%
  filter(n_empty == 2 | drop_plots_with_empty_sample) %>%
  pull(plot_id)

write_csv(sample_richness, check_file("sample_richness.csv"))
write_csv(empty_samples, check_file("empty_samples.csv"))
message(nrow(empty_samples), " empty samples; ", length(empty_plots),
        " plots dropped (both samples): ", paste(empty_plots, collapse = ", "))

if (length(empty_plots) > 0) {
  samples <- samples %>% filter(!plot_id %in% empty_plots)
  Y_pa    <- Y_pa[samples$sample_id, , drop = FALSE]
  Y_list  <- map(Y_list, ~ .x[samples$sample_id, , drop = FALSE])
}

# --- Empty species -----------------------------------------------------------
species_prevalence <- tibble(species = colnames(Y_pa)) %>%
  mutate(n_past     = colSums(Y_pa[samples$period == "past", , drop = FALSE]),
         n_present  = colSums(Y_pa[samples$period == "present", , drop = FALSE]),
         n_total    = n_past + n_present,
         n_abu      = colSums(!is.na(Y_list$abu)))
write_csv(species_prevalence, check_file("species_prevalence.csv"))
species_prevalence %>% arrange(n_total)

empty_species <- species_prevalence %>% filter(n_total == 0) %>% pull(species)
one_period_species <- species_prevalence %>%
  filter(n_total > 0, n_past == 0 | n_present == 0) %>% pull(species)

if (length(empty_species) > 0) {
  message("Species with no occurrence in the retained samples (dropped): ",
          paste(empty_species, collapse = ", "))
  Y_pa    <- Y_pa[, !colnames(Y_pa) %in% empty_species, drop = FALSE]
  Y_list  <- map(Y_list, ~ .x[, !colnames(.x) %in% empty_species, drop = FALSE])
}
if (length(one_period_species) > 0) {
  warning("Species occurring in only one period within the retained samples: ",
          paste(one_period_species, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 4. TrData
# -----------------------------------------------------------------------------
# Trait values are taken from the species pool file when it already contains
# all model traits; otherwise from the layer trait table
trait_names <- traits_model[[stratum]]
traits_source <- if (all(trait_names %in% names(species_pool))) species_pool else traits_raw
message("Trait values taken from: ",
        if (identical(traits_source, species_pool)) "species pool file" else "trait table")
traits <- traits_source %>% mutate(species = sp_key(species))

# Climatic affinities joined from the affinity table if not already present
if (!all(c("MAT_mean", "MAP_mean") %in% names(traits))) {
  affinity <- read_csv(resolve_path("affinity", "species_climatic_affinity.csv"),
                       show_col_types = FALSE) %>%
    mutate(species = sp_key(recode(species, !!!affinity_synonyms))) %>%
    select(species, MAT_mean, MAP_mean)
  traits <- traits %>% select(-any_of(c("MAT_mean", "MAP_mean"))) %>%
    left_join(affinity, by = "species")
}

check_cols(traits, c("species", trait_names), "Trait table")

tr_raw <- tibble(species = colnames(Y_pa)) %>%
  left_join(traits %>% select(species, all_of(trait_names)), by = "species")

species_missing_traits <- tr_raw %>%
  filter(if_any(all_of(trait_names), is.na)) %>%
  pivot_longer(all_of(trait_names), names_to = "trait", values_to = "value",
               values_transform = as.character) %>%
  filter(is.na(value)) %>%
  select(species, trait)
write_csv(species_missing_traits, check_file("species_missing_traits.csv"))
species_missing_traits

if (nrow(species_missing_traits) > 0) {
  if (!drop_species_missing_traits) {
    stop(n_distinct(species_missing_traits$species),
         " species lack trait values (see the check_species_missing_traits file in outputs/).")
  }
  drop_sp <- unique(species_missing_traits$species)
  message("Species dropped for missing traits: ", paste(drop_sp, collapse = ", "))
  Y_pa    <- Y_pa[, !colnames(Y_pa) %in% drop_sp, drop = FALSE]
  Y_list  <- map(Y_list, ~ .x[, !colnames(.x) %in% drop_sp, drop = FALSE])
  tr_raw  <- tr_raw %>% filter(!species %in% drop_sp)
}

# Transformations: log10 where required, then explicit scaling of numeric
# traits; categorical traits as factors
TrData <- tr_raw %>%
  mutate(across(any_of(log_traits), log10),
         across(any_of(categorical_traits), as.factor),
         across(where(is.numeric), ~ as.numeric(scale(.x)))) %>%
  column_to_rownames("species") %>%
  as.data.frame()

TrFormula <- as.formula(paste("~", paste(trait_names, collapse = " + ")))

stopifnot("TrData rows do not match Y columns" =
            identical(rownames(TrData), colnames(Y_pa)))

# -----------------------------------------------------------------------------
# 5. Phylogeny
# -----------------------------------------------------------------------------
tree <- tree_raw
tree$tip.label <- sp_key(tree$tip.label)

species_not_in_tree <- setdiff(colnames(Y_pa), tree$tip.label)
if (length(species_not_in_tree) > 0) {
  stop("Species in Y missing from the phylogeny: ",
       paste(species_not_in_tree, collapse = ", "))
}

tree <- keep.tip(tree, colnames(Y_pa))

stopifnot("Phylogeny tips do not match Y columns" =
            setequal(tree$tip.label, colnames(Y_pa)))
if (is.null(tree$edge.length)) stop("The phylogeny has no branch lengths.")
if (!is.ultrametric(tree)) warning("The pruned phylogeny is not ultrametric.")

# Phylogenetic correlation matrix under Brownian motion, from the tree pruned
# to the modelled species. Hmsc computes the same matrix internally from
# `phyloTree` (vcv(..., model = "Brownian", corr = TRUE)); it is stored here
# for inspection. The C matrix saved with the combined phylogeny is not used
# because it covers all project species and is rooted at the full-tree root.
C_phylo <- vcv(tree, model = "Brownian", corr = TRUE)[colnames(Y_pa), colnames(Y_pa)]

# -----------------------------------------------------------------------------
# 6. XData
# -----------------------------------------------------------------------------
# Climate-change predictors
clim_wide <- climate %>%
  filter(plot_id %in% samples$plot_id) %>%
  select(plot_id, period, all_of(clim_read_vars))

clim_past    <- clim_wide %>% filter(period == "past") %>% arrange(plot_id)
clim_present <- clim_wide %>% filter(period == "present") %>% arrange(plot_id)
stopifnot(identical(clim_past$plot_id, clim_present$plot_id))

if (clim_change_source == "pca") {
  # PCA fitted on historical climate; contemporary climate projected onto the
  # same axes, so that both periods share one coordinate system
  pca_clim <- prcomp(clim_past[, clim_vars], center = TRUE, scale. = TRUE)
  pc1_past    <- predict(pca_clim, clim_past[, clim_vars])[, 1]
  pc1_present <- predict(pca_clim, clim_present[, clim_vars])[, 1]

  # Orient PC1 so that positive values correspond to warmer conditions
  if (cor(pc1_past, clim_past[[clim_temp_var]]) < 0) {
    pc1_past <- -pc1_past; pc1_present <- -pc1_present
    pca_clim$rotation[, 1] <- -pca_clim$rotation[, 1]
  }

  clim_pca_summary <- list(
    loadings      = pca_clim$rotation,
    var_explained = round(100 * pca_clim$sdev^2 / sum(pca_clim$sdev^2), 1)
  )
  clim_pca_summary

  plot_delta <- tibble(plot_id = clim_past$plot_id,
                       delta_PC1 = pc1_present - pc1_past)
} else {
  plot_delta <- tibble(plot_id = clim_past$plot_id) %>%
    bind_cols(map_dfc(set_names(clim_change_vars, paste0("delta_", clim_change_vars)),
                      ~ clim_present[[.x]] - clim_past[[.x]]))
}

delta_cols <- setdiff(names(plot_delta), "plot_id")
clim_terms <- paste0("d_", str_remove(delta_cols, "^delta_"), "_sc")

# Elevation, centred and scaled at the plot level
plot_covariates <- plot_meta %>%
  filter(plot_id %in% samples$plot_id) %>%
  left_join(plot_delta, by = "plot_id") %>%
  mutate(elev_sc = as.numeric(scale(elevation)))

# Climate-change predictors: 0 for historical samples, delta for contemporary
# samples; divided by the SD of delta across plots (not centred, so that 0
# keeps its meaning of no change since the historical survey)
clim_delta_sd <- map_dbl(set_names(delta_cols), ~ sd(plot_covariates[[.x]]))

XData_full <- samples %>%
  left_join(plot_covariates, by = "plot_id") %>%
  mutate(time = period)
for (i in seq_along(delta_cols)) {
  XData_full[[clim_terms[i]]] <- if_else(XData_full$period == "past", 0,
                                         XData_full[[delta_cols[i]]]) / clim_delta_sd[[i]]
}

stopifnot("XData rows do not match Y rows" =
            identical(XData_full$sample_id, rownames(Y_pa)))

XData <- XData_full %>%
  select(time, elev_sc, all_of(clim_terms)) %>%
  as.data.frame()
rownames(XData) <- XData_full$sample_id

# --- Model formulas and collinearity diagnostics -----------------------------
clim_sum  <- paste(clim_terms, collapse = " + ")
XFormulas <- list(
  m0_time            = ~ time,
  m1_time_elev       = ~ time * elev_sc,
  m2_climchange      = as.formula(paste("~ elev_sc +", clim_sum)),
  m3_climchange_elev = as.formula(paste("~ (", clim_sum, ") * elev_sc"))
)
XFormulas

collinearity <- imap_dfr(XFormulas, function(f, name) {
  mm <- model.matrix(f, data = XData)[, -1, drop = FALSE]
  vif <- if (ncol(mm) > 1) diag(solve(cor(mm))) else setNames(1, colnames(mm))
  tibble(model = name, term = names(vif), vif = round(vif, 2))
})
collinearity

# Spatial variation in each climate-change predictor and its relation to
# elevation and to time (near-constant change would make it equivalent to time)
clim_change_diagnostics <- map_dfr(seq_along(delta_cols), function(i) {
  d <- plot_covariates[[delta_cols[i]]]
  tibble(predictor            = clim_terms[i],
         delta_mean           = mean(d),
         delta_sd             = sd(d),
         delta_cv             = sd(d) / abs(mean(d)),
         cor_delta_elevation  = cor(d, plot_covariates$elevation),
         cor_with_time        = cor(as.numeric(XData$time == "present"),
                                    XData[[clim_terms[i]]]))
})
clim_change_diagnostics

# Correlation among climate-change predictors (plot level)
if (length(delta_cols) > 1) {
  cor(plot_covariates[, delta_cols])
}

write_csv(collinearity, check_file("collinearity_vif.csv"))
write_csv(clim_change_diagnostics, check_file("clim_change_diagnostics.csv"))
write_csv(plot_covariates, check_file("plot_covariates.csv"))

# -----------------------------------------------------------------------------
# 7. studyDesign and random levels
# -----------------------------------------------------------------------------
studyDesign <- data.frame(sample = factor(XData_full$sample_id),
                          plot   = factor(XData_full$plot_id))

# Projected plot coordinates (km)
plot_xy <- plot_covariates %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(projected_crs) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  transmute(x = X / 1000, y = Y / 1000)
rownames(plot_xy) <- plot_covariates$plot_id
sData <- plot_xy[levels(studyDesign$plot), , drop = FALSE]

dist_km <- as.matrix(dist(sData))
diag(dist_km) <- NA
min_dist_km <- min(dist_km, na.rm = TRUE)
duplicated_coords <- which(dist_km == 0, arr.ind = TRUE)
if (nrow(duplicated_coords) > 0) {
  stop("Plots with identical coordinates: ",
       paste(unique(rownames(sData)[duplicated_coords[, 1]]), collapse = ", "))
}
message("Minimum distance between plots: ", round(min_dist_km * 1000, 1), " m")

# Random levels for one random-level structure (see rl_structures)
build_random_levels <- function(structure) {
  s <- rl_structures[[structure]]
  rL_sample <- HmscRandomLevel(units = levels(studyDesign$sample)) %>%
    setPriors(nfMax = s$nf_max)
  rL_plot <- if (!s$spatial) {
    HmscRandomLevel(units = levels(studyDesign$plot))
  } else if (spatial_method == "NNGP") {
    HmscRandomLevel(sData = sData, sMethod = "NNGP", nNeighbours = 10)
  } else {
    HmscRandomLevel(sData = sData, sMethod = "Full")
  }
  rL_plot <- setPriors(rL_plot, nfMax = s$nf_max)
  if (s$sample) list(sample = rL_sample, plot = rL_plot) else list(plot = rL_plot)
}

# -----------------------------------------------------------------------------
# 8. Final consistency checks
# -----------------------------------------------------------------------------
checks <- tibble(
  check = c("Y_pa rows == XData rows",
            "Y_pa rows == studyDesign rows",
            "abundance matrices: dimnames == Y_pa dimnames",
            "abundance matrices: NA wherever Y_pa == 0",
            "abundance matrices: no zero or infinite values",
            "relative abundance in (0, 1]",
            "TrData rows == Y columns",
            "phylogeny tips == Y columns",
            "plots in studyDesign == plots in sData",
            "no empty samples (or only allowed ones)",
            "no species without occurrences",
            "no NA in TrData",
            "no NA in XData",
            "Y_pa contains only 0/1 (absences filled)",
            "numeric traits scaled (mean 0, SD 1)",
            "elev_sc scaled at plot level (mean 0, SD 1)",
            "climate-change predictors 0 in all historical samples"),
  passed = c(identical(rownames(Y_pa), rownames(XData)),
             identical(rownames(Y_pa), as.character(studyDesign$sample)),
             all(map_lgl(Y_list, ~ identical(dimnames(.x), dimnames(Y_pa)))),
             all(map_lgl(Y_list, ~ all(is.na(.x[Y_pa == 0])))),
             all(map_lgl(Y_list, ~ all(is.finite(.x[!is.na(.x)])))),
             all(map_lgl(Y_list[keep(names(Y_list), is_relative)], ~ all(.x[!is.na(.x)] <= 1e-12))),
             identical(rownames(TrData), colnames(Y_pa)),
             setequal(tree$tip.label, colnames(Y_pa)),
             setequal(levels(studyDesign$plot), rownames(sData)),
             if (drop_plots_with_empty_sample) all(rowSums(Y_pa) > 0)
             else !any(tapply(rowSums(Y_pa), studyDesign$plot, sum) == 0),
             all(colSums(Y_pa) > 0),
             !anyNA(TrData),
             !anyNA(XData),
             all(Y_pa %in% c(0L, 1L)),
             all(map_lgl(select(TrData, where(is.numeric)),
                         ~ abs(mean(.x)) < 1e-8 && abs(sd(.x) - 1) < 1e-8)),
             abs(mean(plot_covariates$elev_sc)) < 1e-8 &&
               abs(sd(plot_covariates$elev_sc) - 1) < 1e-8,
             all(XData[XData$time == "past", clim_terms] == 0))
)
checks
write_csv(checks, check_file("consistency_checks.csv"))
stopifnot("At least one consistency check failed (see `checks`)" = all(checks$passed))

input_summary <- tibble(
  stratum      = stratum,
  pool         = pool_label,
  n_plots      = nlevels(studyDesign$plot),
  n_samples    = nrow(Y_pa),
  n_species    = ncol(Y_pa),
  prevalence_mean = mean(colMeans(Y_pa)),
  n_abu_obs    = sum(!is.na(Y_list$abu)),
  responses    = paste(responses_by_stratum[[stratum]], collapse = ", "),
  traits       = paste(trait_names, collapse = ", ")
)
input_summary

# -----------------------------------------------------------------------------
# 9. Unfitted models
# -----------------------------------------------------------------------------
# Abundance models: species with fewer than `min_abu_obs` values excluded,
# per response
response_species <- imap(Y_list, function(Y, r) {
  min_obs <- if (stratum == "understory" && r %in% c("abu", "rel")) min_abu_obs_paired else min_abu_obs
  n_obs <- colSums(!is.na(Y))
  sp <- colnames(Y)[n_obs >= min_obs]
  message(r, ": ", ncol(Y) - length(sp), " species excluded (< ", min_obs, " observations).")
  sp
})
few_abundance_values <- imap_dfr(Y_list, function(Y, r) {
  n_obs <- colSums(!is.na(Y))
  tibble(response = r, species = names(n_obs), n_values = n_obs,
         included = species %in% response_species[[r]],
         few_values = included & n_obs < min_abu_obs)
})
write_csv(few_abundance_values, check_file("abundance_values_per_species.csv"))
few_abundance_values %>% group_by(response) %>%
  summarise(n_species = sum(included), n_flagged_few_values = sum(few_values), .groups = "drop")

build_model <- function(f, response, structure) {
  if (response == "pa") {
    Y <- Y_pa; distr <- "probit"; sp <- colnames(Y_pa)
  } else {
    sp <- response_species[[response]]; Y <- Y_list[[response]][, sp, drop = FALSE]
    distr <- "normal"
  }
  Hmsc(Y           = Y,
       XData       = XData,
       XFormula    = f,
       TrData      = TrData[sp, , drop = FALSE],
       TrFormula   = TrFormula,
       phyloTree   = keep.tip(tree, sp),
       studyDesign = if (rl_structures[[structure]]$sample) studyDesign else studyDesign[, "plot", drop = FALSE],
       ranLevels   = build_random_levels(structure),
       distr       = distr,
       # log10 abundance scaled internally so that default priors are on a
       # unit scale; parameters are returned on the log10 scale
       YScale      = (distr == "normal"))
}

models_unfitted <- expand_grid(model     = intersect(names(XFormulas), models_to_build),
                               structure = intersect(names(rl_structures), structures_to_build),
                               response  = responses_by_stratum[[stratum]]) %>%
  mutate(model_id = if_else(structure == "default",
                            paste(model, response, sep = "_"),
                            paste(model, structure, response, sep = "_")),
         hmsc     = pmap(list(model, response, structure),
                         function(m, r, s) build_model(XFormulas[[m]], r, s)))

models_unfitted %>% select(model_id, model, structure, response)

# -----------------------------------------------------------------------------
# 10. Save
# -----------------------------------------------------------------------------
hmsc_inputs <- list(
  Y_pa = Y_pa, Y_list = Y_list, Y_abu = Y_list$abu, Y_rel = Y_list$rel,
  response_vars = response_vars[[stratum]][abundance_responses],
  relative_var = resp_var[["rel"]], abundance_var = resp_var[["abu"]],
  total_abundance_change = total_abundance_change, species_excluded = excluded_species,
  XData = XData, XData_full = XData_full,
  TrData = TrData, TrData_raw = tr_raw, log_traits = log_traits,
  TrFormula = TrFormula, tree = tree, C_phylo = C_phylo,
  plot_xy_km = sData, min_dist_km = min_dist_km,
  studyDesign = studyDesign, sData = sData, XFormulas = XFormulas,
  plot_covariates = plot_covariates, clim_delta_sd = clim_delta_sd,
  clim_change_source = clim_change_source, clim_terms = clim_terms,
  clim_pca = if (clim_change_source == "pca") clim_pca_summary else NULL,
  input_summary = input_summary, checks = checks
)

saveRDS(hmsc_inputs,
        clean_file("hmsc_inputs.rds"))
saveRDS(models_unfitted %>% select(model_id, model, structure, response, hmsc),
        clean_file("hmsc_models_unfitted.rds"))
