# =============================================================================
# 01c_prelaunch_checks.R
# Author: Agostina Torres
#
# Checks run on the local computer after 01_prepare_hmsc_inputs.R (both strata) and
# before copying the unfitted models to the fitting machine. It verifies that
# every model in the launch plan exists and is well posed, and collects the
# diagnostics that support the remaining modelling choices.
#
#   A. Inputs        : pool, excluded species, numbers of plots/samples/species
#   B. Launch plan   : every model in `runs` exists in the unfitted file
#   C. Models        : per model, species and observations; for normal models,
#                      species with zero variance (fatal with YScale = TRUE),
#                      species observed in one survey only (time effect
#                      informed by the prior and by Gamma only) and samples
#                      with no observation
#   D. Traits        : pairwise correlations among model traits
#   E. Phylogeny     : tips == species, ultrametric, C positive definite
#   F. Total abundance change between surveys (absolute vs relative choice),
#      and its relation to elevation
#   G. Overstory stem size (basal area per stem) by survey and species: a
#      general drop in the present indicates a lower minimum stem diameter,
#      which would make stem counts (relative density) non-comparable
#   H. Expected run time of each run from previous medium runs
#   I. MD5 checksums of the unfitted files, to verify the copy on the fitting
#      machine (tools::md5sum() there must return the same values)
#   J. Species flags: species whose estimates need caution in interpretation,
#      from documented data issues (manual table) and automatic criteria (few
#      abundance records, quasi-separation in the occurrence model, patchy
#      distribution in the present-day subplots). All species stay in the
#      models; the flags are used to mark (or optionally omit) species in the
#      interpretation figures (04_hmsc_results.R) and listed in the
#      Supplementary Information.
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/<stratum>_<pool>_hmsc_inputs.rds
#          pipeline/02_analysis/02_hmsc/data_clean/<stratum>_<pool>_hmsc_models_unfitted.rds
#          pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_medium_run_log.csv (optional)
#          overstory_long.csv (section G)
#          <stratum>_<pool>_diag_separation_screen.csv from 01b (section J)
#          composite_reliability_diagnostic.csv (section J, understory)
# Outputs: pipeline/02_analysis/02_hmsc/outputs/<pool>_prelaunch_*.csv
#          pipeline/02_analysis/02_hmsc/outputs/<pool>_species_flags.csv
# =============================================================================

library(here)
library(tidyverse)
library(Hmsc)
library(ape)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
pool <- "strict"

# Launch plan (same rows as in 02b_launch_fits.R, set "plot_only")
runs <- tribble(
  ~stratum,     ~model_id,                      ~run_mode,
  "overstory",  "m0_time_nf2_nos_pa",           "medium",
  "overstory",  "m0_time_nf2_nos_rel",          "medium",
  "overstory",  "m1_time_elev_nf2_nsp_nos_pa",  "medium",
  "overstory",  "m1_time_elev_nf2_nsp_nos_rel", "medium",
  "overstory",  "m1_time_elev_nf2_nos_pa",      "medium",
  "overstory",  "m1_time_elev_nf2_nos_rel",     "medium",
  "understory", "m0_time_nos_pa",               "medium",
  "understory", "m0_time_nos_rel",              "medium",
  "understory", "m1_time_elev_nsp_nos_pa",      "medium",
  "understory", "m1_time_elev_nsp_nos_rel",     "medium",
  "understory", "m1_time_elev_nos_pa",          "medium",
  "understory", "m1_time_elev_nos_rel",         "medium"
)

min_obs_per_period <- 2      # normal models: minimum observations per survey
min_abu_records    <- 20     # species flag: fewer abundance records in total
patchy_threshold   <- 0.4    # species flag: mean subplot occupancy (understory)

# Species flagged from documented data issues (section J). Species names as
# in the model (genus_species).
#   identification : records inconsistent between surveys, suggesting that
#                    the species was not identified consistently
#   non_climatic   : change dominated by a known non-climatic driver
#   stand_dynamics : change dominated by recruitment or self-thinning
#   complex        : records of a congener merged into this species (see
#                    `species_merge` in 01_prepare_hmsc_inputs.R)
manual_flags <- tribble(
  ~stratum,     ~species,             ~flag,            ~affects, ~reason,
  "overstory",  "Abies_grandis",      "identification", "all",
  "Recorded in both surveys in 1 of 14 plots; present-day basal areas of 0.30-0.64% of plot area where absent in the first survey, and reciprocal replacement by A. amabilis in at least one plot",
  "overstory",  "Pinus_monticola",    "non_climatic",   "all",
  "Occupied plots fall from 28 to 10 and basal area per stem drops to 0.23x, consistent with mortality of large trees from white pine blister rust",
  "overstory",  "Alnus_rubra",        "stand_dynamics", "abu, rel, abucnt, relcnt",
  "Stems increase from 78 to 223 while basal area per stem drops to 0.20x, indicating a recruitment pulse of an early-successional species",
  "understory", "Vaccinium_ovalifolium", "complex",        "all",
  "Includes records of V. alaskaense: ovalifolium recorded in the first survey and alaskaense in the second in 22 plots (2 the reverse), indicating inconsistent separation of the two species between surveys",
  "understory", "Viola_sempervirens",    "complex",        "all",
  "Includes records of V. orbiculata: sempervirens recorded in the first survey and orbiculata in the second in 38 plots (0 the reverse), indicating inconsistent separation of the two species between surveys"
)
max_trait_cor      <- 0.7    # |r| above which a trait pair is flagged
iter_factor        <- c(short = 0.1, medium = 1, long = 10)   # thin relative to medium

# -----------------------------------------------------------------------------
# Paths and helpers
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
clean_dir  <- file.path(module_dir, "data_clean")
out_dir    <- file.path(module_dir, "outputs")
out_file   <- function(name) file.path(out_dir, paste0(pool, "_prelaunch_", name))

strata <- unique(runs$stratum)
tag    <- function(s) paste(s, pool, sep = "_")
inputs   <- map(set_names(strata), ~ readRDS(file.path(clean_dir, paste0(tag(.x), "_hmsc_inputs.rds"))))
unfitted <- map(set_names(strata), ~ readRDS(file.path(clean_dir, paste0(tag(.x), "_hmsc_models_unfitted.rds"))))

checklist <- tibble(section = character(), stratum = character(), check = character(),
                    status = character(), detail = character())
add_check <- function(section, stratum, check, ok, detail = "", level = "fail") {
  checklist <<- add_row(checklist, section = section, stratum = stratum, check = check,
                        status = if (isTRUE(ok)) "ok" else level, detail = detail)
}

# -----------------------------------------------------------------------------
# A. Inputs
# -----------------------------------------------------------------------------
input_overview <- imap_dfr(inputs, function(inp, s) {
  inp$input_summary %>%
    mutate(species_excluded = paste(inp$species_excluded, collapse = ", "),
           responses        = paste(names(inp$Y_list %||% list()), collapse = ", "))
})
input_overview

iwalk(inputs, function(inp, s) {
  add_check("A", s, "01 consistency checks passed", all(inp$checks$passed))
  add_check("A", s, "pool label", inp$input_summary$pool == pool, inp$input_summary$pool)
  add_check("A", s, "relative abundance matrix present", !is.null(inp$Y_rel))
})

# -----------------------------------------------------------------------------
# B. Launch plan
# -----------------------------------------------------------------------------
for (s in strata) {
  missing <- setdiff(runs$model_id[runs$stratum == s], unfitted[[s]]$model_id)
  add_check("B", s, "all launch-plan models exist", length(missing) == 0,
            paste(missing, collapse = ", "))
}

# -----------------------------------------------------------------------------
# C. Models
# -----------------------------------------------------------------------------
model_overview <- pmap_dfr(runs, function(stratum, model_id, run_mode) {
  row <- unfitted[[stratum]] %>% filter(model_id == !!model_id)
  if (nrow(row) == 0) return(NULL)
  m      <- row$hmsc[[1]]
  Y      <- m$Y
  period <- inputs[[stratum]]$XData_full$period[match(rownames(Y), inputs[[stratum]]$XData_full$sample_id)]
  obs    <- !is.na(Y) & if (m$distr[1, "family"] == 2) Y == 1 else TRUE
  n_past    <- colSums(obs[period == "past", , drop = FALSE])
  n_present <- colSums(obs[period == "present", , drop = FALSE])
  sd_y      <- apply(Y, 2, sd, na.rm = TRUE)
  is_normal <- m$distr[1, "family"] == 1
  tibble(stratum, model_id, run_mode,
         distr            = if (is_normal) "normal" else "probit",
         fixed            = paste(deparse(m$XFormula), collapse = ""),
         ny = m$ny, ns = m$ns, nc = m$nc, nt = m$nt,
         random_levels    = paste(m$rLNames, collapse = " + "),
         nf_max           = paste(map_dbl(m$rL, "nfMax"), collapse = " / "),
         spatial_plot     = !is.null(m$rL$plot$s),
         n_obs            = sum(!is.na(Y)),
         n_zero_var       = sum(is_normal & (is.na(sd_y) | sd_y == 0)),
         zero_var_species = if (is_normal) paste(names(sd_y)[is.na(sd_y) | sd_y == 0], collapse = ", ") else "",
         n_one_period     = sum(n_past < min_obs_per_period | n_present < min_obs_per_period),
         one_period_species = paste(names(n_past)[n_past < min_obs_per_period |
                                                  n_present < min_obs_per_period], collapse = ", "),
         n_empty_samples  = sum(rowSums(!is.na(Y)) == 0))
})
model_overview %>% select(-ends_with("species")) %>% print(width = Inf)

pwalk(model_overview, function(...) {
  r <- list(...)
  add_check("C", r$stratum, paste(r$model_id, "no zero-variance species"),
            r$n_zero_var == 0, r$zero_var_species)
  add_check("C", r$stratum, paste(r$model_id, "species observed in both surveys (>=",
                                  min_obs_per_period, ")"),
            r$n_one_period == 0, r$one_period_species, level = "warning")
  add_check("C", r$stratum, paste(r$model_id, "plot level non-spatial"),
            !r$spatial_plot, level = "warning")
})

# -----------------------------------------------------------------------------
# D. Traits
# -----------------------------------------------------------------------------
trait_cor <- imap_dfr(inputs, function(inp, s) {
  num <- inp$TrData %>% select(where(is.numeric))
  cm  <- cor(num)
  as_tibble(cm, rownames = "trait1") %>%
    pivot_longer(-trait1, names_to = "trait2", values_to = "r") %>%
    filter(trait1 < trait2) %>%
    mutate(stratum = s)
})
trait_cor %>% arrange(desc(abs(r)))

for (s in strata) {
  high <- trait_cor %>% filter(stratum == s, abs(r) > max_trait_cor)
  add_check("D", s, paste("trait correlations |r| <=", max_trait_cor), nrow(high) == 0,
            paste(high$trait1, high$trait2, round(high$r, 2), collapse = "; "), level = "warning")
}

# -----------------------------------------------------------------------------
# E. Phylogeny
# -----------------------------------------------------------------------------
iwalk(inputs, function(inp, s) {
  add_check("E", s, "tree tips == species", setequal(inp$tree$tip.label, colnames(inp$Y_pa)))
  add_check("E", s, "tree ultrametric", is.ultrametric(inp$tree), level = "warning")
  C <- vcv(keep.tip(inp$tree, colnames(inp$Y_pa)), corr = TRUE)
  min_eig <- min(eigen(C, symmetric = TRUE, only.values = TRUE)$values)
  add_check("E", s, "phylogenetic correlation matrix positive definite", min_eig > 1e-8,
            paste("min eigenvalue =", signif(min_eig, 3)))
})

# -----------------------------------------------------------------------------
# F. Total abundance change between surveys
# -----------------------------------------------------------------------------
# Requires the current version of 01 (total_abundance_change in the inputs)
iwalk(inputs, function(inp, s) {
  add_check("F", s, "inputs built with the current 01 (total abundance table present)",
            !is.null(inp$total_abundance_change), "re-run 01_prepare_hmsc_inputs.R")
})
total_change <- imap_dfr(compact(map(inputs, "total_abundance_change")), function(tac, s) {
  inp <- inputs[[s]]
  tac %>%
    left_join(select(inp$plot_covariates, plot_id, elevation), by = "plot_id") %>%
    mutate(stratum = s, .before = 1)
})
total_change_summary <- total_change %>%
  group_by(stratum, response, variable) %>%
  summarise(n_plots            = sum(!is.na(log10_ratio)),
            median_log10_ratio = median(log10_ratio, na.rm = TRUE),
            median_fold_change = 10^median_log10_ratio,
            prop_decline       = mean(log10_ratio < 0, na.rm = TRUE),
            cor_elevation      = cor(log10_ratio, elevation, use = "complete.obs"),
            .groups = "drop")
total_change_summary

# -----------------------------------------------------------------------------
# G. Overstory stem size by survey
# -----------------------------------------------------------------------------
veg_file <- list.files(here("pipeline"), pattern = "^overstory_long\\.csv$",
                       recursive = TRUE, full.names = TRUE)
veg_file <- veg_file[!str_detect(veg_file, "/tmp_files/")]

if (length(veg_file) == 1 && "overstory" %in% strata) {
  veg_os <- read_csv(veg_file, show_col_types = FALSE) %>%
    filter(count > 0, plot_id %in% inputs$overstory$plot_covariates$plot_id)

  stem_size_period <- veg_os %>%
    group_by(period) %>%
    summarise(stems_total      = sum(count),
              stems_per_sample = sum(count) / n_distinct(plot_id),
              ba_per_stem_m2   = sum(basal_area_spp_m2) / sum(count),
              years            = paste(range(year), collapse = "-"))
  stem_size_period

  stem_size_species <- veg_os %>%
    group_by(species, period) %>%
    summarise(ba_per_stem = sum(basal_area_spp_m2) / sum(count), stems = sum(count),
              .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = c(ba_per_stem, stems)) %>%
    mutate(ratio_present_past = ba_per_stem_present / ba_per_stem_past) %>%
    arrange(ratio_present_past)
  stem_size_species %>% print(n = Inf)
  write_csv(stem_size_species, out_file("overstory_stem_size_species.csv"))

  prop_smaller <- mean(stem_size_species$ratio_present_past < 1, na.rm = TRUE)
  add_check("G", "overstory", "basal area per stem not lower in most species",
            prop_smaller < 0.75,
            paste0(round(100 * prop_smaller), "% of species with smaller stems in the present"),
            level = "warning")
} else {
  add_check("G", "overstory", "overstory_long.csv found (one copy)", FALSE,
            paste(veg_file, collapse = "; "), level = "warning")
}

# -----------------------------------------------------------------------------
# H. Expected run time
# -----------------------------------------------------------------------------
medium_logs <- map_dfr(strata, function(s) {
  f <- file.path(out_dir, paste0(tag(s), "_medium_run_log.csv"))
  if (!file.exists(f)) return(NULL)
  read_csv(f, show_col_types = FALSE) %>%
    group_by(model_id) %>% slice_tail(n = 1) %>% ungroup() %>%
    transmute(stratum = s, model_id, medium_min = elapsed_min)
})

# Runs without a medium fit of their own (e.g. new abundance responses) take
# the time of the matching abu model of the same structure
expected_time <- runs %>%
  mutate(proxy_id = str_replace(model_id, "_(rel|abucnt|relcnt)$", "_abu")) %>%
  left_join(medium_logs, by = c("stratum", "proxy_id" = "model_id")) %>%
  mutate(expected_h = round(medium_min * iter_factor[run_mode] / 60, 1))
expected_time
message("Longest expected run: ", max(expected_time$expected_h, na.rm = TRUE), " h")

# -----------------------------------------------------------------------------
# I. Checksums of the unfitted files
# -----------------------------------------------------------------------------
unfitted_files <- file.path(clean_dir, paste0(tag(strata), "_hmsc_models_unfitted.rds"))
checksums <- tibble(file = basename(unfitted_files), md5 = unname(tools::md5sum(unfitted_files)),
                    modified = file.mtime(unfitted_files))
checksums

# -----------------------------------------------------------------------------
# J. Species flags
# -----------------------------------------------------------------------------
model_species <- imap_dfr(inputs, ~ tibble(stratum = .y, species = colnames(.x$Y_pa)))

# Few abundance records (abundance and relative-abundance models)
flags_records <- imap_dfr(inputs, function(inp, s) {
  tibble(stratum = s, species = colnames(inp$Y_abu),
         n_abu = colSums(!is.na(inp$Y_abu))) %>%
    filter(n_abu < min_abu_records) %>%
    transmute(stratum, species, flag = "few_abundance_records", affects = "abu, rel, abucnt, relcnt",
              reason = paste0(n_abu, " abundance records in total (< ", min_abu_records, ")"))
})

# Quasi-separation in the occurrence model (01b screen)
flags_separation <- map_dfr(strata, function(s) {
  f <- file.path(out_dir, paste0(tag(s), "_diag_separation_screen.csv"))
  if (!file.exists(f)) {
    add_check("J", s, "01b separation screen found", FALSE, basename(f), level = "warning")
    return(NULL)
  }
  read_csv(f, show_col_types = FALSE) %>%
    filter(flag_cells) %>%
    transmute(stratum = s, species, flag = "quasi_separation", affects = "pa",
              reason = paste0("No occurrence in period x elevation-tercile cells: ", cells_zero))
})

# Patchy distribution in the present-day subplots (understory)
rel_file <- list.files(here("pipeline"), pattern = "^composite_reliability_diagnostic\\.csv$",
                       recursive = TRUE, full.names = TRUE)
rel_file <- rel_file[!str_detect(rel_file, "/tmp_files/")]
flags_patchy <- if ("understory" %in% strata && length(rel_file) == 1) {
  read_csv(rel_file, show_col_types = FALSE) %>%
    transmute(species = str_replace_all(str_squish(species), " ", "_"),
              occupancy = mean_subplot_occupancy) %>%
    distinct() %>%
    filter(species %in% colnames(inputs$understory$Y_pa), occupancy < patchy_threshold) %>%
    transmute(stratum = "understory", species, flag = "patchy", affects = "abu, rel",
              reason = paste0("Mean subplot occupancy ", round(occupancy, 2), " (< ",
                              patchy_threshold, "); present-day cover poorly characterised"))
} else {
  add_check("J", "understory", "composite_reliability_diagnostic.csv found (one copy)", FALSE,
            paste(rel_file, collapse = "; "), level = "warning")
  NULL
}

species_flags <- bind_rows(
  manual_flags %>% mutate(source = "manual"),
  flags_records %>% mutate(source = "automatic"),
  flags_separation %>% mutate(source = "automatic"),
  flags_patchy %>% mutate(source = "automatic")
) %>%
  semi_join(model_species, by = c("stratum", "species")) %>%
  arrange(stratum, species, flag)
species_flags %>% print(n = Inf, width = Inf)

# Manual flags for species not in the models (name mismatch or not in the pool)
manual_not_found <- anti_join(manual_flags, model_species, by = c("stratum", "species"))
add_check("J", "both", "manual flags match model species", nrow(manual_not_found) == 0,
          paste(manual_not_found$species, collapse = ", "), level = "warning")

# One row per species, for the Supplementary Information
species_flags_summary <- species_flags %>%
  group_by(stratum, species) %>%
  summarise(flags   = paste(unique(flag), collapse = "; "),
            affects = paste(unique(affects), collapse = "; "),
            reasons = paste(reason, collapse = " | "),
            .groups = "drop")
species_flags_summary %>% count(stratum)

write_csv(species_flags,         file.path(out_dir, paste0(pool, "_species_flags.csv")))
write_csv(species_flags_summary, out_file("species_flags_summary.csv"))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
checklist %>% count(status)
checklist %>% filter(status != "ok") %>% print(n = Inf, width = Inf)

write_csv(checklist,            out_file("checklist.csv"))
write_csv(input_overview,       out_file("inputs.csv"))
write_csv(model_overview,       out_file("models.csv"))
write_csv(trait_cor,            out_file("trait_correlations.csv"))
write_csv(total_change_summary, out_file("total_abundance_change.csv"))
write_csv(expected_time,        out_file("expected_time.csv"))
write_csv(checksums,            out_file("checksums.csv"))

if (any(checklist$status == "fail")) {
  stop("At least one check failed (status 'fail'); do not launch before resolving it.")
}
