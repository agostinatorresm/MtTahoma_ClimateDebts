# =============================================================================
# 11_launch_oct09.R
# Author: Agostina Torres
#
# Understory HMSC runs of 2026-10-09, built and launched on the computing
# machine in one run. All fits are written to outputs/run_oct09/.
#
#   1. Builds, from the unfitted-model file written by 01_prepare_hmsc_inputs.R
#      (plot-level random effect only; abundance paired by plot):
#        - occurrence models on past data reduced to the present sampling
#          effort (04_robustness 05, cover-weighted removal): n_reduced
#          datasets drawn at random from the 100; present data unchanged
#        - occurrence model without sparse species (04_robustness 03)
#        - relative-cover model (paired) without patchy species
#          (04_robustness 01b)
#   2. Writes one lane script per group of 4 cores: each lane fits its models
#      one after the other (02_fit_hmsc.R), so all runs are launched at once
#      and continue unattended.
#   3. Launches the lanes and prints monitoring lines.
#
# Main structure: m1 non-spatial, ~ time * elevation, non-spatial plot level.
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/understory_strict_hmsc_models_unfitted.rds
#          pipeline/04_robustness/outputs/rarefaction/rarefied_past_occurrence.rds
#          pipeline/04_robustness/outputs/occurrence/area_sensitivity_species.csv
#          pipeline/04_robustness/outputs/accumulation_curve/accumulation_curve_summary.csv
# Outputs: models added to the unfitted-model file;
#          pipeline/02_analysis/02_hmsc/data_clean/understory_reduced_datasets_oct09.csv;
#          pipeline/02_analysis/02_hmsc/outputs/run_oct09/ (fits, run logs, lane scripts and logs)
# =============================================================================

library(here)
library(Hmsc)
library(dplyr)
library(ape)

stratum    <- "understory"
pool       <- "strict"
out_subdir <- "run_oct09"
n_reduced  <- 10
launch     <- TRUE    # FALSE builds the models and lane scripts without launching
set.seed(20261009)

module_dir    <- here("pipeline", "02_analysis", "02_hmsc")
unfitted_file <- file.path(module_dir, "data_clean", paste0(stratum, "_", pool, "_hmsc_models_unfitted.rds"))
reduced_file  <- here("pipeline", "04_robustness", "outputs", "rarefaction", "rarefied_past_occurrence.rds")
area_file     <- here("pipeline", "04_robustness", "outputs", "occurrence", "area_sensitivity_species.csv")
patchy_file   <- here("pipeline", "04_robustness", "outputs", "accumulation_curve", "accumulation_curve_summary.csv")
out_dir       <- file.path(module_dir, "outputs", out_subdir)
lane_dir      <- file.path(out_dir, "lanes")
dir.create(lane_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(unfitted_file, reduced_file, area_file, patchy_file)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}

base_pa  <- "m1_time_elev_nsp_nos_pa"
base_rel <- "m1_time_elev_nsp_nos_rel"
reduced_ids <- sprintf("m1_time_elev_nsp_nos_red%02d_pa", seq_len(n_reduced))
area_id     <- "m1_time_elev_nsp_nos_noarea_pa"
patchy_id   <- "m1_time_elev_nsp_nos_nonpatchy_rel"

sp_key <- function(x) gsub("[ _]+", "_", trimws(x))
models <- readRDS(unfitted_file)
needed <- c(base_pa, base_rel, "m1_time_elev_nsp_nos_abu", "m1_time_elev_nsp_nos_relall",
            "m1_time_elev_nos_pa", "m1_time_elev_nos_rel", "m0_time_nos_pa", "m0_time_nos_rel")
missing <- setdiff(needed, models$model_id)
if (length(missing) > 0) stop("Models not found in the unfitted-model file (rerun 01): ", paste(missing, collapse = ", "))
get_model <- function(id) models$hmsc[[which(models$model_id == id)]]

rebuild <- function(m, Y = m$Y, keep_sp = colnames(m$Y)) {
  is_probit <- all(m$distr[, "family"] == 2)
  tree <- keep.tip(m$phyloTree, m$phyloTree$tip.label[m$phyloTree$tip.label %in% keep_sp])
  Hmsc(Y = Y[, keep_sp, drop = FALSE], XData = m$XData, XFormula = m$XFormula,
       TrData = m$TrData[keep_sp, , drop = FALSE], TrFormula = m$TrFormula, phyloTree = tree,
       studyDesign = m$studyDesign, ranLevels = m$ranLevels,
       distr = if (is_probit) "probit" else "normal", YScale = !is_probit)
}
new_row <- function(base_id, id, structure_suffix, hmsc) {
  row <- models[which(models$model_id == base_id), ]
  row$model_id  <- id
  row$structure <- paste0(row$structure, "_", structure_suffix)
  row$hmsc      <- list(hmsc)
  row
}

# -----------------------------------------------------------------------------
# 1a. Occurrence on reduced past data
# -----------------------------------------------------------------------------
m_pa <- get_model(base_pa)
time <- as.character(m_pa$XData$time); plot <- as.character(m_pa$studyDesign$plot)
past_rows <- which(time == "past")
species <- colnames(m_pa$Y)

reduced <- readRDS(reduced_file)
if (!"weighting" %in% names(reduced)) reduced$weighting <- "cover"
reduced <- reduced[reduced$weighting == "cover", ]
reps <- sort(sample(unique(reduced$rep), n_reduced))
write.csv(data.frame(model_id = reduced_ids, reduced_rep = reps),
          file.path(module_dir, "data_clean", "understory_reduced_datasets_oct09.csv"), row.names = FALSE)

reduced_rows <- lapply(seq_along(reps), function(k) {
  rr  <- reduced[reduced$rep == reps[k], ]
  key <- paste(rr$plot_id, sp_key(rr$species))
  Y <- m_pa$Y
  for (j in seq_along(species)) {
    Y[past_rows, j] <- as.integer(paste(plot[past_rows], sp_key(species[j])) %in% key)
  }
  new_row(base_pa, reduced_ids[k], "reduced", rebuild(m_pa, Y))
})
n_past_obs <- sum(m_pa$Y[past_rows, ])
n_past_red <- sapply(reduced_rows, function(r) sum(r$hmsc[[1]]$Y[past_rows, ]))
cat("Past occurrences: observed", n_past_obs, "| reduced (mean of", n_reduced, "):", round(mean(n_past_red)), "\n")
stopifnot("Reduced past occurrences must not exceed the observed ones" = all(n_past_red <= n_past_obs),
          "Present rows must be unchanged" =
            all(sapply(reduced_rows, function(r) identical(r$hmsc[[1]]$Y[-past_rows, ], m_pa$Y[-past_rows, ]))))

# -----------------------------------------------------------------------------
# 1b. Occurrence without sparse species
# -----------------------------------------------------------------------------
area <- read.csv(area_file)
area_sp <- sp_key(area$species[area$area_sensitive %in% c(TRUE, "TRUE")])
keep_area <- species[!sp_key(species) %in% area_sp]
cat("Sparse species removed:", length(species) - length(keep_area), "->",
    paste(setdiff(species, keep_area), collapse = ", "), "\n")
area_row <- new_row(base_pa, area_id, "noarea", rebuild(m_pa, keep_sp = keep_area))

# -----------------------------------------------------------------------------
# 1c. Paired relative cover without patchy species
# -----------------------------------------------------------------------------
m_rel <- get_model(base_rel)
patchy <- read.csv(patchy_file)
patchy_sp <- sp_key(patchy$species[patchy$patchy %in% c(TRUE, "TRUE")])
keep_patchy <- colnames(m_rel$Y)[!sp_key(colnames(m_rel$Y)) %in% patchy_sp]
cat("Patchy species removed from relative cover:", ncol(m_rel$Y) - length(keep_patchy), "\n")
patchy_row <- new_row(base_rel, patchy_id, "nonpatchy", rebuild(m_rel, keep_sp = keep_patchy))

new_rows <- bind_rows(c(reduced_rows, list(area_row, patchy_row)))
models <- bind_rows(filter(models, !model_id %in% new_rows$model_id), new_rows)
saveRDS(models, unfitted_file)
cat("Models added to", basename(unfitted_file), ":", nrow(new_rows), "\n")

# -----------------------------------------------------------------------------
# 2. Lanes (4 cores each; models within a lane run one after the other)
# -----------------------------------------------------------------------------
# Each step: run mode and model ids
lanes <- list(
  lane01 = list(c("long", base_pa), c("long", "m1_time_elev_nsp_nos_relall")),
  lane02 = list(c("long", base_rel), c("long", "m1_time_elev_nsp_nos_abu"))
)
checks_medium <- c(area_id, patchy_id, "m1_time_elev_nos_pa", "m1_time_elev_nos_rel",
                   "m0_time_nos_pa", "m0_time_nos_rel", base_pa, base_rel)
for (k in seq_len(n_reduced)) {
  steps <- list(c("long", reduced_ids[k]))
  if (k <= length(checks_medium)) steps <- c(steps, list(c("medium", checks_medium[k])))
  lanes[[sprintf("lane%02d", k + 2)]] <- steps
}

plan <- bind_rows(lapply(names(lanes), function(l)
  data.frame(lane = l, step = seq_along(lanes[[l]]),
             run_mode = sapply(lanes[[l]], `[`, 1), model_id = sapply(lanes[[l]], `[`, 2))))
write.csv(plan, file.path(out_dir, "launch_plan_oct09.csv"), row.names = FALSE)
cat("\nLaunch plan (", length(lanes), "lanes x 4 cores =", 4 * length(lanes), "cores; available:",
    parallel::detectCores(), ")\n")
print(plan, row.names = FALSE)

fit_script <- file.path(module_dir, "code", "02_fit_hmsc.R")
stopifnot("02_fit_hmsc.R not found" = file.exists(fit_script))
rscript <- file.path(R.home("bin"), "Rscript")
q <- function(x) shQuote(x, type = if (.Platform$OS.type == "windows") "cmd" else "sh")

lane_files <- sapply(names(lanes), function(l) {
  lines <- c("# Lane script written by 11_launch_oct09.R: fits its models one after the other",
             sprintf("rscript <- %s", deparse(rscript)),
             sprintf("fit_script <- %s", deparse(fit_script)))
  for (st in lanes[[l]]) {
    log_file <- file.path(lane_dir, paste0(l, "_", st[1], "_", st[2], ".txt"))
    lines <- c(lines,
               sprintf('cat("Start %s %s:", format(Sys.time()), "\\n")', st[1], st[2]),
               sprintf('system2(rscript, c(shQuote(fit_script), "%s", "%s", "%s", "%s"), stdout = %s, stderr = %s, wait = TRUE)',
                       st[1], stratum, st[2], out_subdir, deparse(log_file), deparse(log_file)),
               sprintf('cat("End %s %s:", format(Sys.time()), "\\n")', st[1], st[2]))
  }
  f <- file.path(lane_dir, paste0(l, ".R"))
  writeLines(lines, f)
  f
})

# -----------------------------------------------------------------------------
# 3. Launch
# -----------------------------------------------------------------------------
if (launch) {
  for (l in names(lane_files)) {
    system2(rscript, q(lane_files[[l]]), stdout = file.path(lane_dir, paste0(l, "_lane_log.txt")),
            stderr = file.path(lane_dir, paste0(l, "_lane_log.txt")), wait = FALSE)
    cat("Launched", l, "\n")
    Sys.sleep(5)   # stagger start-up (reading the unfitted-model file)
  }
}

# -----------------------------------------------------------------------------
# 4. Monitoring (run these lines at any time)
# -----------------------------------------------------------------------------
if (.Platform$OS.type == "windows") system('tasklist /FI "IMAGENAME eq Rscript.exe"')
for (rm_ in c("long", "medium")) {
  f <- file.path(out_dir, paste0(stratum, "_", pool, "_", rm_, "_run_log.csv"))
  if (file.exists(f)) print(read.csv(f)[, c("model_id", "run_mode", "elapsed_min", "finished_at")])
}
for (f in list.files(lane_dir, pattern = "_lane_log\\.txt$", full.names = TRUE)) {
  cat("\n---", basename(f), "---\n")
  cat(tryCatch(readLines(f, warn = FALSE), error = function(e) "(in use)"), sep = "\n")
}
