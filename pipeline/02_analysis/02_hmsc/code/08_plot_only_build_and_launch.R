# =============================================================================
# 08_plot_only_build_and_launch.R
# Author: Agostina Torres
#
# Plot-level-only HMSC models (no sample-level random effect), built and
# launched on the computing machine in one run:
#   1. Builds the unfitted models from models already in the unfitted-model
#      files: each new model copies the data, fixed effects, traits, phylogeny
#      and plot-level random effect of an existing model and drops the sample
#      level (identical to structures nos, nf2_nos, nsp_nos, nf2_nsp_nos of
#      01_prepare_hmsc_inputs.R).
#   2. Launches one background fit per model (02_fit_hmsc.R), medium runs.
#   3. Monitoring lines.
#
# Models, per layer (occurrence and relative abundance):
#   m0 spatial     : ~ time;             spatially explicit plot level
#   m1 non-spatial : ~ time * elevation; plot identity
#   m1 spatial     : ~ time * elevation; spatially explicit plot level
# Traits (Gamma) and phylogeny as in all models: overstory MAT, MAP, seed mass,
# wood density; understory MAT, MAP, SLA, seed mass, woodiness. Latent factors:
# nf_max 2 (overstory), 5 (understory).
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/<stratum>_<pool>_hmsc_models_unfitted.rds
# Outputs: the same files with the new models added (existing models unchanged);
#          fits and logs in pipeline/02_analysis/02_hmsc/outputs/ (see 02_fit_hmsc.R)
# =============================================================================

library(here)
library(Hmsc)
library(dplyr)

# Layers to launch now (both are built in section 1)
strata_to_launch <- c("overstory")
run_mode <- "medium"

# -----------------------------------------------------------------------------
# 1. Build the plot-only models
# -----------------------------------------------------------------------------
pool <- "strict"
module_dir <- here("pipeline", "02_analysis", "02_hmsc")

# Base models (with sample level) from which the plot-only models are built,
# per layer: m0 spatial, m1 non-spatial, m1 spatial; occurrence and relative
# abundance
base_models <- list(
  overstory  = c("m0_time_nf2_pa", "m0_time_nf2_rel",
                 "m1_time_elev_nf2_nsp_pa", "m1_time_elev_nf2_nsp_rel",
                 "m1_time_elev_nf2_pa", "m1_time_elev_nf2_rel"),
  understory = c("m0_time_pa", "m0_time_rel",
                 "m1_time_elev_nsp_pa", "m1_time_elev_nsp_rel",
                 "m1_time_elev_pa", "m1_time_elev_rel")
)

drop_sample_level <- function(m) {
  is_probit <- all(m$distr[, "family"] == 2)
  Hmsc(Y           = m$Y,
       XData       = m$XData,
       XFormula    = m$XFormula,
       TrData      = m$TrData,
       TrFormula   = m$TrFormula,
       phyloTree   = m$phyloTree,
       studyDesign = m$studyDesign[, "plot", drop = FALSE],
       ranLevels   = list(plot = m$ranLevels[["plot"]]),
       distr       = if (is_probit) "probit" else "normal",
       YScale      = !is_probit)
}

for (stratum in names(base_models)) {
  f <- file.path(module_dir, "data_clean", paste0(stratum, "_", pool, "_hmsc_models_unfitted.rds"))
  models <- readRDS(f)
  
  missing <- setdiff(base_models[[stratum]], models$model_id)
  if (length(missing) > 0) stop(stratum, ": base models not found: ", paste(missing, collapse = ", "))
  
  new_rows <- bind_rows(lapply(base_models[[stratum]], function(id) {
    i        <- which(models$model_id == id)
    response <- sub(".*_", "", id)
    new_id   <- paste0(sub(paste0("_", response, "$"), "", id), "_nos_", response)
    row      <- models[i, ]
    row$model_id  <- new_id
    row$structure <- paste0(ifelse(row$structure == "default", "", paste0(row$structure, "_")), "nos")
    row$hmsc      <- list(drop_sample_level(models$hmsc[[i]]))
    row
  }))
  
  models <- bind_rows(filter(models, !model_id %in% new_rows$model_id), new_rows)
  saveRDS(models, f)
  cat(stratum, ": added", nrow(new_rows), "models:\n")
  print(data.frame(model_id = new_rows$model_id,
                   random_levels = sapply(new_rows$hmsc, function(m) paste(names(m$ranLevels), collapse = ", ")),
                   spatial_plot = sapply(new_rows$hmsc, function(m) !is.null(m$ranLevels$plot$s))))
}

# -----------------------------------------------------------------------------
# 2. Launch
# -----------------------------------------------------------------------------
new_ids <- function(ids) sub("_(pa|rel)$", "_nos_\\1", ids)
runs <- do.call(rbind, lapply(names(base_models), function(s)
  data.frame(stratum = s, model_id = new_ids(base_models[[s]]), run_mode = run_mode)))
runs <- runs[runs$stratum %in% strata_to_launch, ]
runs

cores_per_model <- 4   # = n_chains in 02_fit_hmsc.R

# -----------------------------------------------------------------------------
# Checks before launch
# -----------------------------------------------------------------------------
fit_script <- file.path(module_dir, "code", "02_fit_hmsc.R")
out_dir    <- file.path(module_dir, "outputs")
rscript    <- file.path(R.home("bin"), "Rscript")

stopifnot("02_fit_hmsc.R not found" = file.exists(fit_script))

# All requested models must exist in the unfitted-model files
for (s in unique(runs$stratum)) {
  f <- file.path(module_dir, "data_clean", paste0(s, "_", pool, "_hmsc_models_unfitted.rds"))
  if (!file.exists(f)) stop("Unfitted model file not found: ", f)
  available <- readRDS(f)$model_id
  missing   <- setdiff(runs$model_id[runs$stratum == s], available)
  if (length(missing) > 0) {
    stop(s, ": models not found in ", basename(f), ": ", paste(missing, collapse = ", "),
         "\nCheck section 1 of this script.")
  }
}

n_cores_needed <- nrow(runs) * cores_per_model
cat("Processes:", nrow(runs), "| cores needed:", n_cores_needed,
    "| cores available:", parallel::detectCores(), "\n")
if (n_cores_needed > parallel::detectCores()) {
  warning("More cores requested than available; fits will compete for cores.")
}

# -----------------------------------------------------------------------------
# Background processes (one per model)
# -----------------------------------------------------------------------------
quote_arg <- function(x) shQuote(x, type = if (.Platform$OS.type == "windows") "cmd" else "sh")

for (i in seq_len(nrow(runs))) {
  s  <- runs$stratum[i]
  id <- runs$model_id[i]
  mode <- runs$run_mode[i]
  process_log <- file.path(out_dir, paste0(s, "_", pool, "_", mode, "_process_", id, ".txt"))
  system2(rscript,
          args   = c(quote_arg(fit_script), mode, s, id),
          stdout = process_log, stderr = process_log,
          wait   = FALSE)
  cat("Launched:", s, id, mode, "\n")
  Sys.sleep(2)   # stagger start-up (reading the unfitted-model file)
}

# -----------------------------------------------------------------------------
# 3. Monitoring (run these lines at any time)
# -----------------------------------------------------------------------------
# Number of R processes running (on Windows); the launched fits appear as
# Rscript.exe plus their chain workers
if (.Platform$OS.type == "windows") {
  system('tasklist /FI "IMAGENAME eq Rscript.exe"')
}

# Models finished so far (from the run logs)
finished <- do.call(rbind, lapply(split(runs, list(runs$stratum, runs$run_mode), drop = TRUE), function(r) {
  f <- file.path(out_dir, paste0(r$stratum[1], "_", pool, "_", r$run_mode[1], "_run_log.csv"))
  if (file.exists(f)) {
    read.csv(f)[, c("model_id", "stratum", "run_mode", "elapsed_min", "finished_at")] |>
      subset(model_id %in% r$model_id)
  }
}))
finished

# Last lines of each model's console log (written by 02_fit_hmsc.R). The
# process log stays locked while the process runs on Windows; it can be read
# once the process has finished (or failed at start-up).
read_tail <- function(f, n = 3) {
  if (!file.exists(f)) return("(not created yet)")
  tryCatch(suppressWarnings(tail(readLines(f, warn = FALSE), n)),
           error = function(e) "(file in use by the running process)")
}
for (i in seq_len(nrow(runs))) {
  prefix <- paste0(runs$stratum[i], "_", pool, "_", runs$run_mode[i], "_")
  cat("\n---", runs$stratum[i], runs$model_id[i], runs$run_mode[i], "---\n")
  cat(read_tail(file.path(out_dir, paste0(prefix, "console_log_", runs$model_id[i], ".txt"))),
      sep = "\n")
  process_tail <- read_tail(file.path(out_dir, paste0(prefix, "process_", runs$model_id[i], ".txt")))
  if (any(grepl("Error", process_tail))) {
    cat("Process log:\n", process_tail, sep = "\n")
  }
}
