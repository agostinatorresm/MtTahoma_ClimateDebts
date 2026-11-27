# =============================================================================
# 02b_launch_fits.R
# Author: Agostina Torres
#
# Launches several HMSC fits in parallel from an R session. Each model runs in
# its own background R process (Rscript 02_fit_hmsc.R <run_mode> <stratum>
# <model_id>), so the R session stays free and the fits continue if RStudio is
# closed. Within each process the chains run in parallel (4 cores per model).
#
# Each process writes the usual outputs of 02_fit_hmsc.R (fitted model, run
# log, console log) plus a process log with anything printed before the
# console log is opened (e.g. start-up errors).
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/<stratum>_<pool>_hmsc_models_unfitted.rds
# Outputs: pipeline/02_analysis/02_hmsc/outputs/ (see 02_fit_hmsc.R)
#          pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_process_<model_id>.txt
# =============================================================================

library(here)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
pool <- "strict"

# Models to fit: one background process per row, each with its own run mode
# (same rows as in 01c_prelaunch_checks.R and 02c_monitor_fits.R). All models:
# non-spatial plot level, sample + plot random levels (nf_max 2 overstory, 5
# understory). All abundance responses are fitted in parallel:
#   - time x elevation (m1), long: occurrence; relative abundance (relative
#     basal area, relative cover; main abundance responses); absolute
#     abundance (basal area % of plot area, cover); overstory stem counts
#     (relative density, stem density)
#   - time only (m0, robustness), medium: occurrence and relative abundance
runs <- data.frame(
  stratum  = c(rep("overstory", 7), rep("understory", 5)),
  model_id = c("m1_time_elev_nf2_nsp_pa", "m1_time_elev_nf2_nsp_rel", "m1_time_elev_nf2_nsp_abu",
               "m1_time_elev_nf2_nsp_relcnt", "m1_time_elev_nf2_nsp_abucnt",
               "m0_time_nf2_nsp_pa", "m0_time_nf2_nsp_rel",
               "m1_time_elev_nsp_pa", "m1_time_elev_nsp_rel", "m1_time_elev_nsp_abu",
               "m0_time_nsp_pa", "m0_time_nsp_rel"),
  run_mode = c("long", "long", "long", "long", "long", "medium", "medium",
               "long", "long", "long", "medium", "medium")
)

# Layers to launch in this session (e.g. only one layer when the other is
# already running and its inputs have not changed)
strata_to_launch <- c("overstory", "understory")
runs <- runs[runs$stratum %in% strata_to_launch, ]
runs

cores_per_model <- 4   # = n_chains in 02_fit_hmsc.R

# -----------------------------------------------------------------------------
# 1. Checks
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
fit_script <- file.path(module_dir, "code", "02_fit_hmsc.R")
out_dir    <- file.path(module_dir, "outputs")
rscript    <- file.path(R.home("bin"), "Rscript")

stopifnot("02_fit_hmsc.R not found" = file.exists(fit_script))

# All requested models must exist in the unfitted-model files
for (s in unique(runs$stratum)) {
  f <- file.path(module_dir, "data_clean", paste0(s, "_", pool, "_hmsc_models_unfitted.rds"))
  if (!file.exists(f)) stop("Unfitted model file not found: ", f)
  # Must match the checksum written by 01c_prelaunch_checks.R on the local computer
  # (<pool>_prelaunch_checksums.csv), confirming that the copy is complete
  message(basename(f), "  md5: ", unname(tools::md5sum(f)))
  available <- readRDS(f)$model_id
  missing   <- setdiff(runs$model_id[runs$stratum == s], available)
  if (length(missing) > 0) {
    stop(s, ": models not found in ", basename(f), ": ", paste(missing, collapse = ", "),
         "\nRe-run 01_prepare_hmsc_inputs.R and copy the new file.")
  }
}

n_cores_needed <- nrow(runs) * cores_per_model
cat("Processes:", nrow(runs), "| cores needed:", n_cores_needed,
    "| cores available:", parallel::detectCores(), "\n")
if (n_cores_needed > parallel::detectCores()) {
  warning("More cores requested than available; fits will compete for cores.")
}

# -----------------------------------------------------------------------------
# 2. Launch (one background process per model)
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
