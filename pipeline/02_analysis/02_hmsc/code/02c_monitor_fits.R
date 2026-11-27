# =============================================================================
# 02c_monitor_fits.R
# Author: Agostina Torres
#
# Status of the HMSC fits launched with 02b_launch_fits.R. Safe to run at any
# time and from a new R session: it only reads files and lists processes; it
# never launches or stops fits.
#
#   1. Running R processes (each fit appears as Rscript.exe plus one worker
#      per chain)
#   2. Fits finished so far, from the run logs written by 02_fit_hmsc.R
#   3. Last lines of each fit's console log; start-up errors from the process
#      log (readable only once the process has stopped)
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_run_log.csv
#          pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_console_log_<model_id>.txt
#          pipeline/02_analysis/02_hmsc/outputs/<stratum>_<pool>_<run_mode>_process_<model_id>.txt
# =============================================================================

library(here)

# -----------------------------------------------------------------------------
# 0. Settings (same runs as in 02b_launch_fits.R)
# -----------------------------------------------------------------------------
pool <- "strict"

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

out_dir <- here("pipeline", "02_analysis", "02_hmsc", "outputs")

# -----------------------------------------------------------------------------
# 1. Running processes
# -----------------------------------------------------------------------------
if (.Platform$OS.type == "windows") {
  procs <- system('tasklist /FI "IMAGENAME eq Rscript.exe" /FO CSV /NH', intern = TRUE)
  procs <- procs[grepl("Rscript", procs)]
  cat("Rscript processes running:", length(procs),
      "(each fit = 1 main process + 1 per chain)\n")
}

# -----------------------------------------------------------------------------
# 2. Finished fits
# -----------------------------------------------------------------------------
finished <- do.call(rbind, lapply(split(runs, list(runs$stratum, runs$run_mode), drop = TRUE), function(r) {
  f <- file.path(out_dir, paste0(r$stratum[1], "_", pool, "_", r$run_mode[1], "_run_log.csv"))
  if (file.exists(f)) {
    read.csv(f)[, c("model_id", "stratum", "run_mode", "elapsed_min", "finished_at")] |>
      subset(model_id %in% r$model_id)
  }
}))
if (is.null(finished)) {
  finished <- data.frame(model_id = character(), stratum = character(), run_mode = character(),
                         elapsed_min = numeric(), finished_at = character())
}
status <- merge(runs, finished, by = c("stratum", "model_id", "run_mode"), all.x = TRUE)
status$state <- ifelse(is.na(status$finished_at), "running", "finished")
status[order(status$stratum, status$run_mode), ]

# -----------------------------------------------------------------------------
# 3. Console logs and start-up errors
# -----------------------------------------------------------------------------
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
  process_tail <- read_tail(file.path(out_dir, paste0(prefix, "process_", runs$model_id[i], ".txt")), n = 10)
  if (any(grepl("Error", process_tail))) {
    cat("Process log:\n", process_tail, sep = "\n")
  }
}
