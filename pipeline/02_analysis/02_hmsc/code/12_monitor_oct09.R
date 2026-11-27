# =============================================================================
# 12_monitor_oct09.R
# Author: Agostina Torres
#
# Progress of the understory HMSC runs launched by 11_launch_oct09.R.
# With chains in parallel the sampler does not report iterations, so progress
# is followed from the files each step writes:
#   - start of a step : creation time of its log file in run_oct09/lanes/
#   - end of a step   : fitted model file in run_oct09/
# Expected durations are approximate (previous understory runs, m1 non-spatial,
# plot level only): long ~33 h, medium ~3.5 h.
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/run_oct09/launch_plan_oct09.csv
#          pipeline/02_analysis/02_hmsc/outputs/run_oct09/lanes/*.txt
#          pipeline/02_analysis/02_hmsc/outputs/run_oct09/*.rds
# Outputs: console only
# =============================================================================

library(here)

stratum <- "understory"
pool    <- "strict"
expected_h <- c(long = 33, medium = 3.5)

out_dir  <- here("pipeline", "02_analysis", "02_hmsc", "outputs", "run_oct09")
lane_dir <- file.path(out_dir, "lanes")
plan <- read.csv(file.path(out_dir, "launch_plan_oct09.csv"))
now  <- Sys.time()

plan$log_file <- file.path(lane_dir, paste0(plan$lane, "_", plan$run_mode, "_", plan$model_id, ".txt"))
plan$fit_file <- file.path(out_dir, paste0(stratum, "_", pool, "_", plan$run_mode, "_", plan$model_id, ".rds"))
plan$started  <- file.exists(plan$log_file)
plan$finished <- file.exists(plan$fit_file)

info_log <- file.info(plan$log_file)
info_fit <- file.info(plan$fit_file)
plan$start_time <- info_log$ctime
plan$end_time   <- info_fit$mtime

plan$status <- ifelse(plan$finished, "done", ifelse(plan$started, "running", "queued"))
end_ref <- plan$end_time
end_ref[!plan$finished] <- now
plan$elapsed_h <- round(as.numeric(difftime(end_ref, plan$start_time, units = "hours")), 1)
plan$elapsed_h[!plan$started] <- NA
plan$expected_end <- format(plan$start_time + expected_h[plan$run_mode] * 3600, "%a %H:%M")
plan$expected_end[!plan$started | plan$finished] <- ""

# Queued steps: expected end after the previous step of the same lane
for (i in which(plan$status == "queued")) {
  prev <- which(plan$lane == plan$lane[i] & plan$step == plan$step[i] - 1)
  prev_end <- if (plan$finished[prev]) plan$end_time[prev] else plan$start_time[prev] + expected_h[plan$run_mode[prev]] * 3600
  plan$expected_end[i] <- format(prev_end + expected_h[plan$run_mode[i]] * 3600, "%a %H:%M")
}

cat("Status at", format(now, "%a %d %b %H:%M"), "\n\n")
print(plan[, c("lane", "step", "run_mode", "model_id", "status", "elapsed_h", "expected_end")], row.names = FALSE)
cat("\nDone:", sum(plan$finished), "of", nrow(plan), "| running:", sum(plan$status == "running"),
    "| queued:", sum(plan$status == "queued"), "\n")

# Errors: a step whose log exists, whose fit is missing, and whose lane has
# moved on to the next step (or whose log reports an error)
read_safe <- function(f) {
  tmp <- tempfile(fileext = ".txt")
  if (!suppressWarnings(file.copy(f, tmp))) return(character(0))
  readLines(tmp, warn = FALSE)
}
for (i in which(plan$started & !plan$finished)) {
  txt <- read_safe(plan$log_file[i])
  err <- grep("Error|error|Execution halted", txt, value = TRUE)
  if (length(err) > 0) cat("\n!!", plan$model_id[i], ":\n", paste(head(err, 5), collapse = "\n"), "\n")
}

# Running R processes (Windows)
if (.Platform$OS.type == "windows") {
  n_proc <- length(grep("Rscript.exe", system('tasklist /FI "IMAGENAME eq Rscript.exe"', intern = TRUE)))
  cat("\nRscript processes running:", n_proc, "\n")
}
