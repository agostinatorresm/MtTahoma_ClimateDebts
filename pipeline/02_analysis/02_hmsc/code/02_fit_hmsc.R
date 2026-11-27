# =============================================================================
# 02_fit_hmsc.R
# Author: Agostina Torres
#
# Fits the unfitted HMSC models built by 01_prepare_hmsc_inputs.R.
# Fitting only: convergence diagnostics and figures are produced by
# 03_check_convergence.R. Designed to run unattended on a separate machine:
#
#   - The only input is the unfitted-model file
#       data_clean/<stratum>_<pool>_hmsc_models_unfitted.rds
#     copied to the same relative path inside the project on that machine.
#   - Each model is saved as soon as it is fitted, so an interrupted run keeps
#     all models completed so far; existing fits are not overwritten unless
#     `overwrite = TRUE`.
#   - A run log (settings, machine, R and Hmsc versions, run times) is written
#     next to the fitted models, together with a text log of the console.
#   - Can be run interactively or from the command line:
#       Rscript 02_fit_hmsc.R <run_mode> [stratum] [model_ids] [output_subdir]
#     e.g. Rscript 02_fit_hmsc.R short understory
#          Rscript 02_fit_hmsc.R medium overstory m1_time_elev_pa
#     Giving model_ids (comma-separated) allows each model to be fitted in a
#     separate session; the console log then carries the model id in its name.
#
# MCMC runs are organised in increasing thinning levels (pilot -> short ->
# medium -> long); run time scales approximately with thin x samples, so the
# run time of one level is used to anticipate the next.
#
# Parallelisation: models are fitted sequentially; within each model the
# chains run in parallel (nParallel = nChains). Nested parallelisation
# (several models in parallel, each with parallel chains) is avoided.
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/<stratum>_<pool>_hmsc_models_unfitted.rds
# Outputs: pipeline/02_analysis/02_hmsc/outputs/[<output_subdir>/]
#            (output_subdir, e.g. "run_oct09", keeps each set of runs in its
#            own folder; 03 and 04 read it through their fit_subdir argument)
#            <stratum>_<pool>_<run_mode>_<model_id>.rds   fitted Hmsc object
#            <stratum>_<pool>_<run_mode>_run_log.csv      settings, machine, versions, run times
#            <stratum>_<pool>_<run_mode>_console_log[_<model_ids>].txt  console output of the run
# =============================================================================

library(here)
library(Hmsc)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
stratum   <- "overstory"   # "overstory" | "understory"
pool      <- "strict"
run_mode  <- "short"       # "pilot" | "short" | "medium" | "long"
overwrite <- FALSE         # TRUE = refit models already saved for this run mode

# Models to fit per layer (model_id from 01_prepare_hmsc_inputs.R), fitted
# sequentially in this session
models_by_stratum <- list(
  overstory  = c("m1_time_elev_nf2_nsp_pa", "m1_time_elev_nf2_nsp_rel",
                 "m1_time_elev_nf2_nsp_abu", "m1_time_elev_nf2_nsp_relcnt",
                 "m1_time_elev_nf2_nsp_abucnt",
                 "m0_time_nf2_nsp_pa", "m0_time_nf2_nsp_rel"),
  understory = c("m1_time_elev_nsp_pa", "m1_time_elev_nsp_rel", "m1_time_elev_nsp_abu",
                 "m0_time_nsp_pa", "m0_time_nsp_rel")
)

# Run mode, layer and models from the command line, if given (override the
# values above), e.g. Rscript 02_fit_hmsc.R medium overstory m1_time_elev_pa
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) >= 1) run_mode      <- cmd_args[1]
if (length(cmd_args) >= 2) stratum       <- cmd_args[2]
models_to_fit <- models_by_stratum[[stratum]]
if (length(cmd_args) >= 3) models_to_fit <- strsplit(cmd_args[3], ",")[[1]]
output_subdir <- if (length(cmd_args) >= 4) cmd_args[4] else ""

# Console log per session: when models are given on the command line, the
# model ids are added to the log name so parallel sessions do not share a file
log_suffix <- if (length(cmd_args) >= 3) paste0("_", paste(models_to_fit, collapse = "_")) else ""

# MCMC settings per run mode: samples recorded per chain, thinning interval,
# transient (burn-in) iterations = samples x thin / 2
mcmc_settings <- data.frame(
  run_mode = c("pilot", "short", "medium", "long"),
  samples  = c(100, 250, 250, 250),
  thin     = c(1, 10, 100, 1000)
)
mcmc_settings$transient <- mcmc_settings$samples * mcmc_settings$thin / 2

n_chains <- 4
seed     <- 1

if (!run_mode %in% mcmc_settings$run_mode) {
  stop("Unknown run_mode: ", run_mode, ". Use one of: ",
       paste(mcmc_settings$run_mode, collapse = ", "))
}
settings <- mcmc_settings[mcmc_settings$run_mode == run_mode, ]

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
module_dir <- here("pipeline", "02_analysis", "02_hmsc")
run_tag    <- paste(stratum, pool, sep = "_")
input_file <- file.path(module_dir, "data_clean", paste0(run_tag, "_hmsc_models_unfitted.rds"))
model_dir  <- file.path(module_dir, "outputs", output_subdir)

# File names always start with <stratum>_<pool>_<run_mode>
run_file <- function(name) file.path(model_dir, paste0(run_tag, "_", run_mode, "_", name))
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Unfitted model file not found: ", input_file,
       "\nCopy it from data_clean/ of the main project to the same relative path.")
}

# Console output is mirrored to a text log
# (closed at the end of the script; after an interruption run close_log())
log_con <- file(run_file(paste0("console_log", log_suffix, ".txt")), open = "a")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
close_log <- function() {
  if (sink.number(type = "message") > 2) sink(type = "message")
  while (sink.number() > 0) sink()
  if (isOpen(log_con)) close(log_con)
}

cat("\n==== Run started:", format(Sys.time()), "====\n")
cat("Machine:", Sys.info()[["nodename"]], "| R", R.version$version.string,
    "| Hmsc", as.character(packageVersion("Hmsc")), "\n")
cat("Stratum:", stratum, "| pool:", pool, "| run mode:", run_mode,
    "| samples:", settings$samples, "| thin:", settings$thin,
    "| transient:", settings$transient, "| chains:", n_chains, "\n")

# -----------------------------------------------------------------------------
# 2. Load unfitted models
# -----------------------------------------------------------------------------
models_unfitted <- readRDS(input_file)

missing_models <- setdiff(models_to_fit, models_unfitted$model_id)
if (length(missing_models) > 0) {
  stop("Models not found in ", basename(input_file), ": ",
       paste(missing_models, collapse = ", "),
       "\nAvailable: ", paste(models_unfitted$model_id, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 3. Fit models (sequentially; chains in parallel within each model)
# -----------------------------------------------------------------------------
run_log_file <- run_file("run_log.csv")

for (id in models_to_fit) {
  out_file <- run_file(paste0(id, ".rds"))
  if (file.exists(out_file) && !overwrite) {
    cat("\n", id, ": already fitted for this run mode, skipped ",
        "(set overwrite = TRUE to refit)\n", sep = "")
    next
  }

  m <- models_unfitted$hmsc[[which(models_unfitted$model_id == id)]]
  is_probit <- all(m$distr[, "family"] == 2)

  cat("\nFitting", id, "started at", format(Sys.time()), "\n")
  set.seed(seed)
  t_start <- Sys.time()

  m_fit <- sampleMcmc(m,
                      samples   = settings$samples,
                      thin      = settings$thin,
                      transient = settings$transient,
                      nChains   = n_chains,
                      nParallel = n_chains,
                      initPar   = if (is_probit) "fixed effects" else NULL,
                      verbose   = max(1, floor((settings$samples * settings$thin +
                                                  settings$transient) / 10)))

  elapsed_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  saveRDS(m_fit, out_file)
  cat(id, "fitted in", round(elapsed_min, 1), "min; saved to", out_file, "\n")

  log_row <- data.frame(
    model_id    = id,
    stratum     = stratum,
    pool        = pool,
    run_mode    = run_mode,
    samples     = settings$samples,
    thin        = settings$thin,
    transient   = settings$transient,
    n_chains    = n_chains,
    seed        = seed,
    elapsed_min = round(elapsed_min, 2),
    finished_at = format(Sys.time(), "%Y-%m-%d %H:%M"),
    machine     = Sys.info()[["nodename"]],
    r_version   = R.version$version.string,
    hmsc_version = as.character(packageVersion("Hmsc")),
    file        = basename(out_file)
  )
  write.table(log_row, run_log_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(run_log_file), append = file.exists(run_log_file))
}

# -----------------------------------------------------------------------------
# 4. Expected run time of the next run mode
# -----------------------------------------------------------------------------
run_log <- read.csv(run_log_file)
run_log <- run_log[run_log$run_mode == run_mode, ]
print(run_log[, c("model_id", "run_mode", "elapsed_min", "finished_at")])

i_next <- match(run_mode, mcmc_settings$run_mode) + 1
if (i_next <= nrow(mcmc_settings)) {
  nxt <- mcmc_settings[i_next, ]
  iter_now  <- settings$samples * settings$thin + settings$transient
  iter_next <- nxt$samples * nxt$thin + nxt$transient
  expected  <- data.frame(model_id     = run_log$model_id,
                          next_mode    = nxt$run_mode,
                          expected_min = round(run_log$elapsed_min * iter_next / iter_now, 1),
                          expected_h   = round(run_log$elapsed_min * iter_next / iter_now / 60, 2))
  cat("\nExpected run time for the next run mode (", nxt$run_mode, "):\n", sep = "")
  print(expected)
}

cat("\n==== Run finished:", format(Sys.time()), "====\n")
close_log()
