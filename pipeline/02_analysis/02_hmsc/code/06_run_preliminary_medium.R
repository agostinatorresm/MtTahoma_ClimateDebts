# =============================================================================
# 06_run_preliminary_medium.R
# Author: Agostina Torres
#
# Results of several sets of fitted HMSC models in one go. For each set
# (fixed-effect structure, abundance response and folder of the fits):
#   1. convergence diagnostics (03_check_convergence.R), per layer
#   2. results of the occurrence + abundance models (04_hmsc_results.R), per
#      layer
#   3. combined overstory / understory figures (05_hmsc_figures_combined.R)
#
# Each script is sourced in its own environment, with its settings given as
# the command-line arguments the script reads, so the settings written at the
# top of each script are not used and runs do not share objects. Scripts run
# in the current R session, so figures use the same graphics devices and fonts
# as when the scripts are run on their own.
#
# Sets (edit below):
#   model    : fixed-effect structure per layer (model id without response)
#   response : abundance response analysed with the occurrence model
#   subdir   : folder of the fits inside 02_hmsc/outputs ("" = outputs itself)
#   tag      : added to the names of the combined figures
#   fig_subdir : subfolder of 02_hmsc/figures for all figures of the set
# Outputs: those of 03, 04 and 05.
# =============================================================================

library(here)

run_mode        <- "medium"
run_convergence <- TRUE      # FALSE to skip 03 when its outputs are already up to date
strata          <- c("overstory", "understory")

sets <- list(
  # Folder of the fits: on the computing server, m0 fits are in outputs/ (""); on the local computer,
  # in "transfer_medium".
  # Time-only structure, current inputs, relative abundance (main response)
  m0_rel = list(model    = c(overstory = "m0_time_nf2_nsp", understory = "m0_time_nsp"),
                response = "rel", subdir = "transfer_medium", tag = "m0_rel", fig_subdir = "m0")
)

code_dir <- here("pipeline", "02_analysis", "02_hmsc", "code")
out_dir  <- here("pipeline", "02_analysis", "02_hmsc", "outputs")

# Sources a script with the given command-line arguments; an error stops that
# script only and is reported
run_script <- function(script, args) {
  cat("\n#############################################\n",
      "Running ", script, " ", paste(args, collapse = " "),
      "\n#############################################\n", sep = "")
  env <- new.env(parent = globalenv())
  env$commandArgs <- function(trailingOnly = FALSE) args
  ok <- tryCatch({
    sys.source(file.path(code_dir, script), envir = env, keep.source = FALSE)
    TRUE
  }, error = function(e) {
    cat("\nFAILED:", script, paste(args, collapse = " "), "\n  ", conditionMessage(e), "\n")
    FALSE
  })
  graphics.off()
  ok
}

summary_tab <- list()

for (set_name in names(sets)) {
  st_set <- sets[[set_name]]
  options(mora.fig_subdir = if (is.null(st_set$fig_subdir)) "" else st_set$fig_subdir)
  cat("\n\n==================== SET", set_name, "====================\n")
  
  # Fits needed for this set
  for (st in strata) {
    ids   <- paste0(st_set$model[[st]], c("_pa", paste0("_", st_set$response)))
    files <- file.path(out_dir, st_set$subdir, paste0(st, "_strict_", run_mode, "_", ids, ".rds"))
    cat(st, ":", paste0(basename(files), ifelse(file.exists(files), " [found]", " [MISSING]"), collapse = "; "), "\n")
  }
  
  ok_conv <- if (run_convergence) {
    sapply(strata, function(st) run_script("03_check_convergence.R", c(st, run_mode, st_set$subdir)))
  } else setNames(rep(NA, length(strata)), strata)
  
  ok_res <- sapply(strata, function(st)
    run_script("04_hmsc_results.R", c(st, st_set$response, run_mode, st_set$subdir, st_set$model[[st]])))
  
  ok_fig <- if (all(ok_res)) {
    run_script("05_hmsc_figures_combined.R",
               c(run_mode,
                 paste0(st_set$model[["overstory"]], "_", st_set$response),
                 paste0(st_set$model[["understory"]], "_", st_set$response),
                 st_set$tag))
  } else FALSE
  
  summary_tab[[set_name]] <- c(setNames(ok_conv, paste0("03_", strata)),
                               setNames(ok_res, paste0("04_", strata)),
                               `05_combined` = ok_fig)
}

options(mora.fig_subdir = NULL)

# Summary ----------------------------------------------------------------------------
cat("\n==============================\nSUMMARY (", run_mode, ")\n==============================\n", sep = "")
print(do.call(rbind, summary_tab))

for (set_name in names(sets)) {
  for (st in strata) {
    sub <- sets[[set_name]]$subdir
    f <- file.path(out_dir, paste0(st, "_strict_", run_mode, "_", if (nzchar(sub)) paste0(sub, "_") else "",
                                   "conv_model_comparison.csv"))
    if (file.exists(f)) {
      cat("\n--- convergence:", set_name, st, "---\n")
      print(as.data.frame(read.csv(f)), digits = 3)
    }
  }
}
cat("\nFigures:\n")
print(list.files(here("pipeline", "02_analysis", "02_hmsc", "figures"), recursive = TRUE))
