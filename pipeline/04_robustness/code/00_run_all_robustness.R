# =============================================================================
# 00_run_all_robustness.R
# Author: Agostina Torres
#
# Runs the whole robustness module in order, each script in its own
# environment:
#   00 inputs -> 01a accumulation -> 01b patchy -> 02 bootstrap ->
#   03 occurrence sampling -> 04 reliability -> 05 rarefaction ->
#   06 rarefied GLMs -> 07 census-only cover -> 10 raw cover diagnosis ->
#   08 diagnosis figures (Figures S3.1 and S3.3) -> 09 rarefied thermophilization
# 10 needs the HMSC inputs (02_hmsc/code/01_prepare_hmsc_inputs.R) for the
# species set, and 08 reads the panels saved by 01a, 01b, 02, 03 and 10.
# Figures are written to pipeline/04_robustness/figures/.
# The console output of each script is also written to
# pipeline/04_robustness/figures/run_log_<script>.txt.
# =============================================================================

library(here)

code_dir <- here("pipeline", "04_robustness", "code")
log_dir  <- here("pipeline", "04_robustness", "figures")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

scripts <- c("00_robustness_inputs.R",
             "01a_accumulation_curves.R",
             "01b_patchy_species.R",
             "02_bootstrapping.R",
             "03_occurrence_sampling.R",
             "04_reliability_diagnostic.R",
             "05_rarefaction_occurrence.R",
             "06_rarefied_glm_occurrence.R",
             "07_census_only_cover.R",
             "10_raw_cover_diagnosis.R",
             "08_diagnosis_figures.R",
             "09_rarefied_thermophilization.R")

missing <- scripts[!file.exists(file.path(code_dir, scripts))]
if (length(missing) > 0) stop("Scripts not found in ", code_dir, ": ", paste(missing, collapse = ", "))

for (s in scripts) {
  cat("\n=============================", s, "=============================\n")
  t0 <- Sys.time()
  log_file <- file.path(log_dir, paste0("run_log_", sub("\\.R$", "", s), ".txt"))
  sink(log_file, split = TRUE)
  ok <- tryCatch({
    source(file.path(code_dir, s), local = new.env(), echo = FALSE)
    TRUE
  }, error = function(e) {
    cat("ERROR in", s, ":", conditionMessage(e), "\n")
    FALSE
  }, finally = sink())
  graphics.off()
  cat(s, if (ok) "done" else "FAILED", "in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  if (!ok) stop("Stopped at ", s, "; see ", log_file)
}

cat("\nAll scripts finished. Figures and logs in", log_dir, "\n")
print(list.files(log_dir))
