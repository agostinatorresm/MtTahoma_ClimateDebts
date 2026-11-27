# =============================================================================
# 21_launch_reduced100.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine)
# Author: Agostina Torres
#
# Extends the understory occurrence model fitted to past data reduced to the
# present sampling effort from 10 to all 100 cover-weighted reductions
# (04_robustness 05). The 10 reductions already fitted in run_oct09 are kept;
# the remaining 90 are fitted here with the same model structure and run mode
# (m1, ~ time * elevation, non-spatial plot level, long run).
#
# The 90 new reductions are fitted in a random order, so that any set of
# finished fits is a random subset of the 100. Together with the 10 existing
# fits, this allows the pooled posterior to be built at intermediate
# checkpoints (10, 20, 40, 60, 80 and 100 reductions; 22_pool_reduced_checkpoints.R)
# and its stability to be followed as fits accumulate.
#
#   1. Builds the 90 models from the base occurrence model of the unfitted-
#      model file (01_prepare_hmsc_inputs.R) and adds them to that file.
#   2. Writes one lane script per group of 4 cores. Models are assigned to
#      lanes in launch order (round robin), so that fits finish roughly in that
#      order. `n_lanes` can be set below the machine capacity to keep cores free
#      for other runs.
#   3. Launches the lanes.
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/understory_strict_hmsc_models_unfitted.rds
#          pipeline/02_analysis/02_hmsc/data_clean/understory_reduced_datasets_oct09.csv
#          pipeline/04_robustness/outputs/rarefaction/rarefied_past_occurrence.rds
# Outputs: models added to the unfitted-model file;
#          pipeline/02_analysis/02_hmsc/data_clean/understory_reduced_datasets_100.csv
#            (model id, reduction, launch order, fit folder; the 10 earlier fits first)
#          pipeline/02_analysis/02_hmsc/outputs/run_red100/ (fits, run logs, lane scripts and logs)
# =============================================================================

library(here)
library(Hmsc)
library(dplyr)
library(ape)

stratum    <- "understory"
pool       <- "strict"
out_subdir <- "run_red100"
run_mode   <- "long"
n_lanes    <- 10      # lanes of 4 cores; lower it to leave cores for other runs
launch     <- TRUE    # FALSE builds the models and lane scripts without launching
set.seed(20261102)

module_dir    <- here("pipeline", "02_analysis", "02_hmsc")
unfitted_file <- file.path(module_dir, "data_clean", paste0(stratum, "_", pool, "_hmsc_models_unfitted.rds"))
reduced_file  <- here("pipeline", "04_robustness", "outputs", "rarefaction", "rarefied_past_occurrence.rds")
oct09_file    <- file.path(module_dir, "data_clean", "understory_reduced_datasets_oct09.csv")
map_file      <- file.path(module_dir, "data_clean", "understory_reduced_datasets_100.csv")
out_dir       <- file.path(module_dir, "outputs", out_subdir)
lane_dir      <- file.path(out_dir, "lanes")
dir.create(lane_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(unfitted_file, reduced_file, oct09_file)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}

base_pa <- "m1_time_elev_nsp_nos_pa"
sp_key  <- function(x) gsub("[ _]+", "_", trimws(x))

models <- readRDS(unfitted_file)
if (!base_pa %in% models$model_id) stop("Base model not found in the unfitted-model file: ", base_pa)
m_pa <- models$hmsc[[which(models$model_id == base_pa)]]

rebuild <- function(m, Y) {
  Hmsc(Y = Y, XData = m$XData, XFormula = m$XFormula,
       TrData = m$TrData, TrFormula = m$TrFormula, phyloTree = m$phyloTree,
       studyDesign = m$studyDesign, ranLevels = m$ranLevels, distr = "probit")
}
new_row <- function(id, hmsc) {
  row <- models[which(models$model_id == base_pa), ]
  row$model_id  <- id
  row$structure <- paste0(row$structure, "_reduced")
  row$hmsc      <- list(hmsc)
  row
}

# -----------------------------------------------------------------------------
# 1. Reductions: the 10 already fitted, then the other 90 in random order
# -----------------------------------------------------------------------------
reduced <- readRDS(reduced_file)
if (!"weighting" %in% names(reduced)) reduced$weighting <- "cover"
reduced <- reduced[reduced$weighting == "cover", ]
all_reps <- sort(unique(reduced$rep))
stopifnot("Expected 100 cover-weighted reductions" = length(all_reps) == 100)

oct09   <- read.csv(oct09_file)
new_reps <- sample(setdiff(all_reps, oct09$reduced_rep))
new_ids  <- sprintf("m1_time_elev_nsp_nos_red%03d_pa", new_reps)

map <- bind_rows(
  data.frame(model_id = oct09$model_id, reduced_rep = oct09$reduced_rep, fit_subdir = "run_oct09"),
  data.frame(model_id = new_ids, reduced_rep = new_reps, fit_subdir = out_subdir)
) %>% mutate(launch_order = row_number(), .before = 1)
write.csv(map, map_file, row.names = FALSE)
cat("Reductions: 10 fitted in run_oct09 +", length(new_reps), "new =", nrow(map), "\n")

time <- as.character(m_pa$XData$time); plot <- as.character(m_pa$studyDesign$plot)
past_rows <- which(time == "past")
species   <- colnames(m_pa$Y)

new_rows <- lapply(seq_along(new_reps), function(k) {
  rr  <- reduced[reduced$rep == new_reps[k], ]
  key <- paste(rr$plot_id, sp_key(rr$species))
  Y <- m_pa$Y
  for (j in seq_along(species)) {
    Y[past_rows, j] <- as.integer(paste(plot[past_rows], sp_key(species[j])) %in% key)
  }
  new_row(new_ids[k], rebuild(m_pa, Y))
})
n_past_red <- sapply(new_rows, function(r) sum(r$hmsc[[1]]$Y[past_rows, ]))
stopifnot("Reduced past occurrences must not exceed the observed ones" = all(n_past_red <= sum(m_pa$Y[past_rows, ])),
          "Present rows must be unchanged" =
            all(sapply(new_rows, function(r) identical(r$hmsc[[1]]$Y[-past_rows, ], m_pa$Y[-past_rows, ]))))
cat("Past occurrences: observed", sum(m_pa$Y[past_rows, ]), "| reduced (mean):", round(mean(n_past_red)), "\n")

new_rows <- bind_rows(new_rows)
models <- bind_rows(filter(models, !model_id %in% new_rows$model_id), new_rows)
saveRDS(models, unfitted_file)
cat("Models added to", basename(unfitted_file), ":", nrow(new_rows), "\n")

# -----------------------------------------------------------------------------
# 2. Lanes: round robin in launch order, so fits finish roughly in that order
# -----------------------------------------------------------------------------
lane_of <- ((seq_along(new_ids) - 1) %% n_lanes) + 1
lanes <- split(new_ids, sprintf("lane%02d", lane_of))
plan <- bind_rows(lapply(names(lanes), function(l)
  data.frame(lane = l, step = seq_along(lanes[[l]]), run_mode = run_mode, model_id = lanes[[l]])))
write.csv(plan, file.path(out_dir, "launch_plan_red100.csv"), row.names = FALSE)
cat("\nLaunch plan:", n_lanes, "lanes x 4 cores =", 4 * n_lanes, "cores (available:",
    parallel::detectCores(), ") |", max(table(lane_of)), "fits per lane\n")

fit_script <- file.path(module_dir, "code", "02_fit_hmsc.R")
stopifnot("02_fit_hmsc.R not found" = file.exists(fit_script))
rscript <- file.path(R.home("bin"), "Rscript")
q <- function(x) shQuote(x, type = if (.Platform$OS.type == "windows") "cmd" else "sh")

lane_files <- sapply(names(lanes), function(l) {
  lines <- c("# Lane script written by 21_launch_reduced100.R: fits its models one after the other",
             sprintf("rscript <- %s", deparse(rscript)),
             sprintf("fit_script <- %s", deparse(fit_script)))
  for (id in lanes[[l]]) {
    log_file <- file.path(lane_dir, paste0(l, "_", run_mode, "_", id, ".txt"))
    lines <- c(lines,
               sprintf('cat("Start %s:", format(Sys.time()), "\\n")', id),
               sprintf('system2(rscript, c(shQuote(fit_script), "%s", "%s", "%s", "%s"), stdout = %s, stderr = %s, wait = TRUE)',
                       run_mode, stratum, id, out_subdir, deparse(log_file), deparse(log_file)),
               sprintf('cat("End %s:", format(Sys.time()), "\\n")', id))
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
# 4. Progress (run these lines at any time)
# -----------------------------------------------------------------------------
done <- file.exists(file.path(module_dir, "outputs", map$fit_subdir,
                              paste0(stratum, "_", pool, "_", run_mode, "_", map$model_id, ".rds")))
cat("\nFinished fits:", sum(done), "of", nrow(map), "\n")
