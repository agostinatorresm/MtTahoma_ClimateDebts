# =============================================================================
# 27_launch_sparse12_reduced.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine)
# Author: Agostina Torres
#
# Sensitivity of the understory occurrence model to sparsely sampled species.
# The model is refitted without 12 species: the 10 sparse species
# (occurrence depends on sampled area; 04_robustness) plus Veratrum viride and
# Epilobium angustifolium, with at least half of their present occurrences only in the 50-m2 subplot (S3.2).
# As for the main occurrence model, it is fitted to the same 10 past datasets
# reduced to the present sampling effort (run_oct09), so that the check differs
# from the main model only in the species removed. The 10 fits are pooled
# afterwards as for the main model (22_pool_reduced_checkpoints.R).
#
#   1. Builds the 10 models from the base occurrence model of the unfitted-model
#      file: species removed from Y, traits and phylogeny; past rows of Y
#      replaced by each reduced dataset. Models are added to that file.
#   2. Writes one lane script per group of 4 cores (02_fit_hmsc.R, CPU).
#   3. Launches the lanes.
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/understory_strict_hmsc_models_unfitted.rds
#          pipeline/02_analysis/02_hmsc/data_clean/understory_reduced_datasets_oct09.csv
#          pipeline/04_robustness/outputs/rarefaction/rarefied_past_occurrence.rds
#          pipeline/04_robustness/outputs/**/area_sensitivity_species.csv
# Outputs: models added to the unfitted-model file;
#          pipeline/02_analysis/02_hmsc/outputs/run_sparse12/ (fits, lane scripts and logs)
# =============================================================================

library(here)
library(Hmsc)
library(dplyr)
library(ape)

stratum    <- "understory"
pool       <- "strict"
out_subdir <- "run_sparse12"
run_mode   <- "medium"    # as the other robustness checks of S5
n_lanes    <- 10          # lanes of 4 cores (10 lanes: all 10 fits in parallel, ~4 h)
launch     <- TRUE
extra_sparse <- c("Veratrum viride", "Epilobium angustifolium")

module_dir    <- here("pipeline", "02_analysis", "02_hmsc")
unfitted_file <- file.path(module_dir, "data_clean", paste0(stratum, "_", pool, "_hmsc_models_unfitted.rds"))
oct09_file    <- file.path(module_dir, "data_clean", "understory_reduced_datasets_oct09.csv")
reduced_file  <- here("pipeline", "04_robustness", "outputs", "rarefaction", "rarefied_past_occurrence.rds")
area_file     <- list.files(here("pipeline", "04_robustness", "outputs"), "^area_sensitivity_species\\.csv$",
                            recursive = TRUE, full.names = TRUE)[1]
out_dir  <- file.path(module_dir, "outputs", out_subdir)
lane_dir <- file.path(out_dir, "lanes")
dir.create(lane_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c(unfitted_file, oct09_file, reduced_file, area_file)) if (is.na(f) || !file.exists(f)) stop("Input not found: ", f)

sp_key  <- function(x) gsub("[ _]+", "_", trimws(x))
base_pa <- "m1_time_elev_nsp_nos_pa"
models  <- readRDS(unfitted_file)
m_pa    <- models$hmsc[[which(models$model_id == base_pa)]]

# -----------------------------------------------------------------------------
# 1. Species removed and models
# -----------------------------------------------------------------------------
area <- read.csv(area_file)
flag_col <- grep("sensitive", names(area), value = TRUE, ignore.case = TRUE)[1]
sp_col   <- grep("^species$", names(area), value = TRUE, ignore.case = TRUE)[1]
area_sp  <- if (!is.na(flag_col) && is.logical(area[[flag_col]])) area[[sp_col]][area[[flag_col]] %in% TRUE] else area[[sp_col]]
sp_all    <- colnames(m_pa$Y)
# the area-sensitivity screening covers all understory species; keep those in the model
cat("Sparse species in the screening (", length(area_sp), "):", paste(area_sp, collapse = ", "), "\n")
area_sp <- area_sp[sp_key(area_sp) %in% sp_key(sp_all)]
cat("Sparse species in the model (", length(area_sp), "):", paste(area_sp, collapse = ", "), "\n")
if (length(area_sp) != 10) stop("Expected 10 sparse species in the model, found ", length(area_sp))
remove_sp <- unique(c(area_sp, extra_sparse))
drop_cols <- which(sp_key(sp_all) %in% sp_key(remove_sp))
missing   <- remove_sp[!sp_key(remove_sp) %in% sp_key(sp_all)]
if (length(missing)) stop("Species not in the model: ", paste(missing, collapse = ", "))
cat("Species removed (", length(drop_cols), "):", paste(sp_all[drop_cols], collapse = ", "), "\n")
keep_sp <- sp_all[-drop_cols]

oct09   <- read.csv(oct09_file)
reduced <- readRDS(reduced_file)
if (!"weighting" %in% names(reduced)) reduced$weighting <- "cover"
reduced <- reduced[reduced$weighting == "cover", ]

time <- as.character(m_pa$XData$time); plot <- as.character(m_pa$studyDesign$plot)
past_rows <- which(time == "past")
tr   <- m_pa$TrData[match(keep_sp, rownames(m_pa$TrData)), , drop = FALSE]
tree <- drop.tip(m_pa$phyloTree, setdiff(m_pa$phyloTree$tip.label, keep_sp))
stopifnot("Traits and phylogeny must match the kept species" = !anyNA(tr) | nrow(tr) == length(keep_sp),
          setequal(tree$tip.label, keep_sp))
base_row <- models[which(models$model_id == base_pa), ]

ids <- sprintf("m1_time_elev_nsp_nos_sparse12_red%03d_pa", oct09$reduced_rep)
new_rows <- lapply(seq_along(ids), function(k) {
  rr  <- reduced[reduced$rep == oct09$reduced_rep[k], ]
  key <- paste(rr$plot_id, sp_key(rr$species))
  Y <- m_pa$Y[, keep_sp, drop = FALSE]
  for (j in seq_along(keep_sp)) Y[past_rows, j] <- as.integer(paste(plot[past_rows], sp_key(keep_sp[j])) %in% key)
  h <- Hmsc(Y = Y, XData = m_pa$XData, XFormula = m_pa$XFormula, TrData = tr, TrFormula = m_pa$TrFormula,
            phyloTree = tree, studyDesign = m_pa$studyDesign, ranLevels = m_pa$ranLevels, distr = "probit")
  row <- base_row; row$model_id <- ids[k]; row$structure <- paste0(row$structure, "_sparse12_reduced"); row$hmsc <- list(h)
  row
})
cat("Past occurrences per reduced dataset (kept species):",
    paste(sapply(new_rows, function(r) sum(r$hmsc[[1]]$Y[past_rows, ])), collapse = ", "), "\n")
models <- bind_rows(filter(models, !model_id %in% ids), bind_rows(new_rows))
saveRDS(models, unfitted_file)
write.csv(data.frame(model_id = ids, reduced_rep = oct09$reduced_rep, fit_subdir = out_subdir),
          file.path(module_dir, "data_clean", "understory_sparse12_reduced_datasets.csv"), row.names = FALSE)
cat("Models added to", basename(unfitted_file), ":", length(ids), "\n")

# -----------------------------------------------------------------------------
# 2. Lanes (CPU, 4 cores each)
# -----------------------------------------------------------------------------
fit_script <- file.path(module_dir, "code", "02_fit_hmsc.R")
stopifnot("02_fit_hmsc.R not found" = file.exists(fit_script))
rscript <- file.path(R.home("bin"), "Rscript")
lanes <- split(ids, sprintf("lane%02d", ((seq_along(ids) - 1) %% n_lanes) + 1))
lane_files <- sapply(names(lanes), function(l) {
  lines <- c("# Lane script written by 27_launch_sparse12_reduced.R",
             sprintf("rscript <- %s", deparse(rscript)), sprintf("fit_script <- %s", deparse(fit_script)))
  for (id in lanes[[l]]) {
    log_file <- file.path(lane_dir, paste0(l, "_", run_mode, "_", id, ".txt"))
    lines <- c(lines, sprintf('system2(rscript, c(shQuote(fit_script), "%s", "%s", "%s", "%s"), stdout = %s, stderr = %s, wait = TRUE)',
                              run_mode, stratum, id, out_subdir, deparse(log_file), deparse(log_file)))
  }
  f <- file.path(lane_dir, paste0(l, ".R")); writeLines(lines, f); f
})
cat("Lanes:", n_lanes, "x 4 cores |", max(lengths(lanes)), "fits per lane\n")

# -----------------------------------------------------------------------------
# 3. Launch
# -----------------------------------------------------------------------------
if (launch) for (f in lane_files) {
  log <- sub("\\.R$", "_lane_log.txt", f)
  system2(rscript, shQuote(f), stdout = log, stderr = log, wait = FALSE)
  cat("Launched", basename(f), "\n"); Sys.sleep(5)
}

# Progress (run at any time)
done <- file.exists(file.path(out_dir, paste0(stratum, "_", pool, "_", run_mode, "_", ids, ".rds")))
cat("\nFinished fits:", sum(done), "of", length(ids), "\n")
