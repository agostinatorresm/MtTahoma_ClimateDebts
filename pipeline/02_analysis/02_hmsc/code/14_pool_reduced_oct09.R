# =============================================================================
# 14_pool_reduced_oct09.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (run on the computing machine
#           after the runs of 11_launch_oct09.R; then 04_hmsc_results.R)
# Author: Agostina Torres
#
# Pools the posterior of the understory occurrence models fitted to past
# datasets reduced to the present sampling effort (11_launch_oct09.R, cover-
# weighted removal) into one model object, used as the main estimate of
# change in occurrence.
#
# Each reduced dataset gives a posterior conditional on one random reduction;
# pooling the posterior samples of all fits gives a mixture over reductions,
# which carries both the parameter uncertainty and the uncertainty from the
# reduction. Every fit contributes the same number of samples: each chain is
# thinned to `keep_per_chain` samples, so that the pooled posterior has the
# size of a single fit (n_reduced x 4 chains x keep_per_chain samples) and the
# downstream scripts run unchanged.
#
# Quantities computed per posterior sample (species change, Gamma, Omega =
# Lambda'Lambda, Eta x Lambda) are unaffected by sign or rotation differences
# of the latent factors between fits; Eta and Lambda themselves are not
# summarised from the pooled object.
#
# The pooled object keeps the data of the first reduced fit as Y, plus:
#   Y_observed     : observed occurrence matrix (unreduced), used for plot
#                    occupancy in 04_hmsc_results.R
#   Y_reduced_mean : mean occurrence over the reduced datasets, used as the
#                    observed values in the posterior predictive check
#   pooled_from    : model ids of the pooled fits
#
# Inputs : pipeline/02_analysis/02_hmsc/outputs/run_oct09/understory_strict_long_m1_time_elev_nsp_nos_redNN_pa.rds
#          pipeline/02_analysis/02_hmsc/outputs/run_oct09/understory_strict_long_run_oct09_conv_summary.csv (03)
#          pipeline/02_analysis/02_hmsc/data_clean/understory_strict_hmsc_models_unfitted.rds
# Outputs: pipeline/02_analysis/02_hmsc/outputs/run_oct09/understory_strict_long_m1_time_elev_nsp_nos_redpool_pa.rds
# =============================================================================

library(here)
library(Hmsc)

stratum    <- "understory"
pool       <- "strict"
run_mode   <- "long"
fit_subdir <- "run_oct09"
base_pa    <- "m1_time_elev_nsp_nos_pa"
pooled_id  <- "m1_time_elev_nsp_nos_redpool_pa"
keep_per_chain <- 25        # samples kept per chain (of 250)
psrf_target    <- 1.1

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
fit_dir    <- file.path(module_dir, "outputs", fit_subdir)
fit_file   <- function(id) file.path(fit_dir, paste0(stratum, "_", pool, "_", run_mode, "_", id, ".rds"))

prefix      <- paste0(stratum, "_", pool, "_", run_mode, "_")
reduced_ids <- list.files(fit_dir, pattern = paste0("^", prefix, ".*_red[0-9]{2}_pa\\.rds$"))
reduced_ids <- sort(sub("\\.rds$", "", sub(paste0("^", prefix), "", reduced_ids)))
cat("Reduced fits found:", length(reduced_ids), "\n")
stopifnot("No reduced fits found" = length(reduced_ids) > 0)

# Convergence of the fits to be pooled (03_check_convergence.R)
conv_file <- file.path(module_dir, "outputs", paste0(stratum, "_", pool, "_", run_mode, "_", fit_subdir, "_conv_summary.csv"))
if (file.exists(conv_file)) {
  conv <- read.csv(conv_file)
  conv <- conv[conv$model_id %in% reduced_ids & conv$group %in% c("Beta", "Gamma"), ]
  print(conv[, c("model_id", "group", "psrf_max", "prop_psrf_high")], row.names = FALSE)
  if (any(conv$psrf_max > psrf_target)) warning("Some pooled fits have Beta or Gamma PSRF above ", psrf_target)
} else {
  warning("Convergence summary not found (run 03 for run_oct09 first): ", conv_file)
}

# -----------------------------------------------------------------------------
# Pool
# -----------------------------------------------------------------------------
thin_chain <- function(chain) {
  idx <- round(seq(length(chain) / keep_per_chain, length(chain), length.out = keep_per_chain))
  chain[idx]
}

post_all <- list()
Y_list   <- list()
ref      <- NULL
for (id in reduced_ids) {
  m <- readRDS(fit_file(id))
  if (is.null(ref)) {
    ref <- m
  } else {
    stopifnot(identical(m$spNames, ref$spNames), identical(m$covNames, ref$covNames),
              identical(m$trNames, ref$trNames), identical(dim(m$Y), dim(ref$Y)),
              identical(as.character(m$studyDesign$plot), as.character(ref$studyDesign$plot)),
              identical(m$samples, ref$samples), length(m$postList) == length(ref$postList))
  }
  post_all <- c(post_all, lapply(m$postList, thin_chain))
  Y_list[[id]] <- m$Y
  cat(id, ": chains", length(m$postList), "| samples per chain", m$samples,
      "| past occurrences", sum(m$Y[m$XData$time == "past", ]), "\n")
  rm(m); gc()
}

pooled <- ref
pooled$postList <- post_all
pooled$samples  <- keep_per_chain
pooled$thin     <- ref$thin * (ref$samples / keep_per_chain)
pooled$Y_reduced_mean <- Reduce(`+`, Y_list) / length(Y_list)
pooled$pooled_from    <- reduced_ids

# Observed (unreduced) occurrence matrix from the base model
unfitted <- readRDS(file.path(module_dir, "data_clean", paste0(stratum, "_", pool, "_hmsc_models_unfitted.rds")))
Y_obs <- unfitted$hmsc[[which(unfitted$model_id == base_pa)]]$Y
stopifnot(identical(dim(Y_obs), dim(pooled$Y)), identical(colnames(Y_obs), pooled$spNames),
          "Present data must be identical to the reduced fits" =
            identical(Y_obs[pooled$XData$time == "present", ], pooled$Y[pooled$XData$time == "present", ]))
pooled$Y_observed <- Y_obs

cat("\nPooled posterior:", length(pooled$postList), "chains x", pooled$samples, "samples =",
    length(pooled$postList) * pooled$samples, "samples\n")
cat("Past occurrences: observed", sum(Y_obs[pooled$XData$time == "past", ]),
    "| mean over reduced datasets", round(sum(pooled$Y_reduced_mean[pooled$XData$time == "past", ])), "\n")

saveRDS(pooled, fit_file(pooled_id))
cat("Saved:", fit_file(pooled_id), "\n")
