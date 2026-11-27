# =============================================================================
# 34_cover_denominator_paired_check.R
#
# Pipeline: pipeline/02_analysis/02_hmsc/code/ (local computer; medium run, 4 cores)
# Author: Agostina Torres
#
# Sensitivity of the understory relative-cover model to the denominator of
# relative cover. In the main model, relative cover is computed over all pool
# species with cover in the sample and then restricted to the species with
# cover in both surveys of the plot (paired values). Because fewer species had
# cover per plot in the present survey (Table S3.1), this denominator is smaller
# in the present and relative cover increases for most species. Here the
# denominator is computed after pairing, over the paired species only, so that
# past and present relative cover of a plot are computed over the same set of
# species.
#
# Both versions are fitted with identical settings (medium run) so that they
# differ only in the denominator:
#   current : Y of the main relative-cover model (unfitted-model file)
#   paired  : cover of each paired species divided by the total cover of the
#             paired species in the same sample
# Comparison (species response to time, at mean elevation): posterior mean of
# the time effect per species, centred on the mean across species as in the
# main analysis (only differences among species are interpreted); correlation
# between versions, number of species with a positive raw and centred change,
# and number whose supported direction (f >= 0.90) differs.
#
# Inputs : pipeline/02_analysis/02_hmsc/data_clean/understory_strict_hmsc_models_unfitted.rds
#          pipeline/02_analysis/02_hmsc/data_clean/understory_strict_hmsc_inputs.rds (Y_abu: log10 cover)
# Outputs: pipeline/02_analysis/02_hmsc/outputs/run_denominator/
#            understory_strict_medium_<model_id>_{current,paired}.rds
#            denominator_species_change.csv
#            denominator_comparison_summary.csv
#            denominator_species_change.png
# =============================================================================

library(here)
library(Hmsc)
library(tidyverse)

stratum  <- "understory"
pool     <- "strict"
model_id <- "m1_time_elev_nsp_nos_rel"   # main relative-cover model
support_thr <- 0.90

# Medium run, as the other robustness checks of S5
mcmc <- list(samples = 250, thin = 100, transient = 12500, nChains = 4, nParallel = 4)
set.seed(20261119)

module_dir <- here("pipeline", "02_analysis", "02_hmsc")
out_dir    <- file.path(module_dir, "outputs", "run_denominator")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fit_file   <- function(version) file.path(out_dir, paste0(stratum, "_", pool, "_medium_", model_id, "_", version, ".rds"))

# -----------------------------------------------------------------------------
# 1. Main relative-cover model and raw cover
# -----------------------------------------------------------------------------
unfitted <- readRDS(file.path(module_dir, "data_clean", paste0(stratum, "_", pool, "_hmsc_models_unfitted.rds")))
if (!model_id %in% unfitted$model_id) stop("Model not found in the unfitted-model file: ", model_id)
m_cur <- unfitted$hmsc[[which(unfitted$model_id == model_id)]]

inputs <- readRDS(file.path(module_dir, "data_clean", paste0(stratum, "_", pool, "_hmsc_inputs.rds")))
if (is.null(inputs$Y_abu)) stop("Y_abu (log10 cover) not found in the inputs file")
cover <- 10^inputs$Y_abu                       # cover, NA where absent or without cover
samples <- rownames(m_cur$Y)
cover <- cover[samples, , drop = FALSE]
plot_of <- as.character(m_cur$studyDesign$plot)

# -----------------------------------------------------------------------------
# 2. Relative cover with the denominator computed after pairing
# -----------------------------------------------------------------------------
# A species is paired in a plot when it has cover in both samples of the plot.
# The denominator is computed over all paired pool species (also those not
# modelled because of too few cover values), then restricted to the modelled
# species.
has_cover <- !is.na(cover)
paired <- has_cover
for (p in unique(plot_of)) {
  rows <- which(plot_of == p)
  both <- colSums(has_cover[rows, , drop = FALSE]) == length(rows)
  paired[rows, ] <- matrix(both, nrow = length(rows), ncol = ncol(cover), byrow = TRUE) & has_cover[rows, , drop = FALSE]
}
cover_paired <- ifelse(paired, cover, NA)
rel_paired   <- cover_paired / rowSums(cover_paired, na.rm = TRUE)
Y_new <- log10(rel_paired[, colnames(m_cur$Y), drop = FALSE])

# The paired cells must be the same as in the main model; only their values change
same_cells <- identical(is.na(Y_new), is.na(m_cur$Y))
cat("Same paired observations as the main model:", same_cells, "\n")
if (!same_cells) {
  cat("Cells differing in pairing:", sum(is.na(Y_new) != is.na(m_cur$Y)), "of", length(Y_new), "\n")
}
# Shared change in log10 relative cover (present - past) in each version
period <- as.character(m_cur$XData[[grep("^(time|period)", names(m_cur$XData), value = TRUE)[1]]])
shared_change <- function(Y) mean(Y[period == "present", ], na.rm = TRUE) - mean(Y[period == "past", ], na.rm = TRUE)
cat("Mean change in log10 relative cover: current", round(shared_change(m_cur$Y), 3),
    "| paired denominator", round(shared_change(Y_new), 3), "\n")

# -----------------------------------------------------------------------------
# 3. Fits (same structure and settings; only Y differs)
# -----------------------------------------------------------------------------
m_new <- Hmsc(Y = Y_new, XData = m_cur$XData, XFormula = m_cur$XFormula,
              TrData = m_cur$TrData, TrFormula = m_cur$TrFormula,
              phyloTree = m_cur$phyloTree, studyDesign = m_cur$studyDesign,
              ranLevels = m_cur$ranLevels, distr = "normal", YScale = TRUE)

fit <- function(m, version) {
  f <- fit_file(version)
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  m <- sampleMcmc(m, samples = mcmc$samples, thin = mcmc$thin, transient = mcmc$transient,
                  nChains = mcmc$nChains, nParallel = mcmc$nParallel, verbose = 0)
  saveRDS(m, f)
  cat(version, "fitted in", round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 1), "h\n")
  m
}
fits <- list(current = fit(m_cur, "current"), paired = fit(m_new, "paired"))

# -----------------------------------------------------------------------------
# 4. Comparison of species responses to time
# -----------------------------------------------------------------------------
species_time <- function(m, version) {
  post <- poolMcmcChains(m$postList)
  it <- which(m$covNames == "timepresent")
  b <- t(sapply(post, function(s) s$Beta[it, ]))            # samples x species
  b_c <- b - rowMeans(b)                                    # centred on the mean across species
  tibble(species = m$spNames, version = version,
         mean_raw = colMeans(b), p_pos_raw = colMeans(b > 0),
         mean_centred = colMeans(b_c), p_pos_centred = colMeans(b_c > 0))
}
sp <- bind_rows(species_time(fits$current, "current"), species_time(fits$paired, "paired")) %>%
  mutate(direction = case_when(p_pos_centred >= support_thr ~ "increase",
                               1 - p_pos_centred >= support_thr ~ "decrease",
                               TRUE ~ "not supported"))
write_csv(sp, file.path(out_dir, "denominator_species_change.csv"))

w <- sp %>% select(species, version, mean_raw, mean_centred, direction) %>%
  pivot_wider(names_from = version, values_from = c(mean_raw, mean_centred, direction))
summary_tab <- tibble(
  n_species               = nrow(w),
  shared_change_current   = round(shared_change(m_cur$Y), 3),
  shared_change_paired    = round(shared_change(Y_new), 3),
  n_positive_raw_current  = sum(w$mean_raw_current > 0),
  n_positive_raw_paired   = sum(w$mean_raw_paired > 0),
  r_centred               = round(cor(w$mean_centred_current, w$mean_centred_paired), 3),
  n_direction_differs     = sum(w$direction_current != w$direction_paired)
)
print(summary_tab, width = Inf)
write_csv(summary_tab, file.path(out_dir, "denominator_comparison_summary.csv"))
cat("\nSpecies whose supported direction differs:\n")
print(w %>% filter(direction_current != direction_paired) %>% select(species, direction_current, direction_paired))

p <- ggplot(w, aes(mean_centred_current, mean_centred_paired)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(colour = direction_current == direction_paired), size = 2) +
  scale_colour_manual(values = c(`TRUE` = "grey30", `FALSE` = "#F26B5B"),
                      labels = c(`TRUE` = "Same direction", `FALSE` = "Direction differs"), name = NULL) +
  labs(x = "Species change, current denominator (centred)",
       y = "Species change, denominator after pairing (centred)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"), legend.position = "top")
ggsave(file.path(out_dir, "denominator_species_change.png"), p, width = 6, height = 6, dpi = 300, bg = "white")

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(file.path(out_dir, "denominator_species_change.png"))


library(here); library(Hmsc); library(tidyverse)
out_dir <- here("pipeline", "02_analysis", "02_hmsc", "outputs", "run_denominator")
fits <- list(current = readRDS(file.path(out_dir, "understory_strict_medium_m1_time_elev_nsp_nos_rel_current.rds")),
             paired  = readRDS(file.path(out_dir, "understory_strict_medium_m1_time_elev_nsp_nos_rel_paired.rds")))

# 1. Trait effects on the time response (Gamma)
gamma_time <- function(m, version) {
  post <- poolMcmcChains(m$postList)
  it <- which(m$covNames == "timepresent")
  g <- t(sapply(post, function(s) s$Gamma[it, ])); colnames(g) <- m$trNames
  tibble(version = version, trait = m$trNames, mean = colMeans(g),
         lo90 = apply(g, 2, quantile, 0.05), hi90 = apply(g, 2, quantile, 0.95),
         support = pmax(colMeans(g > 0), colMeans(g < 0)))
}
gam <- bind_rows(gamma_time(fits$current, "current"), gamma_time(fits$paired, "paired"))
print(gam %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% arrange(trait, version), n = Inf)

# 2. Paired denominator without centring vs current denominator with centring
d <- read_csv(file.path(out_dir, "denominator_species_change.csv"), show_col_types = FALSE) %>%
  select(species, version, mean_raw, mean_centred) %>%
  pivot_wider(names_from = version, values_from = c(mean_raw, mean_centred))
cat("\nr(paired raw, current centred):", round(cor(d$mean_raw_paired, d$mean_centred_current), 3),
    "\nr(current raw, current centred):", round(cor(d$mean_raw_current, d$mean_centred_current), 3), "\n")
