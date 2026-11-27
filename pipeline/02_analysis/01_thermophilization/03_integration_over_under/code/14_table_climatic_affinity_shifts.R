# =============================================================================
# 14_table_climatic_affinity_shifts.R
#
# Pipeline: pipeline/02_analysis/01_thermophilization/03_integration_over_under/code/
#           (after 07_rate_models_bayes.R)
# Author: Agostina Torres
#
# Table S4.1 of the Supporting Information: rates of change in the community
# temperature (CTI) and moisture (CMI) indices, community tails, expected rates
# and climatic lags, and their relationships with the rate of climate change,
# elevation and overstory basal area, for overstory and understory communities,
# and the differences between layers.
#
# All values are read from the outputs of the Bayesian rate models
# (07_rate_models_bayes.R), so the table is rebuilt whenever those models are
# refitted. Each cell gives the posterior mean, the 95% credible interval and
# the probability of direction (f). The expected rate is the mean observed rate
# plus the mean lag (lag = expected − observed), so it has no interval of its
# own. Slopes are per standard deviation of the predictor.
#
# Inputs : outputs/rate_models_bayes_coefficients.csv (07)
#          outputs/rate_models_bayes_strata_differences.csv (07)
# Outputs: outputs/table_S4_1_climatic_affinity_shifts.csv
#          outputs/table_S4_1_climatic_affinity_shifts.docx
# =============================================================================

library(tidyverse)
library(here)

out_dir <- here("pipeline", "02_analysis", "01_thermophilization", "03_integration_over_under", "outputs")
coef <- read_csv(file.path(out_dir, "rate_models_bayes_coefficients.csv"), show_col_types = FALSE)
diff <- read_csv(file.path(out_dir, "rate_models_bayes_strata_differences.csv"), show_col_types = FALSE)

# 1. Formatting ----------------------------------------------------------------
# CTI in °C per decade (three decimals; expected rate and lag two), CMI in mm
# per decade (one decimal; expected rate and lag none). Positive estimates
# carry a "+" except lags, whose sign is part of their definition.
fmt_num <- function(x, digits, plus = TRUE) {
  s <- formatC(abs(x), format = "f", digits = digits)
  sign <- ifelse(x < 0, "−", ifelse(plus, "+", ""))
  paste0(sign, s)
}
fmt_ci <- function(m, lo, hi, digits, plus = TRUE) {
  paste0(fmt_num(m, digits, plus), " [", fmt_num(lo, digits, FALSE), ", ", fmt_num(hi, digits, FALSE), "]")
}

# One coefficient: model label and term ("Intercept" for mean models, the
# predictor name for slopes)
get_coef <- function(model_label, term_name) {
  r <- coef %>% filter(model == model_label, term == term_name)
  if (nrow(r) != 1) stop("Coefficient not found (or not unique): ", model_label, " / ", term_name)
  r
}

row_ci <- function(layer, index, response, model_label, term_name, digits, plus = TRUE) {
  r <- get_coef(model_label, term_name)
  tibble(Layer = layer, Index = index, Response = response,
         Estimate = fmt_ci(r$mean, r$lower_95, r$upper_95, digits, plus),
         f = formatC(r$pd, format = "f", digits = 2))
}

# 2. Rows per layer and index ----------------------------------------------------
build_block <- function(st, v) {
  layer  <- str_to_sentence(st)
  d_main <- if (v == "CTI") 3 else 1            # rates, tails and slopes
  d_lag  <- if (v == "CTI") 2 else 0            # expected rate and lag
  tails  <- if (v == "CTI") c("cold", "warm") else c("dry", "wet")
  clim   <- if (v == "CTI") "warming rate" else "precipitation trend"
  m      <- function(x) paste(st, v, x)
  
  obs <- get_coef(m("A mean rate"), "Intercept")$mean
  lag <- get_coef(m("A2 mean lag"), "Intercept")$mean
  
  bind_rows(
    row_ci(layer, v, "Rate of change", m("A mean rate"), "Intercept", d_main),
    row_ci(layer, v, paste0(str_to_sentence(tails[1]), " tail (10th percentile)"),
           m(paste0("G ", tails[1], " edge (q10) rate")), "Intercept", d_main),
    row_ci(layer, v, paste0(str_to_sentence(tails[2]), " tail (90th percentile)"),
           m(paste0("G ", tails[2], " edge (q90) rate")), "Intercept", d_main),
    row_ci(layer, v, paste0(str_to_sentence(tails[1]), " − ", tails[2], " tail"),
           m("G q10 - q90"), "Intercept", d_main),
    tibble(Layer = layer, Index = v, Response = "Expected rate",
           Estimate = fmt_num(obs + lag, d_lag), f = ""),
    row_ci(layer, v, "Lag (expected − observed)", m("A2 mean lag"), "Intercept", d_lag, plus = FALSE),
    row_ci(layer, v, paste("Rate ~", clim), m("C climate change rate"), "clim_sc", d_main),
    row_ci(layer, v, paste("Lag ~", clim), m("C2 lag ~ climate rate"), "clim_sc", d_main),
    row_ci(layer, v, "Rate ~ elevation", m("B elevation"), "elev_sc", d_main),
    row_ci(layer, v, "Lag ~ elevation", m("D lag ~ elevation"), "elev_sc", d_main)
  )
}

# Understory only: relationships with overstory basal area
basal_area_rows <- function(v) {
  d_main <- if (v == "CTI") 3 else 1
  bind_rows(
    row_ci("Understory", v, "Rate ~ past basal area",
           paste("understory", v, "F rate ~ past basal area"), "ba_past_sc", d_main),
    row_ci("Understory", v, "Rate ~ change in basal area",
           paste("understory", v, "F rate ~ change in basal area"), "d_ba_sc", d_main)
  )
}

# Differences between layers (understory − overstory)
diff_rows <- map_dfr(c("CTI", "CMI"), function(v) {
  d_main <- if (v == "CTI") 3 else 1
  pick <- function(what) diff %>% filter(comparison == paste(v, what, "(understory - overstory)"))
  r_rate <- pick("mean rate"); r_slope <- pick("elevational slope")
  if (nrow(r_rate) != 1 || nrow(r_slope) != 1) stop("Layer difference not found for ", v)
  tibble(Layer = "Understory − overstory", Index = v,
         Response = c("Rate of change", "Slope with elevation"),
         Estimate = c(fmt_ci(r_rate$mean, r_rate$lower_95, r_rate$upper_95, d_main),
                      fmt_ci(r_slope$mean, r_slope$lower_95, r_slope$upper_95, d_main)),
         f = formatC(c(r_rate$pd, r_slope$pd), format = "f", digits = 2))
})

tab <- bind_rows(
  build_block("overstory", "CTI"),
  build_block("overstory", "CMI"),
  build_block("understory", "CTI"), basal_area_rows("CTI"),
  build_block("understory", "CMI"), basal_area_rows("CMI"),
  diff_rows
)

write_csv(tab, file.path(out_dir, "table_S4_1_climatic_affinity_shifts.csv"))
print(tab, n = Inf, width = Inf)

# 3. Word table (requires flextable and officer) --------------------------------
# Horizontal rules above and below the header, between layers and at the bottom
# of the body; layer and index shown only on the first row of each group.
docx_file <- file.path(out_dir, "table_S4_1_climatic_affinity_shifts.docx")
if (requireNamespace("flextable", quietly = TRUE) && requireNamespace("officer", quietly = TRUE)) {
  rule_thick <- officer::fp_border(color = "black", width = 1)
  rule_thin  <- officer::fp_border(color = "black", width = 0.5)
  layer_ends <- which(tab$Layer != lead(tab$Layer, default = last(tab$Layer)))
  index_ends <- setdiff(which(paste(tab$Layer, tab$Index) != lead(paste(tab$Layer, tab$Index), default = "")),
                        c(layer_ends, nrow(tab)))
  
  ft <- flextable::flextable(tab) %>%
    flextable::merge_v(j = "Layer") %>%
    flextable::merge_v(j = c("Layer", "Index"), target = "Index", combine = TRUE) %>%
    flextable::border_remove() %>%
    flextable::hline_top(border = rule_thick, part = "header") %>%
    flextable::hline_bottom(border = rule_thin, part = "header") %>%
    flextable::hline(i = layer_ends, border = rule_thin, part = "body") %>%
    flextable::hline(i = index_ends, j = c("Index", "Response", "Estimate", "f"),
                     border = officer::fp_border(color = "grey70", width = 0.5), part = "body") %>%
    flextable::hline_bottom(border = rule_thick, part = "body") %>%
    flextable::set_header_labels(Estimate = "Estimate [95% CI]") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = "Layer", part = "body") %>%
    flextable::align(align = "left", part = "all") %>%
    flextable::align(j = c("Estimate", "f"), align = "center", part = "all") %>%
    flextable::valign(valign = "top", part = "body") %>%
    flextable::padding(padding.top = 2, padding.bottom = 2, part = "all") %>%
    flextable::width(j = c("Layer", "Index", "Response", "Estimate", "f"),
                     width = c(1.3, 0.5, 2.2, 1.8, 0.5)) %>%
    flextable::fix_border_issues()
  
  doc <- officer::read_docx() %>% flextable::body_add_flextable(ft)
  print(doc, target = docx_file)
  cat("Word table written to", docx_file, "\n")
  open_file <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
  open_file(docx_file)
} else {
  message("Install flextable and officer to write the Word table; the csv was written")
}