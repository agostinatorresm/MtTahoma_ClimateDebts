# =============================================================================
# 18_table_lag_linear_vs_loess.R
#
# Pipeline: pipeline/02_analysis/01_thermophilization/03_integration_over_under/code/
# Author: Agostina Torres
#
# Table S4.2 of the Supporting Information: climatic lag (expected - observed
# rate) estimated with a linear or a LOESS (span 0.75) space-for-time
# relationship between past community indices and past climate, for overstory
# and understory communities. LOESS cannot predict beyond the range of past
# climate, so plots outside that range are excluded; the number of plots used
# and excluded is given for each estimate. The linear relationship is the one
# used in the main analyses; the LOESS version checks that the lags do not
# depend on assuming a linear relationship.
#
# Values are entered by hand from the lag outputs (02_thermo_xero_elevation_lag.R
# in 01_overstory and 02_understory) and must be updated when those results
# change. The Word table uses the same layout as Table S4.1
# (14_table_climatic_affinity_shifts.R).
#
# Outputs: pipeline/02_analysis/01_thermophilization/03_integration_over_under/outputs/table_s4_2_lag_linear_vs_loess.csv
#          pipeline/02_analysis/01_thermophilization/03_integration_over_under/outputs/table_s4_2_lag_linear_vs_loess.docx
# =============================================================================

library(tidyverse)
library(here)

out_dir <- here("pipeline", "02_analysis", "01_thermophilization", "03_integration_over_under", "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tab <- tribble(
  ~Layer, ~Index, ~Relationship, ~Plots, ~Lag,
  "Overstory",  "CTI (°C decade⁻¹)", "Linear", "324 (0)",   "0.266 [0.234, 0.297]",
  "Overstory",  "CTI (°C decade⁻¹)", "LOESS",  "289 (35)",  "0.236 [0.200, 0.272]",
  "Overstory",  "CMI (mm decade⁻¹)", "Linear", "324 (0)",   "−90.7 [−99.9, −81.6]",
  "Overstory",  "CMI (mm decade⁻¹)", "LOESS",  "205 (119)", "−95.7 [−108, −82.9]",
  "Understory", "CTI (°C decade⁻¹)", "Linear", "316 (0)",   "0.126 [0.092, 0.160]",
  "Understory", "CTI (°C decade⁻¹)", "LOESS",  "283 (33)",  "0.127 [0.092, 0.162]",
  "Understory", "CMI (mm decade⁻¹)", "Linear", "316 (0)",   "−39.4 [−46.5, −32.3]",
  "Understory", "CMI (mm decade⁻¹)", "LOESS",  "200 (116)", "−43.7 [−52.8, −34.6]"
)

write_csv(tab, file.path(out_dir, "table_s4_2_lag_linear_vs_loess.csv"))

# Word table (requires flextable and officer), same layout as Table S4.1:
# horizontal rules above and below the header, between layers and at the
# bottom of the body; layer and index shown only on the first row of each group.
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
    flextable::hline(i = index_ends, j = c("Index", "Relationship", "Plots", "Lag"),
                     border = officer::fp_border(color = "grey70", width = 0.5), part = "body") %>%
    flextable::hline_bottom(border = rule_thick, part = "body") %>%
    flextable::set_header_labels(Plots = "Plots (excluded)", Lag = "Lag, mean [95% CI]") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = "Layer", part = "body") %>%
    flextable::align(align = "left", part = "all") %>%
    flextable::align(j = c("Plots", "Lag"), align = "center", part = "all") %>%
    flextable::valign(valign = "top", part = "body") %>%
    flextable::padding(padding.top = 2, padding.bottom = 2, part = "all") %>%
    flextable::width(j = c("Layer", "Index", "Relationship", "Plots", "Lag"),
                     width = c(1.3, 1.4, 1.0, 1.1, 1.8)) %>%
    flextable::fix_border_issues()
  
  doc <- officer::read_docx() %>% flextable::body_add_flextable(ft)
  print(doc, target = file.path(out_dir, "table_s4_2_lag_linear_vs_loess.docx"))
  cat("Word table written to", out_dir, "\n")
} else {
  message("Install flextable and officer to write the Word table; the csv was written")
}