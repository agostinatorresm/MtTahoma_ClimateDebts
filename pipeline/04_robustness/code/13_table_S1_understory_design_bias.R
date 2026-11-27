# =============================================================================
# 13_table_S1_understory_design_bias.R
#
# Pipeline: pipeline/04_robustness/code/
# Author: Agostina Torres
#
# Table S1 of the Supporting Information: potential biases arising from the
# difference in understory sampling between the past and present surveys, the
# corrections applied to the data, the analyses in which each correction was
# applied, the sensitivity analyses run and the resulting cautions for
# interpretation. One row per source of bias; references to the Supporting
# Information sections are given in each cell. Values are entered by hand from
# the results of the robustness and HMSC scripts and must be updated when those
# results change.
#
# Outputs: pipeline/04_robustness/outputs/table_S1_understory_design_bias.csv
#          pipeline/04_robustness/outputs/table_S1_understory_design_bias.docx
# =============================================================================

library(tidyverse)
library(here)

out_dir <- here("pipeline", "04_robustness", "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tab <- tribble(
  ~`Source of potential bias`, ~Description, ~`Data correction`, ~`Applied in`, ~`Sensitivity analyses`, ~`Cautions for interpretation`,
  
  "(1) Detection of sparse species occurrences",
  "Past occurrence was recorded over the whole plot (250-1000 m2) and present occurrence over a 50-m2 subplot plus four 1-m2 quadrats, so more species, mainly low-cover ones, were recorded per plot in the past (S3.2)",
  "Past data reduced to the present sampling effort, separately for each pool (eight species removed per plot in the shared pool, seven in the restricted pool), with probability inversely related to cover, 100 replicates (Wilson et al., 2024; S3.3)",
  "HMSC (main occurrence model, 10 reduced datasets pooled); climatic-affinity shifts (sensitivity only)",
  "Cover-weighted vs uniform removal (S3.3); observed vs reduced data (CTI and CMI, S4.1; HMSC, S5.3); 10 vs 100 reductions (single-species GLMs, S5.3); without 12 sparse species (S5.3)",
  "Decline shared by all species is removed by the reduction and not interpreted; occurrence-based climatic indices depend on sampling effort and are not interpreted (S4.1)",
  
  "(2) Detection of patchy species cover",
  "Present cover was estimated only in four 1-m2 quadrats, so species with low or patchy cover often had no present cover value or a less precise one, and fewer species had cover per plot, inflating relative cover (S3.4)",
  "HMSC: cover paired by plot (S3.5) and species changes centred on the mean change across species (S5.4). Climatic-affinity shifts: none, because method differences that affect all species alike cancel out in a weighted mean",
  "HMSC (relative cover model)",
  "Paired vs all presences with cover (S5.4); without the 18 patchy species (S5.4)",
  "Overall change in relative cover not interpreted (S5.4); seven species with one or two paired plots not interpreted individually; Moneses uniflora modelled for occurrence only"
)

write_csv(tab, file.path(out_dir, "table_S1_understory_design_bias.csv"))

# Word table (requires flextable and officer); minimal layout: horizontal rules
# only above and below the header and at the bottom of the body
if (requireNamespace("flextable", quietly = TRUE) && requireNamespace("officer", quietly = TRUE)) {
  rule_thick <- officer::fp_border(color = "black", width = 1)
  rule_thin  <- officer::fp_border(color = "black", width = 0.5)
  ft <- flextable::flextable(tab) %>%
    flextable::border_remove() %>%
    flextable::hline_top(border = rule_thick, part = "header") %>%
    flextable::hline_bottom(border = rule_thin, part = "header") %>%
    flextable::hline(i = 1, border = rule_thin, part = "body") %>%
    flextable::hline_bottom(border = rule_thick, part = "body") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 8.5, part = "all") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = "Source of potential bias", part = "body") %>%
    flextable::align(align = "left", part = "all") %>%
    flextable::valign(valign = "top", part = "body") %>%
    flextable::padding(padding.top = 3, padding.bottom = 3, part = "all") %>%
    flextable::width(j = c("Source of potential bias", "Applied in"), width = c(1.1, 1.3)) %>%
    flextable::width(j = c("Description", "Data correction", "Sensitivity analyses", "Cautions for interpretation"),
                     width = c(2.0, 2.0, 2.1, 1.7))
  doc <- officer::read_docx() %>%
    officer::body_set_default_section(officer::prop_section(page_size = officer::page_size(orient = "landscape"))) %>%
    flextable::body_add_flextable(ft)
  print(doc, target = file.path(out_dir, "table_S1_understory_design_bias.docx"))
  cat("Word table written to", out_dir, "\n")
} else {
  message("Install flextable and officer to write the Word table; the csv was written")
}

# Markdown version for copying into a text
if (requireNamespace("knitr", quietly = TRUE)) print(knitr::kable(tab, format = "pipe"))
