# =============================================================================
# 04_species_pools_table.R
#
# Pipeline: pipeline/01_data_prep/02_species_pool/code/ (after 02 and 03)
# Author: Agostina Torres
#
# Table S2.1 of the Supporting Information: species pools used in the
# analyses, the criteria defining each pool, the number of overstory and
# understory species, the analyses using each pool, the reason for each
# criterion and the resulting cautions for interpretation (all
# species, shared pool, frequent set used only in the sampling diagnostics,
# and restricted pool). Numbers are entered by hand from the pool
# tables written by 02 (and summarised in species_table_si.docx by 03) and must be
# updated if the pools change.
#
# Outputs: pipeline/01_data_prep/02_species_pool/outputs/table_S2_1_species_pools.csv
#          pipeline/01_data_prep/02_species_pool/outputs/table_S2_1_species_pools.docx
# =============================================================================

library(tidyverse)
library(here)

out_dir <- here("pipeline", "01_data_prep", "02_species_pool", "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tab <- tribble(
  ~`Species pool`, ~Criteria, ~Overstory, ~Understory, ~`Used in`, ~Rationale, ~`Cautions for interpretation`,
  
  "All",
  "Species recorded in plots sampled at both times, after cleaning and taxonomic resolution",
  "23", "234", "-", "-", "-",
  
  "Shared",
  "Species recorded in at least one plot at each time",
  "15", "122",
  "Climatic-affinity shifts (S4)",
  "Reduces the risk of interpreting differences in detection, sampling or identification between surveys as regional colonizations or extinctions",
  "Regional colonizations and extinctions are excluded; if species recorded only in the present were mostly warm- or dry-affiliated colonizers, climatic-affinity shifts are conservative",
  
  "Frequent",
  "Shared species recorded in more than five plots at each time",
  "12", "69",
  "Diagnostics of the sampling design only (S3.2); intermediate step towards the restricted pool",
  "Ensures enough occurrences to estimate species responses",
  "-",
  
  "Restricted",
  "Frequent species with complete trait data after genus-level imputation",
  "12", "58 (21 woody, 37 non-woody)",
  "HMSC (S5)",
  "HMSC requires a complete trait matrix",
  "Rare species and species with incomplete traits are not modelled (11 of 69 frequent understory species lacked complete traits)"
)

# The table note is given in the caption of Table S2.1
note <- NULL

write_csv(tab, file.path(out_dir, "table_S2_1_species_pools.csv"))

# Word table (requires flextable and officer)
if (requireNamespace("flextable", quietly = TRUE) && requireNamespace("officer", quietly = TRUE)) {
  # Minimal layout: horizontal rules only above and below the header and at
  # the bottom of the body, no vertical lines or shading
  rule_thick <- officer::fp_border(color = "black", width = 1)
  rule_thin  <- officer::fp_border(color = "black", width = 0.5)
  ft <- flextable::flextable(tab) %>%
    flextable::border_remove() %>%
    flextable::hline_top(border = rule_thick, part = "header") %>%
    flextable::hline_bottom(border = rule_thin, part = "header") %>%
    flextable::hline_bottom(border = rule_thick, part = "body") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::fontsize(size = 9, part = "all") %>%
    flextable::bold(part = "header") %>%
    flextable::bold(j = "Species pool", part = "body") %>%
    flextable::align(align = "left", part = "all") %>%
    flextable::align(j = c("Overstory", "Understory"), align = "center", part = "all") %>%
    flextable::valign(valign = "top", part = "body") %>%
    flextable::padding(padding.top = 3, padding.bottom = 3, part = "all") %>%
    flextable::width(j = c("Species pool", "Overstory", "Understory", "Used in"), width = c(0.9, 0.7, 1.0, 1.0)) %>%
    flextable::width(j = c("Criteria", "Rationale", "Cautions for interpretation"), width = c(1.8, 2.0, 2.2)) %>%
    flextable::italic(i = 3, part = "body")     # frequent set: intermediate step, not an analysis pool
  doc <- officer::read_docx() %>%
    officer::body_set_default_section(officer::prop_section(page_size = officer::page_size(orient = "landscape"))) %>%
    flextable::body_add_flextable(ft)
  print(doc, target = file.path(out_dir, "table_S2_1_species_pools.docx"))
  cat("Word table written to", out_dir, "\n")
} else {
  message("Install flextable and officer to write the Word table; the csv was written")
}

# Markdown version for copying into a text
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(tab, format = "pipe"))
}
