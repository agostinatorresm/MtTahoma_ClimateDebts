# ------------------------------------------------------------------------------
# 00_check_cover_definition.R
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: determine, empirically, what kind of cover value is stored in
# understory_long.csv's "cover" column -- needed to write the manuscript's
# methods text accurately (absolute per-species % cover, independently
# estimated per species and allowed to exceed 100% in total due to
# overlapping vegetation layers, vs. relative cover, where all species'
# values sum to 100% within a plot). Diagnostic only; nothing written to disk.
#
# LOGIC: for each plot x period, sum cover across all species present.
#   - If totals cluster near 100 (e.g. mostly 95-105), cover is relative/
#     compositional (species shares of one total).
#   - If totals are frequently well above 100 (and vary a lot plot to plot),
#     cover is absolute/independent per-species estimates (the normal case
#     for field cover estimates across overlapping strata).
# A secondary check also looks at whether any single species exceeds 100%
# cover on its own -- only possible for absolute cover, and if it occurs,
# conclusively rules out the relative/compositional interpretation.
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

cleaning_dir <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)

understory_long <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE)

cat("Columns in understory_long.csv:\n")
print(names(understory_long))

cat("\nRange of the raw 'cover' column (single species x plot x period values):\n")
print(summary(understory_long$cover))

# 1. Any single species already over 100% on its own? -----------------------
over_100_single <- understory_long %>% filter(cover > 100)
cat("\nSingle species x plot x period records with cover > 100:", nrow(over_100_single), "\n")
if (nrow(over_100_single) > 0) {
  cat("(if > 0, cover cannot be a relative/compositional measure -- a single\n",
      " species can't exceed 100% of a total that sums to 100%)\n")
  print(head(over_100_single, 10))
}

# 2. Total cover per plot x period (summed across all species present) ------
total_cover_per_plot <- understory_long %>%
  filter(occurrence == 1) %>%
  group_by(plot_id, period) %>%
  summarise(total_cover = sum(cover, na.rm = TRUE), n_species = n(), .groups = "drop")

cat("\nTotal cover summed across species, per plot x period (n =", nrow(total_cover_per_plot), "plot-periods):\n")
print(summary(total_cover_per_plot$total_cover))

cat("\n% of plot-periods with total cover in different ranges:\n")
total_cover_per_plot %>%
  mutate(range = cut(total_cover,
                     breaks = c(-Inf, 50, 90, 110, 150, 200, Inf),
                     labels = c("<50", "50-90", "90-110 (~100)", "110-150", "150-200", ">200"))) %>%
  count(range) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

cat("\n========================================\n")
cat("INTERPRETATION\n")
cat("========================================\n")
cat("If most plot-periods fall in the ~90-110 bucket and no single species\n",
    "ever exceeds 100%, cover is very likely RELATIVE (species shares summing\n",
    "to ~100% per plot).\n",
    "If totals are widely spread, frequently well above 100, and/or\n",
    "individual species exceed 100%, cover is very likely ABSOLUTE\n",
    "(independent per-species % estimates, not constrained to sum to 100 --\n",
    "the standard convention when cover is assessed per vertical stratum, as\n",
    "the four 1-m2 subplots were).\n")
