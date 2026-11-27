# ------------------------------------------------------------------------------
# 04_check_try_gift_integration.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: check the TRY + GIFT integration (02_build_trait_table_try.R +
# 03_build_trait_table_gift.R) before genus-mean/mode imputation
# (05_impute_trait_by_genus.R). Three checks:
#
#   1. Units: for SLA and seed_mass, the two numeric traits from GIFT, each
#      species' TRY value is compared with its GIFT value after the unit
#      conversion (SLA x0.1, seed_mass x1000; see 03_build_trait_table_gift.R).
#      The ratio is not expected to be near 1 for every species (two
#      independent sources), but it should not be systematically off by
#      ~10x/0.1x or ~1000x/0.001x, which would indicate a wrong conversion
#      factor rather than disagreement between sources.
#
#   2. Source breakdown: for every species x trait, the provenance tag
#      (<trait>_src in trait_species_resolved_try_gift, from
#      03_build_trait_table_gift.R) is grouped into TRY only
#      ("direct"/"synonym:..."), GIFT only ("gift_only"), averaged
#      ("try+gift_avg", or "try+gift:agree"/"...+gift:conflict" for
#      growth_form), or unresolved. There is no "genus imputed" group at
#      this stage.
#
#   3. Coverage and Callitropsis: coverage per trait x layer at this stage
#      (TRY + GIFT only, before genus or phylogenetic imputation and the
#      literature values), and a check of Callitropsis nootkatensis. Its SLA
#      is unresolved here and remains unresolved after
#      05_impute_trait_by_genus.R, since it is the only project species in
#      its genus; its SLA is therefore taken from the literature in
#      08_integration.R (Mitchell 2001 via Brinkman & Associates 2012,
#      10.0 cm2/g -> 1.0 mm2/mg).
#
# Run AFTER 03_build_trait_table_gift.R, BEFORE 05_impute_trait_by_genus.R.
#
# Output:
#   - tmp_files/unit_check_SLA.csv, tmp_files/unit_check_seed_mass.csv
#   - tmp_files/trait_source_breakdown_try_gift_overstory.csv / _understory.csv
#     (species x trait x source_bucket, full detail, pre-genus-imputation)
#   - tmp_files/trait_source_summary_try_gift.csv (counts per trait x layer
#     x source_bucket)
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

obj <- readRDS(this_module_dir("tmp_files", "trait_matching_intermediate.rds"))
list2env(obj, envir = environment())  # overstory_species, understory_species, trait_species_resolved_try_gift, etc.

# ==============================================================================
# 1. UNIT SANITY CHECK -- SLA and seed_mass, the two numeric traits GIFT
#    contributes, before trusting their try+gift_avg values.
# ==============================================================================
gift_matched <- read_csv(this_module_dir("tmp_files", "gift_species_level.csv"), show_col_types = FALSE)
try_resolved <- read_csv(this_module_dir("tmp_files", "try_species_trait_summary_resolved.csv"), show_col_types = FALSE)

unit_conversion <- c(SLA = 0.1, seed_mass = 1000)   # same factors as 03_build_trait_table_gift.R

check_units <- function(trait_label) {
  value_col <- paste0(trait_label, "_value")
  
  gift_tr <- gift_matched %>%
    filter(trait == trait_label) %>%
    distinct(match_species, .keep_all = TRUE) %>%
    transmute(species = match_species, gift_raw = as.numeric(valor))
  
  try_tr <- try_resolved %>%
    select(species, try_value = all_of(value_col)) %>%
    filter(!is.na(try_value))
  
  both <- inner_join(try_tr, gift_tr, by = "species") %>%
    mutate(
      gift_converted = gift_raw * unit_conversion[[trait_label]],
      ratio_try_over_gift = try_value / gift_converted
    )
  
  cat("\n--", trait_label, "-- species with BOTH a TRY and a GIFT value:", nrow(both), "\n")
  if (nrow(both) > 0) {
    cat("  median ratio (TRY / converted GIFT):", round(median(both$ratio_try_over_gift, na.rm = TRUE), 2),
        "-- should be in the right ballpark of 1 if the conversion factor is correct;",
        "a median far from 1 (e.g. ~10x/0.1x or ~1000x/0.001x) points at a wrong",
        "conversion factor, not just species-level disagreement.\n")
    print(both %>% arrange(desc(abs(log(ratio_try_over_gift)))), n = Inf)
  }
  
  write_csv(both, this_module_dir("tmp_files", paste0("unit_check_", trait_label, ".csv")))
  both
}

unit_check_sla       <- check_units("SLA")
unit_check_seed_mass <- check_units("seed_mass")

# ==============================================================================
# 2 + 3. SOURCE BREAKDOWN + COVERAGE, from trait_species_resolved_try_gift
#        (pre-genus-imputation) restricted to each layer's species list.
# ==============================================================================
classify_source <- function(src) {
  case_when(
    is.na(src)                                   ~ "unresolved",
    src == "unresolved"                           ~ "unresolved",
    src %in% c("direct") | str_starts(src, "synonym:")       ~ "TRY only",
    src == "gift_only"                            ~ "GIFT only",
    src %in% c("try+gift_avg", "try+gift:agree") |
      str_detect(src, "\\+gift:conflict$")        ~ "TRY+GIFT (averaged / agreement)",
    TRUE                                           ~ "other (check)"
  )
}

source_breakdown <- function(species_list, stratum_label) {
  src_cols <- names(trait_species_resolved_try_gift) %>% str_subset("_src$")
  species_list %>%
    select(species) %>%
    left_join(trait_species_resolved_try_gift %>% select(species, all_of(src_cols)), by = "species") %>%
    rename_with(~ str_remove(.x, "_src$"), all_of(src_cols)) %>%
    pivot_longer(-species, names_to = "trait", values_to = "src") %>%
    mutate(source_bucket = classify_source(src), stratum = stratum_label)
}

overstory_breakdown  <- source_breakdown(overstory_species,  "overstory")
understory_breakdown <- source_breakdown(understory_species, "understory")

write_csv(overstory_breakdown,  this_module_dir("tmp_files", "trait_source_breakdown_try_gift_overstory.csv"))
write_csv(understory_breakdown, this_module_dir("tmp_files", "trait_source_breakdown_try_gift_understory.csv"))

trait_source_summary <- bind_rows(overstory_breakdown, understory_breakdown) %>%
  count(stratum, trait, source_bucket) %>%
  arrange(trait, stratum, desc(n))

write_csv(trait_source_summary, this_module_dir("tmp_files", "trait_source_summary_try_gift.csv"))

cat("\n========================================\n")
cat("TRAIT SOURCE BREAKDOWN (per trait x stratum, TRY + GIFT only, before genus imputation):\n")
cat("========================================\n")
print(trait_source_summary, n = Inf)

cat("\nCoverage (% resolved by TRY/GIFT alone, before genus/phylogenetic imputation):\n")
coverage_now <- trait_source_summary %>%
  group_by(stratum, trait) %>%
  summarise(n_total = sum(n), n_resolved = sum(n[source_bucket != "unresolved"]),
            pct_resolved = round(100 * n_resolved / n_total, 1), .groups = "drop") %>%
  arrange(trait, stratum)
print(coverage_now, n = Inf)

# ==============================================================================
# 3b. Callitropsis nootkatensis -- explicit check, ahead of genus imputation
#     (05) and the literature override (08_integration.R).
# ==============================================================================
cat("\n========================================\n")
cat("Callitropsis nootkatensis -- current status, all traits (overstory):\n")
cat("========================================\n")
callitropsis_status <- overstory_breakdown %>% filter(species == "Callitropsis nootkatensis")
print(callitropsis_status, n = Inf)

callitropsis_sla_src <- callitropsis_status %>% filter(trait == "SLA") %>% pull(source_bucket)
if (length(callitropsis_sla_src) == 0) {
  cat("\nWARNING: Callitropsis nootkatensis not found in the overstory species list at all --",
      "check it's still in species_pool_broad.csv before relying on 10_integration.R to",
      "patch its SLA value.\n")
} else if (callitropsis_sla_src == "unresolved") {
  cat("\nSLA is 'unresolved' at this stage, as expected -- Callitropsis is the sole project",
      "species in its genus, so 05_impute_trait_by_genus.R's genus-mean imputation will have",
      "nothing to draw on for it either. It should still be 'unresolved' after 05, which is",
      "exactly when 10_integration.R's literature override (Mitchell 2001, 10.0 cm2/g ->",
      "1.0 mm2/mg) is designed to apply. Re-check after running 05 to confirm.\n")
} else {
  cat("\nNOTE: SLA is already resolved here via '", callitropsis_sla_src, "' (TRY and/or GIFT",
      "already cover it) -- 10_integration.R's literature override will NOT replace this value",
      "(it only fills rows still 'unresolved'). If the literature value should apply",
      "regardless, that needs an explicit decision -- flagging rather than silently skipping it.\n", sep = "")
}
