# ------------------------------------------------------------------------------
# 07_impute_phylogenetic.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: last imputation step. Fills species x trait cells still
# "unresolved" after 05_impute_trait_by_genus.R (no congener with TRY/GIFT
# data in the genus) using the project phylogeny (06_build_phylogeny.R),
# borrowing information from more distant relatives through ancestral-state
# estimation.
#
# METHOD (castor package, "hidden state prediction"):
#   - Numeric traits (seed_mass, longevity, SLA, wood_density): squared-change
#     parsimony (hsp_squared_change_parsimony), a fast method equivalent to
#     maximum likelihood under Brownian motion for estimating a tip's
#     continuous trait value from its relatives' known values.
#   - Categorical traits (growth_form, woodiness): an equal-rates Markov
#     (Mk) model of character evolution (hsp_mk_model, rate_model = "ER").
#     The predicted category is the state with the highest posterior
#     probability at that tip; that probability is kept as a confidence
#     value.
#
# INFORMATIVE TIPS: only species-level TRY/GIFT values are given to the
# phylogenetic model as known tips, i.e. every <trait>_src value except
# "genus_mean:*" / "genus_mode:*" (already a within-genus average, not
# independent evidence) and "unresolved". Using genus means as observations
# would overstate what is known about the more distantly related,
# unresolved species. wood_density is overstory-only, as in the previous
# scripts (the column does not exist in understory_traits.csv).
#
# Value hierarchy per species x trait: direct/synonym/gift >
# genus_mean/genus_mode (05) > phylogenetic (this script) > unresolved.
# A row's <trait>_src is changed here only when it is "unresolved" and the
# species' tip is on the phylogeny; all other rows (including those filled
# with a genus mean in 05) are unchanged, so 08_integration.R can still
# complete any row that remains "unresolved".
#
# Output: overwrites outputs/overstory_traits.csv and
# outputs/understory_traits.csv in place, adding a phylogenetic-imputed
# value + "phylogenetic" src tag wherever the phylogeny could supply one,
# plus a <trait>_phylo_confidence column for categorical traits (NA for
# everything that isn't a phylogenetic-imputed categorical cell).
#   - outputs/phylo_imputation_summary.csv: coverage per trait x layer,
#     broken down by src category (observed / genus / phylogenetic /
#     still_unresolved), after this script runs.
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(here)
library(ape)

if (!requireNamespace("castor", quietly = TRUE)) {
  stop("castor is not installed -- install.packages('castor') and re-run this script.")
}
library(castor)

# 0a. Module directories --------------------------------------------------------
this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

NUMERIC_TRAITS     <- c("seed_mass", "longevity", "SLA", "wood_density")
CATEGORICAL_TRAITS <- c("growth_form", "woodiness")
ALL_TRAITS         <- c(NUMERIC_TRAITS, CATEGORICAL_TRAITS)

# 1. Load phylogeny and the two trait tables -------------------------------------
# phy$tip.label is already "Genus species" (space-separated) -- 06_build_
# phylogeny.R restores project accepted names in that form as its very last
# step, so no underscore <-> space conversion is needed here.
load(this_module_dir("outputs", "MORA_combined_phylo_tree.rda"))  # loads C, phy
n_tips <- length(phy$tip.label)
cat("Tree tips:", n_tips, "\n")

overstory_traits  <- read_csv(this_module_dir("outputs", "overstory_traits.csv"),  show_col_types = FALSE) %>%
  mutate(stratum = "overstory")
understory_traits <- read_csv(this_module_dir("outputs", "understory_traits.csv"), show_col_types = FALSE) %>%
  mutate(stratum = "understory")

# wood_density is overstory-only (no such column in understory_traits.csv) --
# every function below checks for column existence per trait rather than
# assuming both tables have the same trait columns.

# 2. src-category classifier ----------------------------------------------------
# Groups the provenance tags written by 02/03 (direct, synonym:*, gift_only,
# try+gift_avg, try+gift:agree, <src>+gift:conflict), 05 (genus_mean:*,
# genus_mode:*, unresolved) and this script (phylogenetic) into four groups.
#
# "phylogenetic" is a separate group, not part of "observed", although the
# input from 05 has no "phylogenetic" rows. If this script is run twice on
# the same overstory/understory_traits.csv without re-running 05 in
# between, values imputed in the first run would otherwise be used as
# independent species-level evidence in the second.
classify_src <- function(src) {
  case_when(
    is.na(src) ~ NA_character_,
    src == "unresolved" ~ "unresolved",
    src == "phylogenetic" ~ "phylogenetic",
    str_starts(src, "genus_mean:") | str_starts(src, "genus_mode:") ~ "genus",
    TRUE ~ "observed"
  )
}

# 3. Numeric trait: hidden-state prediction by squared-change parsimony ---------
impute_numeric_phylo <- function(trait, traits_tbl) {
  src_col <- paste0(trait, "_src")
  if (!trait %in% names(traits_tbl)) return(NULL)
  
  src_class <- classify_src(traits_tbl[[src_col]])
  observed <- traits_tbl %>%
    mutate(src_class = src_class) %>%
    filter(src_class == "observed", !is.na(.data[[trait]])) %>%
    select(species, value = all_of(trait))
  
  tip_states <- rep(NA_real_, n_tips)
  idx <- match(observed$species, phy$tip.label)
  keep <- !is.na(idx)
  tip_states[idx[keep]] <- observed$value[keep]
  
  cat(" ", trait, "-- informative tips (species-level TRY/GIFT, on this tree):",
      sum(keep), "of", n_tips, "\n")
  
  if (sum(keep) < 2) {
    cat("   skipped (fewer than 2 informative tips) --", trait, "cannot be phylogenetically imputed.\n")
    return(NULL)
  }
  
  result <- hsp_squared_change_parsimony(phy, tip_states = tip_states, weighted = TRUE, check_input = TRUE)
  if (trait == NUMERIC_TRAITS[1]) {
    cat("   hsp_squared_change_parsimony() result fields:", paste(names(result), collapse = ", "), "\n")
  }
  
  tibble(
    species          = phy$tip.label,
    trait            = trait,
    phylo_value      = result$states[1:n_tips],
    phylo_confidence = NA_character_
  )
}

# 4. Categorical trait: hidden-state prediction under an equal-rates Mk model ---
impute_categorical_phylo <- function(trait, traits_tbl) {
  src_col <- paste0(trait, "_src")
  if (!trait %in% names(traits_tbl)) return(NULL)
  
  src_class <- classify_src(traits_tbl[[src_col]])
  observed <- traits_tbl %>%
    mutate(src_class = src_class) %>%
    filter(src_class == "observed", !is.na(.data[[trait]])) %>%
    select(species, value = all_of(trait))
  
  levels_vec <- sort(unique(observed$value))
  tip_states <- rep(NA_integer_, n_tips)
  idx <- match(observed$species, phy$tip.label)
  keep <- !is.na(idx)
  tip_states[idx[keep]] <- match(observed$value[keep], levels_vec)
  
  cat(" ", trait, "-- informative tips (species-level TRY/GIFT, on this tree):", sum(keep), "of", n_tips,
      " | categories:", paste(levels_vec, collapse = ", "), "\n")
  
  if (sum(keep) < 2 || length(levels_vec) < 2) {
    cat("   skipped (fewer than 2 informative tips or fewer than 2 categories) --",
        trait, "cannot be phylogenetically imputed.\n")
    return(NULL)
  }
  
  result <- hsp_mk_model(phy, tip_states = tip_states, rate_model = "ER",
                         Nstates = length(levels_vec), check_input = TRUE)
  if (trait == CATEGORICAL_TRAITS[1]) {
    cat("   hsp_mk_model() result fields:", paste(names(result), collapse = ", "), "\n")
  }
  
  tip_likelihoods <- result$likelihoods[1:n_tips, , drop = FALSE]
  predicted_idx   <- apply(tip_likelihoods, 1, which.max)
  predicted_conf  <- apply(tip_likelihoods, 1, max)
  
  tibble(
    species          = phy$tip.label,
    trait            = trait,
    phylo_value      = levels_vec[predicted_idx],
    phylo_confidence = as.character(round(predicted_conf, 3))
  )
}

# 5. Apply phylogenetic imputation to one trait table (overstory or understory) -
apply_phylo_to_table <- function(traits_tbl, stratum_label) {
  cat("\n--- ", stratum_label, " ---\n")
  
  for (trait in ALL_TRAITS) {
    if (!trait %in% names(traits_tbl)) next
    
    phylo_result <- if (trait %in% NUMERIC_TRAITS) {
      impute_numeric_phylo(trait, traits_tbl)
    } else {
      impute_categorical_phylo(trait, traits_tbl)
    }
    if (is.null(phylo_result)) next
    
    src_col   <- paste0(trait, "_src")
    n_col     <- paste0(trait, "_n_records")
    conf_col  <- paste0(trait, "_phylo_confidence")
    
    src_class <- classify_src(traits_tbl[[src_col]])
    
    fill <- traits_tbl %>%
      mutate(src_class = src_class) %>%
      select(species, src_class) %>%
      filter(src_class == "unresolved") %>%
      inner_join(phylo_result, by = "species")
    
    cat("  ", trait, ": phylogenetic value available for", nrow(fill),
        "of the", sum(src_class == "unresolved", na.rm = TRUE), "still-unresolved species\n")
    
    traits_tbl <- traits_tbl %>%
      left_join(fill %>% select(species, phylo_value, phylo_confidence), by = "species") %>%
      mutate(
        was_filled := !is.na(phylo_value) & .data[[src_col]] == "unresolved",
        !!trait    := if_else(was_filled, phylo_value, .data[[trait]]),
        !!src_col  := if_else(was_filled, "phylogenetic", .data[[src_col]]),
        !!n_col    := if_else(was_filled, NA_integer_, .data[[n_col]]),
        !!conf_col := if_else(was_filled, phylo_confidence, NA_character_)
      ) %>%
      select(-phylo_value, -phylo_confidence, -was_filled)
  }
  
  traits_tbl
}

cat("\nRunning phylogenetic hidden-state prediction per trait:\n")
overstory_traits  <- apply_phylo_to_table(overstory_traits,  "overstory")
understory_traits <- apply_phylo_to_table(understory_traits, "understory")

# 6. Save (overwrite in place, same convention as 08_integration.R) -------------
write_csv(overstory_traits  %>% select(-stratum), this_module_dir("outputs", "overstory_traits.csv"))
write_csv(understory_traits %>% select(-stratum), this_module_dir("outputs", "understory_traits.csv"))

# 7. Coverage summary per trait x layer x src category ------------------------
summarize_coverage <- function(traits_tbl, stratum_label) {
  map_dfr(ALL_TRAITS, function(trait) {
    src_col <- paste0(trait, "_src")
    if (!trait %in% names(traits_tbl)) return(NULL)
    
    traits_tbl %>%
      mutate(src_class = case_when(
        .data[[src_col]] == "unresolved" ~ "still_unresolved",
        .data[[src_col]] == "phylogenetic" ~ "phylogenetic",
        str_starts(.data[[src_col]], "genus_mean:") | str_starts(.data[[src_col]], "genus_mode:") ~ "genus",
        TRUE ~ "observed"
      )) %>%
      count(src_class) %>%
      pivot_wider(names_from = src_class, values_from = n, values_fill = 0) %>%
      mutate(trait = trait, stratum = stratum_label) %>%
      relocate(trait, stratum)
  })
}

phylo_summary <- bind_rows(
  summarize_coverage(overstory_traits,  "overstory"),
  summarize_coverage(understory_traits, "understory")
)

write_csv(phylo_summary, this_module_dir("outputs", "phylo_imputation_summary.csv"))

cat("\nFinal coverage (observed + genus + phylogenetic + still_unresolved), by trait x stratum:\n")
print(phylo_summary)

cat("\n========================================\n")
cat("PHYLOGENETIC IMPUTATION -- COMPLETE\n")
cat("========================================\n")
cat("Saved (overwritten in place):\n")
cat(" - outputs/overstory_traits.csv\n")
cat(" - outputs/understory_traits.csv\n")
cat(" - outputs/phylo_imputation_summary.csv\n")
cat("\nFor categorical traits, check <trait>_phylo_confidence for the rows this\n",
    "script filled -- a low top-state probability means the Mk model wasn't\n",
    "confident about that species even using its relatives.\n",
    "Any trait x species cell still 'unresolved' after this script has NO\n",
    "congener anywhere on the tree with species-level data -- 10_integration.R\n",
    "is the last chance to fill it, from a literature value.\n")
