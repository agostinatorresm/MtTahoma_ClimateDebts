# ------------------------------------------------------------------------------
# 08_integration.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Last step of the trait pipeline: runs after 05 (genus imputation) and
# before 09 (validation) and 10 (figures). Two steps:
#
#   1. Phylogenetic imputation is not used: any value tagged "phylogenetic*"
#      (07_impute_phylogenetic.R) is set back to "unresolved". Leave-one-out
#      validation showed genus imputation predicts better, and a value that
#      only phylogenetic imputation could fill has no genus-level evidence.
#   2. Literature values: species x trait cells still "unresolved" are filled
#      with values from published sources (numeric: literature_overrides_raw;
#      categorical: literature_overrides_categorical). Literature values never
#      replace a TRY, GIFT, BIEN or genus value.
#
# Literature values are kept in this separate step so they can be added or
# edited without rerunning the database scripts (02-05). Units are converted
# to the TRY scale on load (SLA cm2/g -> mm2/mg, x0.1; seed mass g -> mg,
# x1000).
#
# Provenance: filled values are tagged "literature:<citation>" with
# _n_records = NA.
#
# Inputs : outputs/overstory_traits.csv, outputs/understory_traits.csv (05)
# Outputs: the same files, overwritten with the patched tables; these are the
#          trait tables read by the species pool and HMSC scripts. Because they
#          are overwritten, 05 must be rerun before rerunning this script.
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

overstory_traits  <- read_csv(this_module_dir("outputs", "overstory_traits.csv"),  show_col_types = FALSE)
understory_traits <- read_csv(this_module_dir("outputs", "understory_traits.csv"), show_col_types = FALSE)

# ==============================================================================
# Step 1: revert phylogenetic imputation (see header). Applied to every trait
# column that has one, not only the traits shown in the plotting scripts,
# since this is a change of the data, not a display choice.
# ==============================================================================
trait_cols_to_revert <- c("seed_mass", "longevity", "SLA", "wood_density", "growth_form", "woodiness")

revert_phylogenetic_imputation <- function(traits_df, trait_label, stratum_label) {
  if (!trait_label %in% names(traits_df)) return(traits_df)
  
  src_col <- paste0(trait_label, "_src")
  n_col   <- paste0(trait_label, "_n_records")
  
  is_phylo <- str_starts(traits_df[[src_col]], "phylogenetic")
  is_phylo[is.na(is_phylo)] <- FALSE
  
  if (any(is_phylo)) {
    cat(" ", stratum_label, "-", trait_label, ": reverting", sum(is_phylo),
        "phylogenetically-imputed value(s) to unresolved -- candidates for",
        "eventual drop from the species pool (", paste(traits_df$species[is_phylo], collapse = ", "), ")\n")
  }
  
  traits_df[[trait_label]] <- ifelse(is_phylo, NA, traits_df[[trait_label]])
  traits_df[[src_col]]     <- ifelse(is_phylo, "unresolved", traits_df[[src_col]])
  traits_df[[n_col]]       <- ifelse(is_phylo, NA_integer_, traits_df[[n_col]])
  traits_df
}

cat("Step 1 -- reverting phylogenetic imputation to unresolved:\n")
cat("Overstory:\n")
for (tr in trait_cols_to_revert) overstory_traits <- revert_phylogenetic_imputation(overstory_traits, tr, "overstory")
cat("Understory:\n")
for (tr in trait_cols_to_revert) understory_traits <- revert_phylogenetic_imputation(understory_traits, tr, "understory")

# ==============================================================================
# Step 2: literature values (one row per species x trait)
# ==============================================================================
literature_overrides_raw <- tribble(
  ~species,                    ~trait, ~raw_value, ~raw_unit,
  "Callitropsis nootkatensis", "SLA",  10.0,       "cm2/g"
) %>%
  mutate(source = case_when(
    species == "Callitropsis nootkatensis" & trait == "SLA" ~
      "Mitchell 2001 (avg of sunlight/shaded values), cited in Brinkman & Associates 2012, FFESC Chapter 6 Appendices, Table A6.1 (https://ffesc.brinkman.ca/wp-content/uploads/2012/03/Chapter-6-Appendices.pdf)"
  ))

# Categorical literature values (woodiness, growth form): species-level
# descriptions from floras, used where TRY, GIFT and congeners give no value
literature_overrides_categorical <- tribble(
  ~species,                      ~trait,      ~value,      ~source,
  "Trautvetteria caroliniensis", "woodiness", "non-woody", "Flora of North America (perennial herb)",
  "Lysichitum americanum",       "woodiness", "non-woody", "Flora of North America (perennial herb)",
  "Anemonoides grayi",           "woodiness", "non-woody", "Flora of North America (perennial herb)"
)

sla_unit_conversion       <- 0.1   # cm2/g -> mm2/mg
seed_mass_unit_conversion <- 1000  # g -> mg

literature_overrides <- literature_overrides_raw %>%
  mutate(value = case_when(
    trait == "SLA"       & raw_unit == "cm2/g" ~ raw_value * sla_unit_conversion,
    trait == "seed_mass" & raw_unit == "g"      ~ raw_value * seed_mass_unit_conversion,
    TRUE ~ raw_value
  ))

# Applied to whichever of the two trait tables has that trait column
# (wood_density does not exist in understory_traits, so the join finds
# nothing to fill there).
apply_literature_override <- function(traits_df, trait_label) {
  if (!trait_label %in% names(traits_df)) return(traits_df)
  
  src_col <- paste0(trait_label, "_src")
  n_col   <- paste0(trait_label, "_n_records")
  
  ov <- bind_rows(
    literature_overrides %>% transmute(species, trait, value = as.character(value), source),
    literature_overrides_categorical
  ) %>% filter(trait == trait_label)
  if (is.numeric(traits_df[[trait_label]])) ov <- ov %>% mutate(value = as.numeric(value))
  if (nrow(ov) == 0) return(traits_df)
  
  fill <- traits_df %>%
    select(species, !!src_col) %>%
    filter(.data[[src_col]] == "unresolved") %>%
    inner_join(ov %>% select(species, lit_value = value, source), by = "species")
  
  if (nrow(fill) > 0) {
    cat(" ", trait_label, ": literature override applied for", nrow(fill), "species (",
        paste(fill$species, collapse = ", "), ")\n")
  }
  
  traits_df %>%
    left_join(fill %>% select(species, lit_value, source), by = "species") %>%
    mutate(
      !!src_col       := if_else(!is.na(lit_value), paste0("literature:", source), .data[[src_col]]),
      !!n_col         := if_else(!is.na(lit_value), NA_integer_, .data[[n_col]]),
      !!trait_label   := coalesce(.data[[trait_label]], lit_value)
    ) %>%
    select(-lit_value, -source)
}

trait_labels <- unique(c(literature_overrides$trait, literature_overrides_categorical$trait))

cat("Step 2 -- literature overrides -- overstory_traits.csv:\n")
for (tr in trait_labels) overstory_traits <- apply_literature_override(overstory_traits, tr)

cat("\nStep 2 -- literature overrides -- understory_traits.csv:\n")
for (tr in trait_labels) understory_traits <- apply_literature_override(understory_traits, tr)

write_csv(overstory_traits,  this_module_dir("outputs", "overstory_traits.csv"))
write_csv(understory_traits, this_module_dir("outputs", "understory_traits.csv"))

# Still-unresolved species x trait rows after both steps -- the pool to
# review for eventual dropping from the species pools (not done here).
still_unresolved <- bind_rows(
  overstory_traits  %>% select(species, ends_with("_src")) %>% mutate(stratum = "overstory"),
  understory_traits %>% select(species, ends_with("_src")) %>% mutate(stratum = "understory")
) %>%
  pivot_longer(cols = ends_with("_src"), names_to = "trait_src", values_to = "src") %>%
  mutate(trait = str_remove(trait_src, "_src$")) %>%
  filter(src == "unresolved", trait %in% trait_cols_to_revert) %>%
  select(species, stratum, trait)

cat("\n========================================\n")
cat("Trait tables patched (phylogenetic imputation reverted, literature overrides applied) and re-written:\n")
cat(" outputs/overstory_traits.csv\n")
cat(" outputs/understory_traits.csv\n")
cat("\nStill unresolved after both steps (", nrow(still_unresolved), "species x trait rows -- review for eventual drop):\n")
print(still_unresolved, n = Inf)
cat("========================================\n")
