# ------------------------------------------------------------------------------
# 03_build_trait_table_gift.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: extract SLA, seed mass and growth form from GIFT (Global
# Inventory of Floras and Traits) as a second trait source and combine it
# with the TRY table from 02_build_trait_table_try.R wherever both sources
# have data (not only as a fallback):
#
#   - TRY value (direct or synonym) and GIFT value: mean of the two (after
#     converting GIFT's value to TRY's units) -- "try+gift_avg"
#   - no TRY value, GIFT value: GIFT alone (unit-converted) -- "gift_only"
#   - TRY value, no GIFT value: TRY value unchanged
#   - neither: "unresolved", passed on to genus-mean/mode imputation
#     (05_impute_trait_by_genus.R)
#
# This applies to the 3 traits in GIFT (SLA, seed_mass, growth_form).
# Longevity, wood density and woodiness come from TRY only.
#
# growth_form is categorical:
#   - only one source has a value -> that value is used
#   - both sources agree -> kept ("try+gift:agree")
#   - the sources disagree -> the TRY value is kept (a per-species mode over
#     possibly many records, compared with a single aggregated GIFT value)
#     and the conflict is written out for manual review
#
# GIFT nomenclature follows WCVP, so species are matched with the same
# direct -> WFO-synonym sequence as the TRY matching
# (02_build_trait_table_try.R section 3b), rather than GIFT's own fuzzy
# matching, to keep species identity consistent across the pipeline.
#
# UNIT CONVERSION -- GIFT units differ from TRY StdValue units (to be
# checked against GIFT_traits_meta()$Units with any new GIFT version):
#   SLA:       GIFT cm2/g -> TRY mm2/mg  factor = 0.1   (cm2 = 100 mm2, g = 1000 mg)
#   seed_mass: GIFT g     -> TRY mg      factor = 1000  (TRY TraitID 26
#              StdValue is in mg, the TRY standard unit for seed dry mass)
#
# Output:
#   - tmp_files/gift_species_level.csv: GIFT trait values matched to project
#     species (direct + WFO synonym), same shape as TRY's resolved table
#   - tmp_files/gift_vs_try_genus_coverage.csv: genera covered by GIFT but
#     not by TRY, per trait, flagged for whether a project species belongs
#     to them (informational)
#   - tmp_files/gift_try_growth_form_conflicts.csv: species where TRY and
#     GIFT growth_form disagree -- TRY's value is kept, these are flagged for
#     manual review
#   - tmp_files/gift_try_numeric_conflicts.csv: species x trait (SLA,
#     seed_mass) where the unit-converted GIFT value differs from TRY by
#     more than NUMERIC_CONFLICT_RATIO_FOLD-fold -- the TRY value is kept
#     (not averaged with GIFT) and the case is flagged for manual review
#   - tmp_files/try_gift_species_trait_summary_resolved.csv: full project
#     pool, TRY+GIFT-combined trait table (input to 05_impute_trait_by_genus.R)
#   - tmp_files/trait_matching_intermediate.rds: re-saved with
#     trait_species_resolved_try_gift added alongside everything 02 already
#     saved (raw TRY tables, TRY-only trait_species_resolved, species pools)
# ------------------------------------------------------------------------------

library(GIFT)
library(tidyverse)
library(here)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

dir_raw <- this_module_dir("data_raw", "GIFT")
dir.create(dir_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("tmp_files"), recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(this_module_dir("tmp_files", "trait_matching_intermediate.rds"))
list2env(obj, envir = environment())  # project_species, spp_synonyms, raw TRY tables, trait_species_resolved, etc.

# Defined locally rather than relying on 02's own global-environment variable
# of the same name (not part of the .rds) -- this script should give the same
# result whether it's run right after 02 in the same R session or started
# fresh from the .rds alone.
trait_labels <- c("seed_mass", "longevity", "SLA", "wood_density", "growth_form", "woodiness")

# WFO synonym table plus the alternative names added in 02 (extra_synonyms,
# carried in spp_synonyms from the .rds: other genus treatments and spelling
# variants of the epithet)
spp_synonyms_all <- read_csv(this_module_dir("outputs", "spp_synonyms.csv"), show_col_types = FALSE) %>%
  select(final_species, synonym) %>%
  bind_rows(spp_synonyms %>% select(final_species, synonym)) %>%
  distinct(final_species, synonym)

# --- 1. Trait catalog match (SLA, seed mass, growth form) ----------------------
gift_meta <- GIFT_traits_meta() %>% as_tibble()

# All GIFT trait IDs whose Trait2 contains "SLA" (any case) are printed:
# GIFT lists several SLA-related variables (e.g. min/max/mean), and the
# list makes any renaming in a new GIFT version visible.
cat("All GIFT trait codes with 'SLA' in Trait2 (for visibility -- see which gets picked below):\n")
print(gift_meta %>% filter(str_detect(Trait2, regex("SLA", ignore_case = TRUE))) %>%
        select(Lvl3, Trait2, Units))

mapa <- tribble(
  ~trait,        ~patron,            ~tipo,
  "SLA",         "^SLA",             "num",
  "seed_mass",   "^Seed_mass_mean",  "num",
  "growth_form", "^Growth_form_1",   "cual"
)

# When a pattern matches more than one GIFT trait (as "^SLA" does), the one
# whose name contains "mean" is preferred, so that a per-species mean is
# used rather than, e.g., SLA_min, which has few records; otherwise the
# first match is used.
# Written in base R rather than with dplyr verbs, to avoid an error in
# exists() raised inside mutate(map()) when library(GIFT) is loaded
# together with the tidyverse.
pick_gift_trait <- function(pattern) {
  idx <- which(grepl(pattern, gift_meta$Trait2, ignore.case = TRUE))
  if (length(idx) == 0) return(gift_meta[0, , drop = FALSE])
  if (length(idx) > 1) {
    mean_idx <- idx[grepl("mean", gift_meta$Trait2[idx], ignore.case = TRUE)]
    if (length(mean_idx) > 0) idx <- mean_idx
  }
  gift_meta[idx[1], , drop = FALSE]
}

ids_rows <- lapply(seq_len(nrow(mapa)), function(i) {
  found <- pick_gift_trait(mapa$patron[i])
  if (nrow(found) == 0) {
    warning("No GIFT trait matched pattern '", mapa$patron[i], "' for trait '",
            mapa$trait[i], "' -- check gift_meta$Trait2 naming (GIFT version change?).")
    return(NULL)
  }
  tibble(trait = mapa$trait[i], tipo = mapa$tipo[i],
         Lvl3 = found$Lvl3, Trait2 = found$Trait2, Units = found$Units)
})
ids <- bind_rows(ids_rows)

cat("\nGIFT trait codes matched (after mean-preference):\n"); print(ids)

# --- 2. Download (cached) -------------------------------------------------------
# The cache is keyed on the trait IDs requested (stored with the data), so
# a change in the selected traits forces a new download instead of reusing
# a cache built for other IDs.
path_cache    <- file.path(dir_raw, "GIFT_traits.rds")
path_cache_ids <- file.path(dir_raw, "GIFT_traits_ids.rds")

cached_ids_match <- file.exists(path_cache) && file.exists(path_cache_ids) &&
  identical(sort(readRDS(path_cache_ids)), sort(ids$Lvl3))

if (!cached_ids_match) {
  if (file.exists(path_cache)) {
    cat("Trait IDs requested have changed since the cache was built -- re-downloading GIFT_traits.\n")
  }
  GIFT_traits(trait_IDs = ids$Lvl3, agreement = 0.66, bias_ref = FALSE, bias_deriv = FALSE) %>%
    saveRDS(path_cache)
  saveRDS(ids$Lvl3, path_cache_ids)
}
gift_raw <- readRDS(path_cache) %>% as_tibble()

gift_long <- gift_raw %>%
  select(work_species, starts_with("trait_value_")) %>%
  pivot_longer(-work_species, names_to = "Lvl3", values_to = "valor",
               names_prefix = "trait_value_", values_transform = as.character) %>%
  filter(!is.na(valor)) %>%
  left_join(ids %>% select(Lvl3, trait, tipo), by = "Lvl3") %>%
  rename(species = work_species)

# --- 3. Match to project species: direct, then WFO synonym fallback -----------
# Mirrors the TRY cascade (02_build_trait_table_try.R section 3b) exactly.
project_sp <- project_species %>% select(species)

direct_match <- gift_long %>%
  semi_join(project_sp, by = "species") %>%
  transmute(match_species = species, species, trait, tipo, valor, src = "direct")

gift_via_synonym <- project_sp %>%
  anti_join(direct_match, by = c("species" = "match_species")) %>%
  inner_join(spp_synonyms_all, by = c("species" = "final_species")) %>%
  inner_join(gift_long, by = c("synonym" = "species")) %>%
  transmute(match_species = species, species = synonym, trait, tipo, valor,
            src = paste0("synonym:", synonym))

gift_matched <- bind_rows(direct_match, gift_via_synonym)

write_csv(gift_matched, this_module_dir("tmp_files", "gift_species_level.csv"))
cat("\nGIFT species-level matches:", n_distinct(gift_matched$match_species), "of",
    nrow(project_sp), "project species (any of the 3 traits)\n")

# --- 4. Genus-level comparison: TRY vs GIFT (informational) --------------------
# Diagnostic only: shows how much additional genus coverage GIFT provides.
extract_genus <- function(x) word(x, 1)

try_genera <- list(
  SLA = unique(extract_genus(sla_raw$species)),
  seed_mass = unique(extract_genus(seed_mass_raw$species)),
  growth_form = unique(extract_genus(growth_form_raw$species))
)

gift_genus_coverage <- gift_long %>%
  mutate(genus = extract_genus(species)) %>%
  distinct(trait, genus)

genus_comparison <- gift_genus_coverage %>%
  rowwise() %>%
  mutate(in_try = genus %in% try_genera[[trait]]) %>%
  ungroup()

gift_only_genera <- genus_comparison %>% filter(!in_try)

project_genera <- project_species %>% mutate(genus = extract_genus(species)) %>% distinct(genus)
relevant_gain <- gift_only_genera %>% semi_join(project_genera, by = "genus")

write_csv(relevant_gain, this_module_dir("tmp_files", "gift_vs_try_genus_coverage.csv"))

cat("\nGenera covered by GIFT but not TRY, relevant to a project genus (informational only):\n")
print(count(relevant_gain, trait, name = "n_project_genera_gained"))

# --- 5. Combine TRY + GIFT per species x trait ---------------------------------
# Numeric (SLA, seed_mass): average when both present (unit-converted GIFT) --
# UNLESS the two disagree by more than NUMERIC_CONFLICT_RATIO_FOLD-fold (see
# below), in which case this is treated the same way as a categorical
# growth_form conflict: TRY is kept (it already reflects a per-species mean
# over potentially many records), GIFT's value is discarded for averaging
# purposes and logged for manual review, rather than silently pulled into
# the average.
# Categorical (growth_form): agreement kept as-is; disagreement keeps TRY,
# logged as a conflict for review; GIFT used alone only when TRY has nothing.
gift_unit_conversion <- c(SLA = 0.1, seed_mass = 1000)

# Numeric TRY-vs-GIFT disagreement threshold. For Pseudotsuga menziesii,
# TRY gives ~11 mg seed mass (28 records, consistent with the species'
# biology), whereas the GIFT value converts to ~510 mg (unit "g" in
# GIFT_traits_meta()$Units, so the x1000 conversion is correct), a ~46-fold
# difference most likely due to an error in the GIFT source data. Species x
# trait pairs where the converted GIFT value and the TRY value differ by
# more than this fold-change (either direction) are therefore not averaged:
# the TRY value is kept and the case is written to
# tmp_files/gift_try_numeric_conflicts.csv for manual review, as for
# growth_form conflicts below.
NUMERIC_CONFLICT_RATIO_FOLD <- 5

harmonize_gift_growth_form <- function(v) {
  v <- str_to_lower(str_trim(v))
  if_else(v %in% c("herb", "shrub", "tree"), v, NA_character_)  # "other" excluded, ambiguous
}

growth_form_conflicts <- tibble(species = character(), try_value = character(), gift_value = character())
numeric_conflicts <- tibble(species = character(), trait = character(),
                            try_value = double(), gift_value = double(), ratio = double())

combine_try_gift <- function(resolved_tbl, trait_label, is_categorical) {
  value_col <- paste0(trait_label, "_value")
  n_col     <- paste0(trait_label, "_n_records")
  src_col   <- paste0(trait_label, "_src")
  
  gift_tr <- gift_matched %>% filter(trait == trait_label)
  if (nrow(gift_tr) == 0) return(resolved_tbl)
  
  if (is_categorical) {
    gift_tr <- gift_tr %>% mutate(valor = harmonize_gift_growth_form(valor))
  } else {
    gift_tr <- gift_tr %>% mutate(valor = as.numeric(valor) * gift_unit_conversion[[trait_label]])
  }
  
  gift_tr <- gift_tr %>%
    filter(!is.na(valor)) %>%
    distinct(match_species, .keep_all = TRUE) %>%
    transmute(species = match_species, gift_value = valor, gift_src = src)
  
  joined <- resolved_tbl %>%
    select(species, !!value_col, !!n_col, !!src_col) %>%
    left_join(gift_tr, by = "species") %>%
    mutate(has_try = .data[[src_col]] != "unresolved", has_gift = !is.na(gift_value))
  
  if (is_categorical) {
    n_agree    <- joined %>% filter(has_try, has_gift, .data[[value_col]] == gift_value) %>% nrow()
    n_conflict <- joined %>% filter(has_try, has_gift, .data[[value_col]] != gift_value) %>% nrow()
    n_gift_only <- joined %>% filter(!has_try, has_gift) %>% nrow()
    
    conflicts <- joined %>%
      filter(has_try, has_gift, .data[[value_col]] != gift_value) %>%
      transmute(species, try_value = .data[[value_col]], gift_value)
    growth_form_conflicts <<- bind_rows(growth_form_conflicts, conflicts)
    
    cat(" ", trait_label, ": ", n_agree, "agree, ", n_conflict,
        "conflict (TRY kept, logged for review), ", n_gift_only, "GIFT-only\n", sep = "")
    
    out <- joined %>%
      mutate(
        !!value_col := case_when(
          has_try ~ .data[[value_col]],                      # TRY present: keep it (agree or conflict)
          has_gift ~ gift_value,                               # TRY missing, GIFT present
          TRUE ~ .data[[value_col]]
        ),
        !!src_col := case_when(
          has_try & has_gift & .data[[value_col]] == gift_value ~ "try+gift:agree",
          has_try & has_gift                                    ~ paste0(.data[[src_col]], "+gift:conflict"),
          has_try                                               ~ .data[[src_col]],
          has_gift                                              ~ "gift_only",
          TRUE                                                  ~ .data[[src_col]]
        ),
        !!n_col := if_else(!has_try & has_gift, 1L, .data[[n_col]])
      )
  } else {
    # Fold-change disagreement between TRY and GIFT (unit-converted), each
    # direction treated symmetrically (10x too high and 10x too low are
    # both a 10-fold disagreement) -- see NUMERIC_CONFLICT_RATIO_FOLD above.
    joined <- joined %>%
      mutate(ratio = if_else(has_try & has_gift, gift_value / .data[[value_col]], NA_real_),
             is_conflict = has_try & has_gift &
               (ratio > NUMERIC_CONFLICT_RATIO_FOLD | ratio < 1 / NUMERIC_CONFLICT_RATIO_FOLD))
    
    n_avg       <- joined %>% filter(has_try, has_gift, !is_conflict) %>% nrow()
    n_conflict  <- joined %>% filter(is_conflict) %>% nrow()
    n_gift_only <- joined %>% filter(!has_try, has_gift) %>% nrow()
    cat(" ", trait_label, ": ", n_avg, "averaged (TRY+GIFT), ", n_conflict,
        "conflict (>", NUMERIC_CONFLICT_RATIO_FOLD, "-fold disagreement, TRY kept, logged), ",
        n_gift_only, "GIFT-only\n", sep = "")
    
    conflicts <- joined %>%
      filter(is_conflict) %>%
      transmute(species, trait = trait_label, try_value = .data[[value_col]], gift_value, ratio)
    numeric_conflicts <<- bind_rows(numeric_conflicts, conflicts)
    
    out <- joined %>%
      mutate(
        !!value_col := case_when(
          has_try & has_gift & !is_conflict ~ (.data[[value_col]] + gift_value) / 2,
          has_try                           ~ .data[[value_col]],   # includes conflicts -- TRY kept
          has_gift                          ~ gift_value,
          TRUE                               ~ .data[[value_col]]
        ),
        !!src_col := case_when(
          has_try & has_gift & is_conflict  ~ paste0(.data[[src_col]], "+gift:conflict"),
          has_try & has_gift                ~ "try+gift_avg",
          has_try                           ~ .data[[src_col]],
          has_gift                          ~ "gift_only",
          TRUE                               ~ .data[[src_col]]
        ),
        !!n_col := case_when(
          has_try & has_gift & !is_conflict ~ .data[[n_col]] + 1L,  # +1 for the GIFT value folded in
          has_try                           ~ .data[[n_col]],       # unchanged for conflicts -- GIFT not folded in
          has_gift                          ~ 1L,
          TRUE                               ~ .data[[n_col]]
        )
      )
  }
  
  out <- out %>% select(species, !!value_col, !!n_col, !!src_col)
  resolved_tbl %>% select(-!!value_col, -!!n_col, -!!src_col) %>% left_join(out, by = "species")
}

cat("\nTRY + GIFT combination, per applicable trait:\n")
trait_species_resolved_try_gift <- trait_species_resolved
for (tr in c("SLA", "seed_mass")) {
  trait_species_resolved_try_gift <- combine_try_gift(trait_species_resolved_try_gift, tr, is_categorical = FALSE)
}
trait_species_resolved_try_gift <- combine_try_gift(trait_species_resolved_try_gift, "growth_form", is_categorical = TRUE)

if (nrow(growth_form_conflicts) > 0) {
  write_csv(growth_form_conflicts, this_module_dir("tmp_files", "gift_try_growth_form_conflicts.csv"))
  cat("\nTRY/GIFT growth_form conflicts (TRY kept, written for manual review):",
      nrow(growth_form_conflicts), "species\n")
}

if (nrow(numeric_conflicts) > 0) {
  write_csv(numeric_conflicts, this_module_dir("tmp_files", "gift_try_numeric_conflicts.csv"))
  cat("\nTRY/GIFT numeric conflicts (>", NUMERIC_CONFLICT_RATIO_FOLD, "-fold disagreement, TRY kept, ",
      "written for manual review):\n", sep = "")
  print(numeric_conflicts)
}

write_csv(trait_species_resolved_try_gift, this_module_dir("tmp_files", "try_gift_species_trait_summary_resolved.csv"))

# --- 6b. GIFT raw records, unit-converted/harmonized, for genus-level pooling --
# Genus-level imputation (05_impute_trait_by_genus.R) pools TRY + GIFT raw
# records when computing a genus mean/mode for the 3 traits in GIFT (SLA,
# seed_mass, growth_form; see that script's header). Built from the full
# GIFT database (gift_long), not only the records matched to project
# species (gift_matched), so that a genus mean uses every GIFT record of
# the genus, as the TRY raw tables (seed_mass_raw, sla_raw,
# growth_form_raw) are built from the full TRY database.
# Numeric and categorical traits are handled in two separate filter +
# mutate steps, so that as.numeric() is never applied to growth_form text
# values ("herb", "shrub", ...), which would raise "NAs introduced by
# coercion" warnings (if_else()/case_when() evaluate both branches across
# the whole vector).
gift_numeric_raw <- gift_long %>%
  filter(trait %in% names(gift_unit_conversion)) %>%
  mutate(value = as.numeric(valor) * gift_unit_conversion[trait]) %>%
  filter(!is.na(value))

gift_categorical_raw <- gift_long %>%
  filter(trait == "growth_form") %>%
  mutate(value = harmonize_gift_growth_form(valor)) %>%
  filter(!is.na(value))

gift_raw_by_trait <- list(
  seed_mass   = gift_numeric_raw     %>% filter(trait == "seed_mass") %>% transmute(species, value),
  SLA         = gift_numeric_raw     %>% filter(trait == "SLA")       %>% transmute(species, value),
  growth_form = gift_categorical_raw %>% transmute(species, value)
)

cat("\nGIFT raw records available for genus-level pooling (full GIFT database, pre project-species match):\n")
for (tr in names(gift_raw_by_trait)) {
  cat(" ", tr, ":", nrow(gift_raw_by_trait[[tr]]), "records,",
      n_distinct(gift_raw_by_trait[[tr]]$species), "species,",
      n_distinct(word(gift_raw_by_trait[[tr]]$species, 1)), "genera\n")
}

cat("\nUnresolved after TRY + GIFT (candidates for genus imputation), per trait:\n")
for (tr in trait_labels) {
  n_unresolved <- sum(trait_species_resolved_try_gift[[paste0(tr, "_src")]] == "unresolved", na.rm = TRUE)
  cat(" ", tr, ":", n_unresolved, "species\n")
}

# --- 6. Re-save the intermediate .rds with the TRY+GIFT table added -----------
saveRDS(
  list(
    seed_mass_raw = seed_mass_raw, longevity_raw = longevity_raw,
    sla_raw = sla_raw, wood_dens_raw = wood_dens_raw,
    growth_form_raw = growth_form_raw, woodiness_raw = woodiness_raw,
    gift_raw_by_trait = gift_raw_by_trait,  # GIFT raw records (SLA/seed_mass/growth_form), for genus pooling (05)
    trait_species_resolved = trait_species_resolved,                    # TRY only (kept for reference/validation)
    trait_species_resolved_try_gift = trait_species_resolved_try_gift,  # TRY + GIFT combined -- used from here on
    project_species = project_species,
    overstory_species = overstory_species,
    understory_species = understory_species,
    spp_synonyms = spp_synonyms
  ),
  this_module_dir("tmp_files", "trait_matching_intermediate.rds")
)
cat("\nSaved tmp_files/trait_matching_intermediate.rds with trait_species_resolved_try_gift and gift_raw_by_trait added.\n")
