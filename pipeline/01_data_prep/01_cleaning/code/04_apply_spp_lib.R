# -------------------------------------------------------------------------
# Author: Agostina Torres
# Project: MORAclimshift
# Script: 04_apply_spp_lib.R
# Description:
# Gathers every unique raw species code/name found in the Tree and Veg
# dataframes (Legacy + Resurvey), looks each one up against spp_lib_legacy
# (past_code) or spp_lib_current (prst_code) as appropriate, and diagnoses
# whether it resolved to a final_species or not.
#
# Uses tree_R_clean/tree_L_clean/veg_R_clean/veg_L_clean (from
# 02_apply_plot_ids.R, in tree_veg_plotid_clean.rda) rather than the raw
# _Org.csv files, so that species are gathered from the same 332 confirmed
# plots used in the rest of the pipeline, with Plot_ID already resolved.
#
# Species are not restricted to those recorded in both periods: every
# species found in either period is included. The purpose of this step is
# to check that the raw codes of the survey data are represented in the
# species list, not to filter by temporal coverage.
#
# Format differences per dataset:
#   tree_R, veg_R : "Species" column, format "CODE-Scientific name" ->
#                   the code before the first hyphen is looked up via
#                   prst_code in spp_lib_current
#   tree_L        : "Species" column is the 4-letter code -> looked up via
#                   past_code in spp_lib_legacy
#   veg_L         : wide format, species are column names -> each column
#                   name is a past_code, looked up via spp_lib_legacy
#
# Outputs (pipeline/01_data_prep/01_cleaning/outputs/):
#   spp_diagnostic.csv - every unique raw code found, per dataset x period,
#                         with the resolved final_species (or NA)
#   spp_needs_review.csv - codes that did not resolve
#   spp_lib.csv - one row per raw-code variant, columns split by dataset x
#                 period (species_org_legacy_tree, _veg, _resurvey_tree, _veg)
#   spp_list_final.csv - unique final_species found across everything
# -------------------------------------------------------------------------

library(tidyverse)
library(here)
library(TNRS)

outputs_dir <- here("pipeline", "01_data_prep", "01_cleaning", "outputs")
tmp_dir     <- here("pipeline", "01_data_prep", "01_cleaning", "tmp_files")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# STEP 1 — Load the plot-id-cleaned data (built by 02_apply_plot_ids.R) +
# species lookup tables
# ==============================================================================
# The plot-id-cleaned data are used (not the raw _Org.csv files) so that
# species are gathered from the same 332 confirmed plots, with Plot_ID
# already resolved.

load(file.path(outputs_dir, "tree_veg_plotid_clean.rda"))
# brings in: tree_R_clean, tree_L_clean, veg_R_clean, veg_L_clean

spp_lib_current <- read_csv(file.path(outputs_dir, "spp_lib_current.csv"), show_col_types = FALSE)
spp_lib_legacy  <- read_csv(file.path(outputs_dir, "spp_lib_legacy.csv"),  show_col_types = FALSE)


# ==============================================================================
# STEP 2 — Gather unique raw species codes per dataset x period
# ==============================================================================

# tree_R / veg_R: "Species" = "CODE-Scientific name" -> extract the code
extract_code <- function(x) str_squish(str_extract(x, "^[^-]+"))

present_tree_codes <- tree_R_clean |>
  distinct(species_raw = Species) |>
  filter(!is.na(species_raw)) |>
  mutate(species_code = extract_code(species_raw), period = "present", source_dataset = "tree")

present_veg_codes <- veg_R_clean |>
  distinct(species_raw = Species) |>
  filter(!is.na(species_raw)) |>
  mutate(species_code = extract_code(species_raw), period = "present", source_dataset = "veg")

# tree_L: "Species" is the code
past_tree_codes <- tree_L_clean |>
  distinct(species_raw = Species) |>
  filter(!is.na(species_raw)) |>
  mutate(species_code = str_squish(species_raw), period = "past", source_dataset = "tree")

# veg_L: wide format, column names are the codes (excluding non-species
# columns). Many veg_L column names have leading/trailing spaces in the CSV
# header (e.g. "  ACCI  "); str_squish() removes them, which resolves 221 of
# 267 otherwise unresolved codes.
non_species_cols <- c("SiteName", "Plotno", "plot_id_raw", "plot_id")
veg_L_species_cols <- setdiff(names(veg_L_clean), non_species_cols)

past_veg_codes <- tibble(
  species_raw   = veg_L_species_cols,
  species_code  = str_squish(veg_L_species_cols),
  period        = "past",
  source_dataset = "veg"
)

raw_codes <- bind_rows(present_tree_codes, present_veg_codes, past_tree_codes, past_veg_codes)

cat("Unique raw species codes found:\n")
cat("  Tree present:", sum(raw_codes$period == "present" & raw_codes$source_dataset == "tree"), "\n")
cat("  Veg present :", sum(raw_codes$period == "present" & raw_codes$source_dataset == "veg"), "\n")
cat("  Tree past   :", sum(raw_codes$period == "past" & raw_codes$source_dataset == "tree"), "\n")
cat("  Veg past    :", sum(raw_codes$period == "past" & raw_codes$source_dataset == "veg"), "\n")


# ==============================================================================
# STEP 2b — Pre-correction of typos in VegData Resurvey species codes,
# identified from the review of spp_needs_review.csv
#
# Batch 1: single-letter typos and transpositions, unambiguous.

species_name_corrections <- tribble(
  ~species_code_raw, ~corrected_code, ~reason,
  "OLPHOR",  "OPLHOR", "letters transposed (L<->P)",
  "PACMYR",  "PAXMYR", "letter C->X",
  "PAIMYR",  "PAXMYR", "letter I->X",
  "TSHUET",  "TSUHET", "letters transposed (H<->U)",
  "TSUMEN",  "TSUMER", "last letter N->R",
  "VACMEN",  "VACMEM", "last letter N->M",
  "VIOGIA",  "VIOGLA", "letters transposed (I<->L)",
  "POLMAN",  "POLMUN", "letter A->U",
  "PSUMEN",  "PSEMEN", "letters transposed (U<->E)",
  "THUPLE",  "THUPLI", "last letter E->I",
  "BELSPI",  "BLESPI", "letters transposed (E<->L)",
  "BLESPE",  "BLESPI", "last letter E->I",
  "MEMFER",  "MENFER", "letters transposed (M<->N)",
  "MENFEK",  "MENFER", "letter K->R",
  "MENFIR",  "MENFER", "letter I->E",
  "TEATRI",  "TIATRI", "letters transposed (E<->I)",
  "TIATARI", "TIATRI", "extra letter A",
  "TRIATRI", "TIATRI", "letters transposed at start",
  # Batch 2: checked against the valid 2023 prst_codes
  "ACTTRI",  "ACHTRI", "letters transposed (C<->H)",
  "CHIMAP",  "CHIMEN", "truncated/altered ending",
  "CLIUMI",  "CLIUNI", "letter M->N",
  "HEIALB",  "HIEALB", "letters transposed (E<->I)",
  "LUMMAR",  "LOMMAR", "letter U->O",
  "MUEMAC",  "MOEMAC", "letter U->O",
  "PETPAL",  "PETFRI", "different ending altogether (name updated)",
  "STRALP",  "STRAMP", "letters transposed (L<->M... P kept)",
  "STRALT",  "STRAMP", "different ending altogether",
  # Batch 3: genus reclassification Listera -> Neottia (POWO). "Neottia
  # caurina" is not an accepted name but a synonym of Neottia banksiana, so
  # LISCAU maps to NEOBAN.
  "LISCAU",  "NEOBAN",  "genus reclassified Listera->Neottia; caurina further merged into banksiana",
  "LISCON",  "NEOCON",  "genus reclassified Listera->Neottia (convallarioides)",
  "LISCOR",  "NEOCOR",  "genus reclassified Listera->Neottia (cordata)"
)

raw_codes <- raw_codes |>
  left_join(species_name_corrections, by = c("species_code" = "species_code_raw")) |>
  mutate(species_code_lookup = coalesce(corrected_code, species_code))


# ==============================================================================
# STEP 3 — Look up final_species: present -> spp_lib_current, past -> spp_lib_legacy
# ==============================================================================

present_diagnostic <- raw_codes |>
  filter(period == "present") |>
  left_join(spp_lib_current, by = c("species_code_lookup" = "prst_code"))

past_diagnostic <- raw_codes |>
  filter(period == "past") |>
  left_join(spp_lib_legacy, by = c("species_code_lookup" = "past_code"))

# Fallback for rows still unresolved after the code-based lookup: in some
# rows the raw text is the full scientific name instead of a code (no
# hyphen, e.g. "Polystichum munitum" instead of "POLMUN-Polystichum
# munitum"), so extract_code() takes the whole phrase as the code. If that
# phrase matches a known final_species (case-insensitive), it is resolved
# directly.
known_species_names <- spp_lib_current |> distinct(final_species, final_family)

resolve_by_full_name <- function(df) {
  df |>
    mutate(name_match_key = str_to_lower(str_squish(species_raw))) |>
    left_join(
      known_species_names |> mutate(name_match_key = str_to_lower(final_species)),
      by = "name_match_key", suffix = c("", "_by_name")
    ) |>
    mutate(
      final_species = coalesce(final_species, final_species_by_name),
      final_family  = coalesce(final_family, final_family_by_name)
    ) |>
    select(-name_match_key, -final_species_by_name, -final_family_by_name)
}

present_diagnostic <- resolve_by_full_name(present_diagnostic)
past_diagnostic    <- resolve_by_full_name(past_diagnostic)

spp_diagnostic_pre_tnrs <- bind_rows(present_diagnostic, past_diagnostic)

# ==============================================================================
# STEP 3b — TNRS resolution: rows still unresolved after the code lookup
# and the full-name fallback, but with a name attached to the raw code
# (e.g. "DOEJEF-Dodecatheon jeffryi"), are checked against TNRS to see
# whether the name resolves to a species already in the catalog under
# another name or code (e.g. Dodecatheon jeffryi -> Primula jeffreyi, a
# genus reclassification). The TNRS-accepted name is applied when it
# resolves and matches a known species. This generalizes the Listera ->
# Neottia corrections of Step 2b.
# ==============================================================================

known_species_for_tnrs_fix <- spp_lib_current |> distinct(final_species, final_family)

still_unresolved_with_name <- spp_diagnostic_pre_tnrs |>
  filter(is.na(final_species), str_detect(species_raw, "-")) |>
  mutate(name_part = str_squish(str_remove(species_raw, "^[^-]+-"))) |>
  distinct(name_part)

if (nrow(still_unresolved_with_name) > 0) {
  tnrs_prepass <- TNRS(still_unresolved_with_name$name_part)
  
  tnrs_fix_lookup <- still_unresolved_with_name |>
    left_join(select(tnrs_prepass, Name_submitted, Accepted_name),
              by = c("name_part" = "Name_submitted")) |>
    inner_join(known_species_for_tnrs_fix, by = c("Accepted_name" = "final_species")) |>
    transmute(name_part, tnrs_final_species = Accepted_name, tnrs_final_family = final_family)
  
  cat("\nApplying", nrow(tnrs_fix_lookup), "TNRS-confirmed corrections",
      "(name resolves to a species already in the catalog):\n")
  print(tnrs_fix_lookup)
  
  spp_diagnostic_pre_tnrs <- spp_diagnostic_pre_tnrs |>
    mutate(name_part = if_else(str_detect(species_raw, "-"),
                               str_squish(str_remove(species_raw, "^[^-]+-")), NA_character_)) |>
    left_join(tnrs_fix_lookup, by = "name_part") |>
    mutate(
      final_species = coalesce(final_species, tnrs_final_species),
      final_family  = coalesce(final_family, tnrs_final_family)
    ) |>
    select(-name_part, -tnrs_final_species, -tnrs_final_family)
}

spp_diagnostic <- spp_diagnostic_pre_tnrs |>
  mutate(resolved = !is.na(final_species)) |>
  select(species_raw, species_code, corrected_code, correction_reason = reason,
         period, source_dataset, final_species, final_family, resolved)

write_csv(spp_diagnostic, file.path(tmp_dir, "spp_diagnostic.csv"))

cat("\nSaved spp_diagnostic.csv:", nrow(spp_diagnostic), "rows\n")
cat("  Resolved    :", sum(spp_diagnostic$resolved), "\n")
cat("  NOT resolved:", sum(!spp_diagnostic$resolved), "\n")


# ==============================================================================
# STEP 4 — Categorize the unresolved codes: substrate / unidentified / genus
# level / real typos needing review
# ==============================================================================
# Three categories of unresolved codes are expected never to resolve to a
# species:
#   - substrate/non-plant: not a species (moss, litter, log, rock...)
#   - unidentified in the field: a plant the surveyor could not identify
#     (UNK* codes)
#   - genus level: excluded from spp_lib by design (id_level != "Species"
#     in 03_clean_spp_names.R)
# The remaining codes are the candidates for a species-name correction
# lookup (as plot_name_corrections.csv for plots).

substrate_patterns <- c(
  "moss", "moss spp\\.?", "litter", "log", "log/litter", "rock", "stump",
  "bare soil", "bareground", "soil", "water", "duff", "germinant", "germinants",
  "none", "n/?a", "\\(blank\\)"
)
substrate_regex <- str_c("^(", str_c(substrate_patterns, collapse = "|"), ")$")

unresolved <- spp_diagnostic |> filter(!resolved)

is_substrate     <- str_detect(str_to_lower(unresolved$species_code), regex(substrate_regex, ignore_case = TRUE))
is_unidentified  <- str_detect(unresolved$species_code, regex("^UNK", ignore_case = TRUE))

spp_substrate    <- unresolved[is_substrate, ]
spp_unidentified <- unresolved[is_unidentified & !is_substrate, ]

# Genus level: a code can indicate "genus known, species unknown" in three
# ways:
#   (a) tagged as such in PlantList (Id Level != "Species")
#   (b) code ending in "SPP" (e.g. CORSPP, VERSPP - genus abbreviation +
#       "species unknown" suffix), with no PlantList entry of its own
#   (c) code that is a bare genus name or abbreviation matching the first
#       word of a known species binomial in spp_list (e.g. VIOLA, CAREX,
#       TSUGA, RIBES), also with no PlantList entry of its own
plant_list <- read_csv(file.path(here("pipeline", "01_data_prep", "01_cleaning",
                                      "data", "data_raw", "01_raw_surveys"),
                                 "MORA_PlantList_AllSurvey_Org.csv"), show_col_types = FALSE)
genus_level_codes_current <- plant_list |>
  filter(`Id Level` != "Species") |>
  pull(SpeciesCode2023) |>
  na.omit()

# Same check for the past code system (1975 4 Letter Code): a code can be
# genus level under the 1975 code even if its current 2023 code is species
# level, since the two systems do not always align 1:1.
genus_level_codes_past <- plant_list |>
  filter(`Id Level` != "Species") |>
  pull(`1975 4 Letter Code`) |>
  na.omit()

genus_level_codes <- union(genus_level_codes_current, genus_level_codes_past)

known_genera <- spp_lib_current |>
  mutate(genus = word(final_species, 1)) |>
  pull(genus) |>
  unique() |>
  str_to_upper()

is_tagged_genus <- unresolved$species_code %in% genus_level_codes
is_spp_suffix   <- str_detect(unresolved$species_code, regex("SPP$", ignore_case = TRUE))
is_bare_genus   <- map_lgl(unresolved$species_code, function(code) {
  code_upper <- str_to_upper(code)
  if (nchar(code_upper) < 4) return(FALSE)
  any(str_starts(known_genera, fixed(code_upper))) || code_upper %in% known_genera
})

is_genus_level <- (is_tagged_genus | is_spp_suffix | is_bare_genus) & !is_substrate & !is_unidentified
spp_genus_level <- unresolved[is_genus_level, ]

# Ambiguity exclusion: MITEL and VIOL are 1975 past_codes that correspond
# to more than one current species in the project's species reference
# table (MITEL to 3 species, Mitella breweri/pentandra/trifida; VIOL to 2,
# Viola adunca/orbiculata). Unlike COCA/PYEL/MADI (single-destination
# corrections; see 03_clean_spp_names.R Step 6), these two cannot be
# resolved to one species in the 1975 data, so their raw codes are removed
# from past_code_lookup. This applies only to the past: present-day data
# for these species (MITBRE/MITPEN/MITTRI/VIOADU/VIOORB) resolve normally
# and are not excluded.
ambiguity_complex_codes <- c("MITEL", "VIOL")
is_ambiguity_exclusion <- unresolved$species_code %in% ambiguity_complex_codes &
  !is_substrate & !is_unidentified & !is_genus_level
spp_excluded_ambiguous <- unresolved[is_ambiguity_exclusion, ]

spp_needs_review <- unresolved[!is_substrate & !is_unidentified & !is_genus_level & !is_ambiguity_exclusion, ]

write_csv(spp_substrate,          file.path(tmp_dir, "spp_substrate.csv"))
write_csv(spp_unidentified,       file.path(tmp_dir, "spp_unidentified.csv"))
write_csv(spp_genus_level,        file.path(tmp_dir, "spp_genus_level.csv"))
write_csv(spp_excluded_ambiguous, file.path(tmp_dir, "spp_excluded_ambiguous.csv"))
write_csv(spp_needs_review,       file.path(tmp_dir, "spp_needs_review.csv"))

cat("\nUnresolved codes broken down:\n")
cat("  Substrate/non-plant           :", nrow(spp_substrate), "\n")
cat("  Unidentified in field         :", nrow(spp_unidentified), "\n")
cat("  Genus-level (excluded)        :", nrow(spp_genus_level), "\n")
cat("  MITEL/VIOL ambiguity exclusion:", nrow(spp_excluded_ambiguous), "\n")
cat("  Real typos needing review     :", nrow(spp_needs_review), "\n")


# ==============================================================================
# Codes with a name attached that remain in spp_needs_review after the
# TNRS resolution of Step 3b did not match an existing catalog species.


# ==============================================================================
# STEP 4b — spp_master: every code in the catalog (whether or not it
# appears in Tree/Veg), with its observations where it does
#
# Two layers, joined together:
#   1. Catalog layer (spp_catalog_current.csv / spp_catalog_legacy.csv,
#      from 03_clean_spp_names.R): every code in the species list,
#      including genus-level entries. The two ambiguous 1975 codes
#      (MITEL/VIOL) are absent from spp_catalog_legacy, since they were
#      removed at the source.
#   2. Observation layer (spp_diagnostic, built above): for codes found in
#      the Tree/Veg data, the raw spelling, dataset/period and resolution.
# A code with no observation row is still in spp_master. This table is the
# reference for checking any other dataset (seedlings, presence/absence,
# etc.).

spp_catalog_current <- read_csv(file.path(outputs_dir, "spp_catalog_current.csv"), show_col_types = FALSE)
spp_catalog_legacy  <- read_csv(file.path(outputs_dir, "spp_catalog_legacy.csv"),  show_col_types = FALSE)

catalog_long <- bind_rows(
  spp_catalog_current |> transmute(code = prst_code, code_system = "current",
                                   catalog_species = final_species, catalog_family = final_family,
                                   catalog_status),
  spp_catalog_legacy |> transmute(code = past_code, code_system = "past",
                                  catalog_species = final_species, catalog_family = final_family,
                                  catalog_status)
)

decision_lookup <- bind_rows(
  spp_substrate          |> mutate(decision = "excluded_substrate"),
  spp_unidentified       |> mutate(decision = "excluded_unidentified_in_field"),
  spp_genus_level        |> mutate(decision = "excluded_genus_level"),
  spp_excluded_ambiguous |> mutate(decision = "excluded_ambiguity_complex"),
  spp_needs_review       |> mutate(decision = "needs_review")
) |>
  select(species_raw, species_code, corrected_code, period, source_dataset, decision)

observation_layer <- spp_diagnostic |>
  left_join(decision_lookup, by = c("species_raw", "species_code", "corrected_code", "period", "source_dataset")) |>
  mutate(
    decision = if_else(resolved & is.na(decision), "included", decision),
    veg_type = recode(source_dataset, "tree" = "overstory", "veg" = "understory"),
    # code used for the catalog join: the corrected code if there is one,
    # otherwise the raw code as observed
    code = coalesce(corrected_code, species_code),
    code_system = if_else(period == "present", "current", "past")
  ) |>
  select(code, code_system, species_raw, veg_type, period, decision,
         observed_species = final_species, observed_family = final_family)

spp_master <- catalog_long |>
  full_join(observation_layer, by = c("code", "code_system"), relationship = "many-to-many") |>
  mutate(
    # the resolution of spp_diagnostic (the one used in the pipeline)
    # takes priority over the resolution of the original species list;
    # the catalog match is a fallback for codes never observed.
    final_species = coalesce(observed_species, catalog_species),
    final_family  = coalesce(observed_family, catalog_family),
    decision = if_else(is.na(decision) & !is.na(catalog_status),
                       str_c("in_catalog_not_observed_in_tree_veg (", catalog_status, ")"),
                       decision)
  ) |>
  select(code, code_system, catalog_status, final_species, final_family,
         species_raw, veg_type, period, decision)

# MITEL/VIOL 1975-code ambiguity (see 03_clean_spp_names.R Step 6):
# these 5 species resolve in the present-day data but have no resolvable
# past_code, since MITEL/VIOL were removed from past_code_lookup rather
# than assigned to one species. They are flagged here so that spp_master
# documents why they are excluded from spp_lib/spp_list_final below
# (otherwise they would appear as new species in a temporal comparison).
species_excluded_temporal_asymmetry <- c(
  "Mitella breweri", "Mitella pentandra", "Mitella trifida",  # MITEL complex
  "Viola adunca", "Viola orbiculata"                          # VIOL complex
)

spp_master <- spp_master |>
  mutate(decision = if_else(
    final_species %in% species_excluded_temporal_asymmetry & decision == "included",
    "included_present_only_excluded_from_final (1975 MITEL/VIOL ambiguity)",
    decision
  ))


# ==============================================================================
# STEP 4c — tree_size_category: 1975 tree codes with a trailing digit
# (ABAM, ABAM1, ABAM2, ABAM3, ABAM4) encode size class (seedling/sapling/
# pole/mature), not a different species. The digit is split into its own
# column, and it is checked that the base code (without the digit) also
# exists in the catalog and maps to the same species.
# ==============================================================================

base_code_species <- spp_master |>
  filter(code_system == "past", catalog_status == "species_level") |>
  distinct(code, final_species) |>
  rename(base_code = code, base_species = final_species)

spp_master <- spp_master |>
  mutate(
    has_size_suffix = code_system == "past" & str_detect(code, "^[A-Z]{3,6}[1-4]$"),
    code_base = if_else(has_size_suffix, str_remove(code, "[1-4]$"), NA_character_),
    size_digit = if_else(has_size_suffix, str_extract(code, "[1-4]$"), NA_character_)
  ) |>
  left_join(base_code_species, by = c("code_base" = "base_code")) |>
  mutate(
    tree_size_category = if_else(
      has_size_suffix & !is.na(base_species) & base_species == final_species,
      size_digit, NA_character_
    )
  ) |>
  select(-has_size_suffix, -code_base, -size_digit, -base_species)

cat("\nCodes classified with a tree_size_category (verified against base code):\n")
print(spp_master |> filter(!is.na(tree_size_category)) |> distinct(code, final_species, tree_size_category))


# ==============================================================================
# STEP 4d — TNRS audit: species_to_check / species_tnrs / genus_tnrs /
# family_tnrs, for every row with a name
#
# Every species_raw with a name part (after the code and hyphen) is checked
# against TNRS, so that an outdated name in an already included row is also
# detected. These 4 columns are kept separate from final_species/
# final_family (the code-based resolution used in the pipeline): they are
# an independent cross-check, so that disagreements are visible rather than
# overwritten.

spp_master <- spp_master |>
  mutate(species_to_check = case_when(
    # a name is attached to the raw code (present-period "CODE-Name" style)
    str_detect(species_raw, "-") ~ str_squish(str_remove(species_raw, "^[^-]+-")),
    # no name attached (e.g. past-period bare codes like "ABAM") but the
    # code was resolved: the resolved name is checked, so that past-only rows
    # are also audited
    !is.na(final_species) ~ final_species,
    TRUE ~ NA_character_
  ))

names_to_check <- unique(na.omit(spp_master$species_to_check))

if (length(names_to_check) > 0) {
  tnrs_audit <- TNRS(names_to_check) |>
    select(Name_submitted, Accepted_name, Accepted_family) |>
    distinct(Name_submitted, .keep_all = TRUE) |>
    mutate(
      # A species-level match has a genus and a specific epithet (2+ words).
      # A genus-only or failed match (empty, or a single word) is not species
      # level: species_tnrs is NA in that case, while the genus/family found
      # are still reported.
      is_species_level = !is.na(Accepted_name) & str_count(str_squish(Accepted_name), "\\S+") >= 2,
      species_tnrs = if_else(is_species_level, Accepted_name, NA_character_),
      # Binomial-only version: genus + specific epithet, dropping any
      # variety/subspecies/form qualifier (e.g. "Alnus alnobetula subsp.
      # sinuata" -> "Alnus alnobetula"). The first two words; used whenever
      # varieties/subspecies are collapsed to their parent species.
      species_tnrs_binomial = if_else(is_species_level, word(species_tnrs, 1, 2), NA_character_),
      genus_tnrs   = if_else(!is.na(Accepted_name), word(Accepted_name, 1), NA_character_),
      family_tnrs  = Accepted_family
    ) |>
    select(Name_submitted, species_tnrs, species_tnrs_binomial, genus_tnrs, family_tnrs)
  
  spp_master <- spp_master |>
    left_join(tnrs_audit, by = c("species_to_check" = "Name_submitted"))
} else {
  spp_master <- spp_master |> mutate(species_tnrs = NA_character_, species_tnrs_binomial = NA_character_,
                                     family_tnrs = NA_character_, genus_tnrs = NA_character_)
}

cat("\nTNRS audit: rows where species_tnrs DISAGREES with final_species (worth a look):\n")
print(spp_master |>
        filter(!is.na(species_tnrs), !is.na(final_species), species_tnrs != final_species) |>
        distinct(code, final_species, species_tnrs))

spp_master <- spp_master |> arrange(final_species, code_system, code)
write_csv(spp_master, file.path(tmp_dir, "spp_master.csv"))

cat("\nSaved spp_master.csv:", nrow(spp_master), "rows (full catalog + all observations)\n")
cat("Decision breakdown:\n")
print(spp_master |> count(decision))
cat("\nRows still missing final_species despite being 'included':",
    sum(spp_master$decision == "included" & is.na(spp_master$final_species)), "\n")


# ==============================================================================
# STEP 5 — spp_lib: one row per raw-code variant, columns split by dataset x period
# ==============================================================================
# Same structure as plot_lib: when a species has more than one raw
# spelling within the same dataset, each variant gets its own row via
# full_join instead of being combined in one cell.

get_variants <- function(ds, per, colname) {
  spp_diagnostic |>
    filter(source_dataset == ds, period == per, resolved) |>
    distinct(final_species, !!colname := species_raw)
}

legacy_tree   <- get_variants("tree", "past",    "species_org_legacy_tree")
legacy_veg    <- get_variants("veg",  "past",    "species_org_legacy_veg")
resurvey_tree <- get_variants("tree", "present", "species_org_resurvey_tree")
resurvey_veg  <- get_variants("veg",  "present", "species_org_resurvey_veg")

all_species <- spp_diagnostic |>
  filter(resolved) |>
  distinct(final_species)

# spp_info (Step 5c below) uses the full resolved species set, including
# MITEL/VIOL: those species get an info row (with presence_coverage/
# tree_resurvey showing the asymmetry), although they are excluded from
# spp_lib/spp_list_final.
all_species_before_exclusion <- all_species

# Excluded from spp_lib/spp_list_final (but kept, flagged, in spp_master):
# MITEL/VIOL are ambiguous in 1975, so these 5 species resolve in the
# present-day data but have no resolvable past_code. Kept as they are,
# they would appear as 3 + 2 new species in any temporal comparison.
all_species <- all_species |>
  filter(!final_species %in% species_excluded_temporal_asymmetry)

cat("\nExcluded from spp_lib/spp_list_final (present-only resolved due to\n")
cat("1975 MITEL/VIOL ambiguity - would create a false temporal signal):\n")
print(species_excluded_temporal_asymmetry)

# TNRS audit columns from spp_master (a species can appear via several
# codes/rows in spp_master, so the first non-NA value per species is used).
tnrs_by_species <- spp_master |>
  filter(!is.na(final_species)) |>
  group_by(final_species) |>
  summarise(
    species_tnrs          = first(na.omit(species_tnrs)),
    species_tnrs_binomial = first(na.omit(species_tnrs_binomial)),
    .groups = "drop"
  )

# canopy_broad lives on spp_lib_current (attached in 03_clean_spp_names.R).
canopy_by_species <- spp_lib_current |>
  distinct(final_species, canopy_broad)


# ==============================================================================
# STEP 5b — species_presence_by_period: whether each species has at least
# one plot with a presence in each period (a resolved code with no
# occurrence, or an explicit zero-cover record, does not count)
#
# final_species is attached to every row of the 4 cleaned datasets, using
# the species_code_lookup/spp_lib tables built above, and presence is
# checked per dataset:
#   Tree (present and past, long format): every row is a measured tree, so
#     any resolved row at a plot counts as presence.
#   Veg present (long format, "Percent Cover"): a literal "0" is an
#     explicit zero-cover record, not presence. Any other non-NA value
#     ("T", a positive number, "present"/"Present") counts.
#   Veg past (wide format, one column per past_code): presence requires a
#     cover value > 0 in that species' column.

attach_species_long <- function(df, lookup_col_join, resolved_lookup) {
  df |>
    mutate(species_code = extract_code(Species)) |>
    left_join(species_name_corrections, by = c("species_code" = "species_code_raw")) |>
    mutate(species_code_lookup = coalesce(corrected_code, species_code)) |>
    left_join(resolved_lookup, by = setNames(lookup_col_join, "species_code_lookup"))
}

# Tree: every resolved row is a presence.
tree_present_species <- attach_species_long(tree_R_clean, "prst_code", spp_lib_current) |>
  filter(!is.na(final_species)) |>
  distinct(plot_id, final_species) |>
  mutate(period = "present", source_dataset = "tree")

tree_past_species <- tree_L_clean |>
  mutate(species_code = str_squish(Species)) |>
  left_join(species_name_corrections, by = c("species_code" = "species_code_raw")) |>
  mutate(species_code_lookup = coalesce(corrected_code, species_code)) |>
  left_join(spp_lib_legacy, by = c("species_code_lookup" = "past_code")) |>
  filter(!is.na(final_species)) |>
  distinct(plot_id, final_species) |>
  mutate(period = "past", source_dataset = "tree")

# Veg present: explicit zero cover ("0") and missing values are excluded;
# any other value ("T", a positive number as text, "present"/"Present") counts.
veg_present_species <- veg_R_clean |>
  filter(!is.na(`Percent Cover`), str_squish(`Percent Cover`) != "0") |>
  mutate(species_code = extract_code(Species)) |>
  left_join(species_name_corrections, by = c("species_code" = "species_code_raw")) |>
  mutate(species_code_lookup = coalesce(corrected_code, species_code)) |>
  left_join(spp_lib_current, by = c("species_code_lookup" = "prst_code")) |>
  filter(!is.na(final_species)) |>
  distinct(plot_id, final_species) |>
  mutate(period = "present", source_dataset = "veg")

# Veg past: wide format, pivoted longer (using veg_L_species_cols from
# Step 2), keeping value > 0; each column name (a past_code) is then
# resolved to its species.
veg_past_species <- veg_L_clean |>
  select(plot_id, all_of(veg_L_species_cols)) |>
  pivot_longer(cols = -plot_id, names_to = "species_code_raw", values_to = "value") |>
  filter(!is.na(value), value > 0) |>
  mutate(species_code_lookup = str_squish(species_code_raw)) |>
  left_join(spp_lib_legacy, by = c("species_code_lookup" = "past_code")) |>
  filter(!is.na(final_species)) |>
  distinct(plot_id, final_species) |>
  mutate(period = "past", source_dataset = "veg")

species_presence_all <- bind_rows(tree_present_species, tree_past_species,
                                  veg_present_species, veg_past_species)

species_presence_by_period <- species_presence_all |>
  distinct(final_species, period) |>
  group_by(final_species) |>
  summarise(
    has_present = any(period == "present"),
    has_past    = any(period == "past"),
    .groups = "drop"
  ) |>
  mutate(presence_coverage = case_when(
    has_present & has_past ~ "both",
    has_present             ~ "only_present",
    has_past                ~ "only_past"
  )) |>
  select(final_species, presence_coverage)

cat("\nSpecies presence by period (>=1 plot with REAL presence, not just a resolved code):\n")
print(species_presence_by_period |> count(presence_coverage))

# tree_resurvey: as above, restricted to the Tree data (the relevant
# dataset for overstory species).
tree_coverage <- bind_rows(tree_present_species, tree_past_species) |>
  distinct(final_species, period) |>
  group_by(final_species) |>
  summarise(
    has_present = any(period == "present"),
    has_past    = any(period == "past"),
    .groups = "drop"
  ) |>
  mutate(tree_resurvey = case_when(
    has_present & has_past ~ "both",
    has_present             ~ "only_present",
    has_past                ~ "only_past"
  )) |>
  select(final_species, tree_resurvey)

cat("\nTree resurvey coverage for overstory species (dynamic, replaces old hardcoded tiers):\n")
print(
  canopy_by_species |>
    filter(canopy_broad == "overstory") |>
    left_join(tree_coverage, by = "final_species") |>
    arrange(final_species)
)

spp_lib <- all_species |>
  full_join(legacy_tree,   by = "final_species", relationship = "many-to-many") |>
  full_join(legacy_veg,    by = "final_species", relationship = "many-to-many") |>
  full_join(resurvey_tree, by = "final_species", relationship = "many-to-many") |>
  full_join(resurvey_veg,  by = "final_species", relationship = "many-to-many") |>
  left_join(tnrs_by_species, by = "final_species") |>
  arrange(final_species)

write_csv(spp_lib, file.path(outputs_dir, "spp_lib.csv"))
cat("\nSaved spp_lib.csv:", nrow(spp_lib), "rows\n")


# ==============================================================================
# STEP 5c — spp_info: rebuilt here, once tree_resurvey/presence_coverage
# are available. One row per resolved species (same species set as
# all_species, before the MITEL/VIOL exclusion; those 5 species get their
# own row, flagged via presence_coverage/tree_resurvey).
#   canopy_type: over / under (from canopy_broad)
#   tree_resurvey: both / only_past / only_present (Tree data only, the
#     relevant one for overstory species)
#   presence_coverage: both / only_past / only_present (Tree + Veg
#     combined, the relevant one for understory species)
# ==============================================================================

spp_info_base <- spp_lib_current |>
  distinct(species = final_species, family = final_family, life_form, tree_type, canopy_broad) |>
  rename(canopy_type = canopy_broad) |>
  mutate(canopy_type = if_else(canopy_type == "overstory", "over", "under"))

spp_info <- all_species_before_exclusion |>
  rename(species = final_species) |>
  left_join(spp_info_base, by = "species") |>
  left_join(tree_coverage |> rename(species = final_species), by = "species") |>
  left_join(species_presence_by_period |> rename(species = final_species), by = "species") |>
  arrange(species)

write_csv(spp_info, file.path(outputs_dir, "spp_info.csv"))
cat("\nSaved spp_info.csv:", nrow(spp_info), "rows (one per resolved species)\n")


# ==============================================================================
# STEP 5d — overstory_filter.csv: filtering steps for overstory species
#   Filter 1: species with tree_resurvey == "both" (Tree data in both
#     periods).
#   n_plots_past/n_plots_present: number of plots with >= 1 individual of
#     the species per period (from tree_past_species/tree_present_species,
#     Step 5b).
#   Filter 2: species in more than 5 plots in each period, since a
#     temporal comparison based on very few plots per period is unreliable.
#   Filter 3 (manual): a blank column for an ecological judgement,
#     independent of the two numeric filters.
# Saved as a separate reference file.
# ==============================================================================

n_plots_past_by_species <- tree_past_species |>
  count(final_species, name = "n_plots_past")

n_plots_present_by_species <- tree_present_species |>
  count(final_species, name = "n_plots_present")

overstory_filter <- canopy_by_species |>
  filter(canopy_broad == "overstory") |>
  left_join(tree_coverage,               by = "final_species") |>
  left_join(n_plots_past_by_species,     by = "final_species") |>
  left_join(n_plots_present_by_species,  by = "final_species") |>
  mutate(
    n_plots_past    = coalesce(n_plots_past, 0L),
    n_plots_present = coalesce(n_plots_present, 0L),
    filter_1_both_periods = if_else(coalesce(tree_resurvey, "") == "both", "keep", "drop"),
    # Threshold of >5 plots (rather than >=10), so that species such as
    # Abies grandis (9 past / 8 present) and Picea engelmannii (17 past /
    # 9 present) are kept.
    filter_2_min_plots     = if_else(n_plots_past > 5 & n_plots_present > 5, "keep", "drop"),
    filter_3_manual_comment = NA_character_
  ) |>
  select(species = final_species, tree_resurvey, n_plots_past, n_plots_present,
         filter_1_both_periods, filter_2_min_plots, filter_3_manual_comment) |>
  arrange(species)

write_csv(overstory_filter, file.path(outputs_dir, "overstory_filter.csv"))
cat("\nSaved overstory_filter.csv:", nrow(overstory_filter), "rows\n")
cat("  Filter 1 (both periods)   - keep:", sum(overstory_filter$filter_1_both_periods == "keep"), "\n")
cat("  Filter 1 (both periods)   - drop:", sum(overstory_filter$filter_1_both_periods == "drop"), "\n")
cat("  Filter 2 (>5 plots/period) - keep:", sum(overstory_filter$filter_2_min_plots == "keep"), "\n")
cat("  Filter 2 (>5 plots/period) - drop:", sum(overstory_filter$filter_2_min_plots == "drop"), "\n")
cat("  Filter 3: left blank for manual review\n")


# ==============================================================================
# STEP 6 — spp_list_final: unique final_species across all datasets
# (with the TNRS audit; canopy_type/tree_resurvey/presence_coverage are
# in spp_info.csv)
# ==============================================================================

spp_list_final <- all_species |>
  left_join(tnrs_by_species, by = "final_species") |>
  arrange(final_species)
write_csv(spp_list_final, file.path(outputs_dir, "spp_list_final.csv"))
cat("Saved spp_list_final.csv:", nrow(spp_list_final), "unique species\n")
cat("presence_coverage breakdown (from spp_info.csv):\n")
print(spp_info |> count(presence_coverage))


# ==============================================================================
# STEP 7 — Consistency check: a species present in both periods should
# resolve to the same final_species via both code systems (true by
# construction, since spp_lib_legacy is built from the final_species of
# spp_lib_current)
# ==============================================================================

species_present <- spp_diagnostic |> filter(period == "present", resolved) |> distinct(final_species)
species_past    <- spp_diagnostic |> filter(period == "past", resolved) |> distinct(final_species)

cat("\nSpecies found in present:", nrow(species_present), "\n")
cat("Species found in past   :", nrow(species_past), "\n")
cat("Species found in BOTH   :", nrow(inner_join(species_present, species_past, by = "final_species")), "\n")
cat("Species ONLY in present :", nrow(anti_join(species_present, species_past, by = "final_species")), "\n")
cat("Species ONLY in past    :", nrow(anti_join(species_past, species_present, by = "final_species")), "\n")
