# ------------------------------------------------------------------------------
# 03b_build_trait_table_bien.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Third trait source, after TRY (02) and GIFT (03): the Botanical Information
# and Ecology Network (BIEN; Maitner et al. 2018, R package BIEN), for the
# numeric traits SLA and seed mass. Every trait database added to the
# pipeline goes through the same sequence:
#   1. retrieve records for the project species, their synonyms, the genera of
#      the species and the genera of the synonyms;
#   2. check units against the TRY scale;
#   3. combine at species level with the sources already integrated;
#   4. pool the records into the genus-level evidence used for imputation (05);
#   5. validate genus imputation by leave-one-out (09).
# Literature values for single species are not a database and are added only
# in 08_integration.R, to cells still unresolved.
#
# Steps in this script:
#   1. Genera to download: the genus of every project species plus the genera
#      of all their synonyms (e.g. Mahonia for Berberis nervosa). All BIEN
#      records of these genera are downloaded once and cached
#      (data_raw/bien_traits_genus_raw.csv); only genera missing from the
#      cache are queried, unless REFRESH_BIEN.
#   2. Units: BIEN unit strings are converted to the TRY scale (SLA mm2/mg =
#      m2/kg; seed mass mg). Dataset-level corrections documented in
#      dataset_unit_corrections are applied.
#   3. Unit check against TRY: for every BIEN dataset (project PI), species
#      means are compared with the TRY species means of the same species.
#      Datasets whose median ratio departs from 1 by MAX_FOLD or more (with at
#      least MIN_OVERLAP shared species) are excluded and listed for review.
#      Also excluded: non-numeric values, records outside a plausible range
#      for the trait, and records that duplicate a TRY record (same species
#      and value to four significant digits; BIEN includes datasets that are
#      also in TRY, such as LEDA and the Kew Seed Information Database), so
#      the same measurement is not counted twice.
#   4. Species level: the BIEN species mean (records under the vegetation
#      name; if none, records under its synonyms) is combined with the
#      TRY/GIFT value by the rule used for TRY and GIFT in 03:
#        - no TRY/GIFT value: BIEN fills the cell ("bien:<dataset PIs>");
#        - BIEN within MAX_FOLD of the TRY/GIFT value: equal-weight mean of
#          the source means, i.e. (k * current + BIEN) / (k + 1), where k is
#          the number of sources already averaged (2 for "try+gift_avg",
#          otherwise 1); tag "<previous tag>+bien_avg";
#        - larger disagreement: TRY/GIFT value kept, tag
#          "<previous tag>+bien:conflict", logged for review.
#   5. All retained records are saved as bien_raw_by_trait for genus pooling.
#
# Dataset unit corrections:
#   Stevens JT (doi:10.5061/dryad.q2n8p): SLA stored as m2/kg but 70-115 times
#   lower than other datasets for the same species (Pinus ponderosa,
#   Pseudotsuga menziesii, Arctostaphylos patula, Rubus parviflorus): values
#   are cm2/mg and are multiplied by 100.
#
# The script always starts from the TRY + GIFT table written by 03
# (trait_species_resolved_try_gift) and writes its result to a separate
# object (trait_species_resolved_bien), so it can be rerun without combining
# BIEN twice. 05 uses trait_species_resolved_bien.
#
# Inputs : tmp_files/trait_matching_intermediate.rds (02, 03)
# Outputs: data_raw/bien_traits_genus_raw.csv        raw BIEN records (cache)
#          tmp_files/bien_unit_check_datasets.csv    dataset vs TRY ratios
#          tmp_files/bien_records_used.csv           records retained
#          tmp_files/bien_records_excluded.csv       records excluded, with reason
#          tmp_files/bien_try_gift_conflicts.csv     species where BIEN disagrees
#                                                    with TRY/GIFT by > MAX_FOLD
#          tmp_files/trait_matching_intermediate.rds trait_species_resolved_bien
#                                                    and bien_raw_by_trait added
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(BIEN)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

REFRESH_BIEN <- FALSE   # TRUE: query BIEN again for all genera
MIN_OVERLAP  <- 3       # shared species with TRY needed to test a dataset's scale
MAX_FOLD     <- 5       # fold difference used for the dataset check and the species-level combination

bien_traits <- c(SLA = "leaf area per leaf dry mass", seed_mass = "seed mass")

unit_factor <- tribble(
  ~trait,      ~unit,      ~factor,
  "SLA",       "m2.kg-1",  1,        # m2/kg = mm2/mg
  "seed_mass", "mg",       1,
  "seed_mass", "g",        1000
)

dataset_unit_corrections <- tribble(
  ~trait, ~project_pi,   ~dataset_factor,
  "SLA",  "Stevens JT",  100
)

plausible_range <- tribble(
  ~trait,      ~min,   ~max,
  "SLA",       1,      150,     # mm2/mg
  "seed_mass", 0.001,  10000    # mg
)

genus_of <- function(x) word(x, 1)

# ------------------------------------------------------------------------------
# 1. Genera to download
# ------------------------------------------------------------------------------
obj <- readRDS(this_module_dir("tmp_files", "trait_matching_intermediate.rds"))
resolved     <- obj$trait_species_resolved_try_gift
spp_synonyms <- obj$spp_synonyms
if (is.null(resolved)) stop("Run 03_build_trait_table_gift.R first.")
# The TRY + GIFT table must not already contain BIEN values
if (any(str_detect(unlist(select(resolved, ends_with("_src"))), "bien"), na.rm = TRUE)) {
  stop("trait_species_resolved_try_gift already contains BIEN values; ",
       "re-run 03_build_trait_table_gift.R before this script.")
}

unresolved_before <- map_dfr(names(bien_traits), function(tr) {
  resolved %>% filter(.data[[paste0(tr, "_src")]] == "unresolved") %>% transmute(species, trait = tr)
})
cat("Species x trait unresolved after TRY + GIFT:\n")
print(count(unresolved_before, trait))

# Names of each project species: vegetation name and synonyms
species_names <- obj$project_species %>%
  select(species) %>%
  left_join(spp_synonyms %>% select(species = final_species, name = synonym), by = "species") %>%
  bind_rows(obj$project_species %>% transmute(species, name = species)) %>%
  filter(!is.na(name)) %>%
  distinct()

genera <- species_names %>%
  mutate(genus = genus_of(name)) %>%
  distinct(genus) %>% pull(genus) %>% sort()
cat("\nGenera queried in BIEN (project species and their synonyms):", length(genera), "\n")

# ------------------------------------------------------------------------------
# 2. BIEN download by genus (cached)
# ------------------------------------------------------------------------------
cache_path <- this_module_dir("data_raw", "bien_traits_genus_raw.csv")
cache <- if (file.exists(cache_path) && !REFRESH_BIEN) {
  read_csv(cache_path, col_types = cols(.default = col_character()))
} else {
  tibble(query_genus = character(), query_trait = character())
}

to_download <- crossing(genus = genera, trait = names(bien_traits)) %>%
  anti_join(cache %>% distinct(genus = query_genus, trait = query_trait), by = c("genus", "trait"))

if (nrow(to_download) > 0) {
  missing_traits <- setdiff(bien_traits, BIEN_trait_list()$trait_name)
  if (length(missing_traits) > 0) stop("BIEN trait names not found: ", paste(missing_traits, collapse = ", "))
  new_records <- map_dfr(names(bien_traits), function(tr) {
    gn <- to_download$genus[to_download$trait == tr]
    if (length(gn) == 0) return(NULL)
    res <- BIEN_trait_traitbygenus(genus = gn, trait = bien_traits[[tr]])
    queried <- tibble(query_genus = gn, query_trait = tr)
    if (is.null(res) || nrow(res) == 0) return(queried)
    res <- as_tibble(lapply(res, as.character)) %>%
      mutate(query_genus = genus_of(scrubbed_species_binomial), query_trait = tr)
    bind_rows(res, queried %>% anti_join(res, by = c("query_genus", "query_trait")))
  })
  cache <- bind_rows(cache, new_records)
  write_csv(cache, cache_path)
}
for (col in c("scrubbed_species_binomial", "trait_value", "unit", "project_pi", "url_source")) {
  if (!col %in% names(cache)) cache[[col]] <- NA_character_
}

# ------------------------------------------------------------------------------
# 3. Units, dataset check against TRY, plausibility
# ------------------------------------------------------------------------------
records <- cache %>%
  filter(!is.na(trait_value), query_genus %in% genera) %>%
  transmute(bien_species = scrubbed_species_binomial, trait = query_trait,
            raw_value = suppressWarnings(as.numeric(trait_value)), unit, project_pi, url_source) %>%
  left_join(unit_factor, by = c("trait", "unit")) %>%
  left_join(dataset_unit_corrections, by = c("trait", "project_pi")) %>%
  mutate(value = raw_value * factor * coalesce(dataset_factor, 1))

# TRY records on the same scale, used for the dataset check and to detect
# BIEN records that duplicate a TRY record
try_records <- bind_rows(
  obj$sla_raw       %>% transmute(trait = "SLA",       bien_species = species, value),
  obj$seed_mass_raw %>% transmute(trait = "seed_mass", bien_species = species, value)
)
try_keys <- try_records %>% transmute(trait, bien_species, key = signif(value, 4)) %>% distinct()
records <- records %>%
  mutate(key = signif(value, 4)) %>%
  left_join(try_keys %>% mutate(dup_try = TRUE), by = c("trait", "bien_species", "key")) %>%
  mutate(dup_try = coalesce(dup_try, FALSE)) %>%
  select(-key)

# TRY species means (same scale) as the reference for each dataset
try_means <- try_records %>%
  group_by(trait, bien_species) %>%
  summarise(try_mean = mean(value), .groups = "drop")

dataset_check <- records %>%
  filter(!is.na(value)) %>%
  group_by(trait, project_pi, bien_species) %>%
  summarise(bien_mean = mean(value), .groups = "drop") %>%
  inner_join(try_means, by = c("trait", "bien_species")) %>%
  group_by(trait, project_pi) %>%
  summarise(n_overlap = n(), median_ratio = median(bien_mean / try_mean), .groups = "drop") %>%
  mutate(scale_ok = n_overlap < MIN_OVERLAP | abs(log10(median_ratio)) < log10(MAX_FOLD))
write_csv(dataset_check, this_module_dir("tmp_files", "bien_unit_check_datasets.csv"))
cat("\nBIEN datasets compared with TRY (median ratio of species means):\n")
print(dataset_check %>% arrange(trait, desc(n_overlap)), n = Inf, width = Inf)

records <- records %>%
  left_join(dataset_check %>% select(trait, project_pi, scale_ok), by = c("trait", "project_pi")) %>%
  left_join(plausible_range, by = "trait") %>%
  mutate(exclusion = case_when(
    is.na(raw_value)           ~ "non-numeric value",
    is.na(factor)              ~ "unit not recognised",
    !coalesce(scale_ok, TRUE)  ~ "dataset scale differs from TRY",
    value < min | value > max  ~ "outside plausible range",
    dup_try                    ~ "duplicate of TRY record",
    TRUE                       ~ NA_character_
  ))

write_csv(records %>% filter(is.na(exclusion)),  this_module_dir("tmp_files", "bien_records_used.csv"))
write_csv(records %>% filter(!is.na(exclusion)), this_module_dir("tmp_files", "bien_records_excluded.csv"))
cat("\nBIEN records retained:", sum(is.na(records$exclusion)), "| excluded:\n")
print(count(records %>% filter(!is.na(exclusion)), trait, exclusion))

records_ok <- records %>% filter(is.na(exclusion))

# ------------------------------------------------------------------------------
# 4. Species level: combine with TRY/GIFT
# ------------------------------------------------------------------------------
# Records under the vegetation name are used when present; otherwise records
# under its synonyms (same order of preference as TRY and GIFT)
bien_values <- records_ok %>%
  inner_join(species_names, by = c("bien_species" = "name"), relationship = "many-to-many") %>%
  group_by(species, trait) %>%
  filter(if (any(bien_species == species)) bien_species == species else TRUE) %>%
  summarise(bien_value = mean(value), bien_n = n(),
            bien_src = paste0("bien:", paste(sort(unique(project_pi)), collapse = "; ")),
            .groups = "drop")

combined_log <- list()
for (tr in names(bien_traits)) {
  v <- bien_values %>% filter(trait == tr) %>% select(-trait)
  if (nrow(v) == 0) next
  value_col <- paste0(tr, "_value"); n_col <- paste0(tr, "_n_records"); src_col <- paste0(tr, "_src")
  joined <- resolved %>%
    left_join(v, by = "species") %>%
    mutate(
      cur   = as.numeric(.data[[value_col]]),
      has_cur  = .data[[src_col]] != "unresolved" & !is.na(cur),
      has_bien = !is.na(bien_value),
      k     = if_else(.data[[src_col]] == "try+gift_avg", 2, 1),
      ratio = if_else(has_cur & has_bien, bien_value / cur, NA_real_),
      action = case_when(
        !has_bien                                   ~ "none",
        !has_cur                                    ~ "bien_only",
        ratio <= MAX_FOLD & ratio >= 1 / MAX_FOLD   ~ "averaged",
        TRUE                                        ~ "conflict"
      )
    )
  combined_log[[tr]] <- joined %>% filter(action != "none") %>%
    transmute(species, trait = tr, previous_src = .data[[src_col]], previous_value = cur,
              bien_value, bien_n, bien_src, ratio, action)
  resolved <- joined %>%
    mutate(
      !!value_col := case_when(action == "bien_only" ~ bien_value,
                               action == "averaged"  ~ (k * cur + bien_value) / (k + 1),
                               TRUE                  ~ cur),
      !!n_col     := case_when(action == "bien_only" ~ as.integer(bien_n),
                               action == "averaged"  ~ as.integer(.data[[n_col]]) + as.integer(bien_n),
                               TRUE                  ~ as.integer(.data[[n_col]])),
      !!src_col   := case_when(action == "bien_only" ~ bien_src,
                               action == "averaged"  ~ paste0(.data[[src_col]], "+bien_avg"),
                               action == "conflict"  ~ paste0(.data[[src_col]], "+bien:conflict"),
                               TRUE                  ~ .data[[src_col]])
    ) %>%
    select(-bien_value, -bien_n, -bien_src, -cur, -has_cur, -has_bien, -k, -ratio, -action)
}
combined_log <- bind_rows(combined_log)

cat("\nBIEN combined with TRY/GIFT at species level:\n")
print(count(combined_log, trait, action))
cat("\nCells filled by BIEN alone:\n")
print(as.data.frame(combined_log %>% filter(action == "bien_only") %>%
                      select(species, trait, bien_value, bien_n, bien_src)))
conflicts <- combined_log %>% filter(action == "conflict")
write_csv(conflicts, this_module_dir("tmp_files", "bien_try_gift_conflicts.csv"))
if (nrow(conflicts) > 0) {
  cat("\nBIEN vs TRY/GIFT conflicts (> ", MAX_FOLD, "-fold; TRY/GIFT kept):\n", sep = "")
  print(as.data.frame(conflicts %>% select(species, trait, previous_value, bien_value, ratio, bien_src)))
}

cat("\nStill unresolved after TRY + GIFT + BIEN (species level):\n")
for (tr in names(bien_traits)) {
  cat(" ", tr, ":", sum(resolved[[paste0(tr, "_src")]] == "unresolved", na.rm = TRUE), "species\n")
}

# ------------------------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------------------------
obj$trait_species_resolved_bien <- resolved
obj$bien_raw_by_trait <- map(set_names(names(bien_traits)), ~
                               records_ok %>% filter(trait == .x) %>% transmute(species = bien_species, value))
saveRDS(obj, this_module_dir("tmp_files", "trait_matching_intermediate.rds"))
write_csv(resolved, this_module_dir("tmp_files", "try_gift_bien_species_trait_summary_resolved.csv"))
