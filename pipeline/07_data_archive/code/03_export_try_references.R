
# ------------------------------------------------------------------------------
# 03_export_try_references.R
# Module: pipeline/07_data_archive
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: list the original TRY datasets that contributed trait records to
# the species-level trait values in species_traits.csv, so that they can be
# cited in the data archive. Public TRY data are released under CC BY 4.0,
# which requires citing TRY (Kattge et al. 2020) and the original datasets.
#
# A TRY record is counted when it was used for a value in species_traits.csv:
#   - species-level values from TRY (source tag "direct", "synonym:<name>" or
#     "try+gift..."): records of the species under its accepted name or one
#     of its synonyms;
#   - genus-level values ("genus_mean...:<Genus>" or "genus_mode...:<Genus>"):
#     records of any species of that genus, since the genus mean pools them.
# The same record filters as in 02_build_trait_table_try.R are applied
# (trait IDs, ErrorRisk <= 4, non-missing value).
#
# Input:  raw TRY exports (pipeline/01_data_prep/04_traits/data_raw/),
#         archive/species_traits.csv, 04_traits/outputs/spp_synonyms.csv
# Output: archive/try_references.csv, one row per TRY dataset, with the
#         dataset reference, the traits it contributed and the number of
#         records used.
# Run after 01_export_data_archive.R.
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(data.table)

p      <- function(...) here("pipeline", ...)
traits <- function(...) p("01_data_prep", "04_traits", ...)
out    <- function(...) p("07_data_archive", "archive", ...)

# 1. TRY files and trait IDs, as in 02_build_trait_table_try.R ------------------
try_sources <- tribble(
  ~trait,         ~file,       ~trait_id,
  "seed_mass",    "38172.txt", list(26),
  "longevity",    "40000.txt", list(59),
  "SLA",          "40000.txt", list(c(3115, 3116, 3117)),
  "wood_density", "52607.txt", list(4),
  "growth_form",  "52732.txt", list(42),
  "woodiness",    "52732.txt", list(38)
)

# 2. Species and genera whose TRY records were used, per trait --------------------
species_traits <- read_csv(out("species_traits.csv"), show_col_types = FALSE)

extra_synonyms <- tribble(
  ~final_species,                ~synonym,
  "Berberis aquifolium",         "Mahonia aquifolium",
  "Berberis nervosa",            "Mahonia nervosa",
  "Monotropa hypopitys",         "Hypopitys monotropa",
  "Lactuca muralis",             "Mycelis muralis",
  "Trautvetteria caroliniensis", "Trautvetteria carolinensis",
  "Mitella trifida",             "Ozomelis trifida"
)
synonyms <- read_csv(traits("outputs", "spp_synonyms.csv"), show_col_types = FALSE) %>%
  select(final_species, synonym) %>%
  bind_rows(extra_synonyms) %>%
  distinct()

sources_long <- species_traits %>%
  select(species, ends_with("_src")) %>%
  pivot_longer(-species, names_to = "trait", values_to = "src") %>%
  mutate(trait = str_remove(trait, "_src$")) %>%
  filter(!is.na(src)) %>%
  distinct()

# Species-level TRY values: accepted name plus synonyms
used_species <- sources_long %>%
  filter(str_detect(src, "^(direct|synonym|try)")) %>%
  select(trait, species) %>%
  left_join(synonyms, by = c("species" = "final_species"), relationship = "many-to-many") %>%
  transmute(trait, name = species, synonym) %>%
  pivot_longer(c(name, synonym), names_to = "source", values_to = "name") %>%
  filter(!is.na(name)) %>%
  distinct(trait, name)

# Genus-level values: every species of the genus
used_genera <- sources_long %>%
  filter(str_detect(src, "^genus")) %>%
  transmute(trait, genus = str_extract(src, "[^:]+$")) %>%
  distinct()

# 3. Read the TRY records and keep those used ------------------------------------
read_used_records <- function(trait_name, file, trait_id) {
  path <- traits("data_raw", file)
  if (!file.exists(path)) { cat("File not found:", path, "\n"); return(NULL) }
  dt <- fread(path, sep = "\t", header = TRUE, quote = "", encoding = "Latin-1",
              showProgress = FALSE)
  dt <- dt[TraitID %in% trait_id]
  if ("ErrorRisk" %in% names(dt)) dt <- dt[is.na(ErrorRisk) | ErrorRisk <= 4]
  value_ok <- if (trait_name %in% c("growth_form", "woodiness")) {
    !is.na(dt$OrigValueStr) & dt$OrigValueStr != ""
  } else {
    !is.na(dt$StdValue)
  }
  dt <- as_tibble(dt[value_ok & !is.na(AccSpeciesName)])
  sp  <- used_species %>% filter(trait == trait_name) %>% pull(name)
  gen <- used_genera  %>% filter(trait == trait_name) %>% pull(genus)
  dt %>%
    filter(AccSpeciesName %in% sp | word(AccSpeciesName, 1) %in% gen) %>%
    transmute(trait = trait_name, DatasetID, Dataset, Reference,
              custodian = str_squish(paste(FirstName, LastName)))
}

records <- pmap_dfr(try_sources, function(trait, file, trait_id)
  read_used_records(trait, file, trait_id[[1]]))

# 4. One row per dataset ----------------------------------------------------------
try_refs <- records %>%
  group_by(DatasetID, Dataset) %>%
  summarise(
    reference = paste(sort(unique(na.omit(Reference[Reference != ""]))), collapse = " | "),
    custodian = paste(sort(unique(custodian)), collapse = "; "),
    traits    = paste(sort(unique(trait)), collapse = ", "),
    n_records = n(),
    .groups   = "drop"
  ) %>%
  arrange(Dataset)

write_csv(try_refs, out("try_references.csv"))

cat("TRY records used:", nrow(records), "| datasets:", nrow(try_refs), "\n")
print(count(records, trait, name = "n_records"))
cat("Datasets without a reference:", sum(try_refs$reference == ""), "\n")
cat("Written to", out("try_references.csv"), "\n")
