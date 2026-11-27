# ------------------------------------------------------------------------------
# 06_build_phylogeny.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Purpose: build one phylogeny for all project species with occurrence
# (overstory + understory together) with V.PhyloMaker2, which grafts the
# species list onto GBOTB.extended.WP, a mega-phylogeny backbone of ~74,000
# seed plant species. The tree covers the full species catalog (broader
# than any analysis species pool); subsets for each layer or pool
# (pipeline/01_data_prep/02_species_pool) are obtained downstream by
# pruning tips, so the tree does not need to be rebuilt when a pool changes.
#
# Species names are matched with the same synonym fallback as in
# 02_build_trait_table_try.R. V.PhyloMaker2 places a species at its tip when
# its binomial is in the mega-tree, within its genus when not (grafted at a
# random point within the genus clade), or at its family stem when the genus
# is not recognized. Family-level placement carries much less phylogenetic
# information than genus-level placement, so a species whose TNRS-accepted
# name is not recognized by the backbone, but whose WFO synonym genus is,
# is submitted under that synonym:
#   pass 1: all species under their TNRS-accepted name.
#   pass 2: species placed only at family level in pass 1 are resubmitted
#           under the first WFO synonym (spp_synonyms.csv, from
#           01_build_species_synonyms.R) whose genus is in the mega-tree.
# Tip labels and C matrix dimnames are set back to the project's accepted
# names whatever name was submitted: the synonym is only used for matching.
#
# Requires V.PhyloMaker2 (GitHub-only, not on CRAN):
#   devtools::install_github("jinyizju/V.PhyloMaker2")
#
# V.PhyloMaker2 output names differ slightly between versions (e.g. the
# wording of the species-list status column, or whether the tree is
# $scenario.3 or $tree). str(tree_pass1, max.level = 1) and
# names(tree_pass1$species.list) show the names in the installed version;
# the two places that depend on them are marked below.
#
# Output:
#   - outputs/MORA_combined_phylo_tree.rda: C (phylogenetic correlation
#     matrix, cov2cor(ape::vcv(phy))) + phy (ape phylo object), covering all
#     project species with occurrence (overstory + understory together).
#   - tmp_files/phylomaker_placement_log.csv: per species, which name was
#     submitted (accepted vs. synonym) and the placement level achieved
#     (species/genus/family) -- family-level placements are the ones to
#     treat cautiously in any downstream phylogenetic trait imputation.
# ------------------------------------------------------------------------------

if (!requireNamespace("V.PhyloMaker2", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("jinyizju/V.PhyloMaker2")
}

library(V.PhyloMaker2)
library(tidyverse)
library(ape)
library(here)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)
cleaning_dir     <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)

dir.create(this_module_dir("outputs"),   recursive = TRUE, showWarnings = FALSE)
dir.create(this_module_dir("tmp_files"), recursive = TRUE, showWarnings = FALSE)

# 1. Combined species catalog (real occurrence, overstory + understory, union) -
# Same occurrence filter as 02_build_trait_table_try.R, but without that
# script's species_name_corrections step (the Berberis -> Mahonia-type
# understory renames applied before the climatic-affinity join) and without
# a species-pool restriction (species_pool_broad.csv /
# species_pool_strict.csv), which is applied downstream by pruning tips.
#
# species_name_corrections matches the names of understory_long.csv to the
# naming of the GBIF/CHELSA climatic-affinity table; spp_info.csv (the
# source of `family`, used below) uses the understory names as recorded
# (e.g. Berberis aquifolium). Applying the corrections here would leave the
# 6 renamed species (Mahonia aquifolium, Mahonia nervosa, Mycelis muralis,
# Ozomelis trifida, Hypopitys monotropa, Trautvetteria carolinensis)
# without a family, so family is looked up on the names recorded in
# understory_long.csv.
drop_zero_occurrence_species <- function(long_df) {
  long_df %>%
    group_by(species) %>%
    summarise(total_occ = sum(occurrence, na.rm = TRUE), .groups = "drop") %>%
    filter(total_occ > 0) %>%
    pull(species)
}

overstory_long  <- read_csv(cleaning_dir("outputs", "overstory_long.csv"),  show_col_types = FALSE)
understory_long <- read_csv(cleaning_dir("outputs", "understory_long.csv"), show_col_types = FALSE)

overstory_real_species  <- drop_zero_occurrence_species(overstory_long)
understory_real_species <- drop_zero_occurrence_species(understory_long)

all_species <- union(overstory_real_species, understory_real_species)
cat("Combined species pool (overstory + understory, real occurrence, deduplicated):",
    length(all_species), "\n")

# 2. Family lookup ---------------------------------------------------------------
# spp_info.csv has "species" and "family" columns (with life_form,
# tree_type, canopy_type, tree_resurvey, presence_coverage).
spp_info <- read_csv(cleaning_dir("outputs", "spp_info.csv"), show_col_types = FALSE)
cat("spp_info.csv columns:", paste(names(spp_info), collapse = ", "), "\n")

family_lookup <- spp_info %>%
  select(species, family) %>%
  distinct(species, .keep_all = TRUE)

sp_df <- tibble(species = all_species) %>%
  left_join(family_lookup, by = "species") %>%
  mutate(genus = word(species, 1))

n_no_family <- sum(is.na(sp_df$family))
cat("Species with no family match in spp_info.csv:", n_no_family,
    "-- these cannot be placed by V.PhyloMaker2 at all.\n")
if (n_no_family > 0) print(sp_df %>% filter(is.na(family)) %>% pull(species))

# 3. Pass 1: submit all species under their TNRS-accepted name -------------------
make_phylomaker_input <- function(df) {
  df %>% transmute(species = str_replace_all(species, " ", "_"), genus, family)
}

sp_list_pass1 <- make_phylomaker_input(sp_df %>% filter(!is.na(family)))

tree_pass1 <- phylo.maker(
  sp.list = sp_list_pass1,
  tree = GBOTB.extended.WP,
  nodes = nodes.info.1.WP,
  scenario = "S3"
)

placement_pass1 <- tree_pass1$species.list %>%
  as_tibble() %>%
  mutate(species = str_replace_all(species, "_", " "))

cat("\nPass 1 placement levels:\n")
print(count(placement_pass1, status, sort = TRUE))  # column name may differ between V.PhyloMaker2 versions

# 4. Pass 2: retry WFO synonyms for species that only got family-level binding --
spp_synonyms <- read_csv(this_module_dir("outputs", "spp_synonyms.csv"), show_col_types = FALSE)

# The string V.PhyloMaker2 uses for family-level placement depends on the
# version (e.g. "bind to family", "F"); see count(placement_pass1, status)
# above.
family_level_species <- placement_pass1 %>%
  filter(str_detect(str_to_lower(status), "famil")) %>%
  pull(species)

cat("\nSpecies with only family-level placement in pass 1:", length(family_level_species), "\n")

backbone_genera <- unique(word(GBOTB.extended.WP$tip.label, 1))

retry_candidates <- tibble(species = family_level_species) %>%
  inner_join(spp_synonyms, by = c("species" = "final_species")) %>%
  left_join(family_lookup, by = "species") %>%
  mutate(synonym_genus = word(synonym, 1)) %>%
  filter(synonym_genus %in% backbone_genera) %>%
  distinct(species, .keep_all = TRUE)  # first matching synonym per species

cat("Species rescued to genus-level via a WFO synonym:", nrow(retry_candidates), "\n")
if (nrow(retry_candidates) > 0) print(retry_candidates %>% select(species, synonym))

if (nrow(retry_candidates) > 0) {
  sp_list_pass2_additions <- retry_candidates %>%
    transmute(species = str_replace_all(synonym, " ", "_"), genus = synonym_genus, family)
  
  sp_list_final <- sp_list_pass1 %>%
    filter(!(str_replace_all(species, "_", " ") %in% retry_candidates$species)) %>%
    bind_rows(sp_list_pass2_additions)
} else {
  sp_list_final <- sp_list_pass1
}

# 5. Final tree build -------------------------------------------------------------
tree_final <- phylo.maker(
  sp.list = sp_list_final,
  tree = GBOTB.extended.WP,
  nodes = nodes.info.1.WP,
  scenario = "S3"
)

phy <- tree_final$scenario.3  # name of the result tree may differ between V.PhyloMaker2 versions

# Restore project accepted names as tip labels (undo synonym substitution and
# the underscore formatting V.PhyloMaker2 requires) -----------------------------
tip_rename <- tibble(tip_submitted = phy$tip.label) %>%
  mutate(tip_species = str_replace_all(tip_submitted, "_", " ")) %>%
  left_join(retry_candidates %>% select(species, synonym), by = c("tip_species" = "synonym")) %>%
  mutate(final_name = coalesce(species, tip_species))

phy$tip.label <- tip_rename$final_name

# 6. Phylogenetic correlation matrix (same construction as the existing
#    overstory/understory .rda files: correlation, not raw covariance) ---------
C <- cov2cor(ape::vcv(phy))

# 7. Save ---------------------------------------------------------------------------
save(C, phy, file = this_module_dir("outputs", "MORA_combined_phylo_tree.rda"))

placement_log <- tree_final$species.list %>%
  as_tibble() %>%
  mutate(tip_species = str_replace_all(species, "_", " ")) %>%
  left_join(tip_rename %>% select(tip_species, final_name), by = "tip_species") %>%
  select(submitted_name = tip_species, final_name, status)

write_csv(placement_log, this_module_dir("tmp_files", "phylomaker_placement_log.csv"))

cat("\n========================================\n")
cat("Combined phylogeny built:", length(phy$tip.label), "species\n")
print(count(placement_log, status, sort = TRUE))
cat("Saved: outputs/MORA_combined_phylo_tree.rda\n")
cat("Saved: tmp_files/phylomaker_placement_log.csv\n")
