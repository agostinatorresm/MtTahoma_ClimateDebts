# ------------------------------------------------------------------------------
# 03_species_table_si.R
# Module: pipeline/01_data_prep/02_species_pool
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Supporting Information S2 tables of the species pools.
#
# Table A (summary, per stratum): number of species after cleaning, in the
#   shared and frequent pools, and with climatic-affinity and functional-trait
#   coverage within each pool; the two analysis sets are
#     thermophilization: shared pool with climatic affinity
#     HMSC:              frequent pool with climatic affinity and functional traits
#
# Table B (one row per species after cleaning, one table per stratum):
#   - plots occupied at each time;
#   - pool membership (shared, frequent);
#   - climatic affinity: available (✓) or not (–);
#   - functional traits used in the HMSC models of the layer, with the
#     source of each value:
#       O = species-level record (TRY / GIFT / BIEN, including synonyms)
#       G = imputed from congeners (genus mean or mode)
#       L = literature value
#       – = missing
#   - analysis sets: thermophilization and HMSC;
#   - highlights (HMSC species): species retained in the models whose
#     estimates are interpreted with caution, with the reason:
#       documented data issues (manual flags of 01c_prelaunch_checks.R:
#         identification, non-climatic change, recruitment pulse, species
#         complex) and, for the overstory, few abundance records;
#       understory, from the past vs present design diagnostics (04_robustness):
#         rare: fewer than 10 plots on average across times (occurrence not
#           interpreted individually);
#         one or two plots with cover at both times (cover not interpreted
#           individually); no paired cover (occurrence only);
#         sparse: past occurrence increases with plot area;
#         patchy: mean present subplot occupancy below 0.4.
#   Coverage is assessed for shared-pool species only (02_full_coverage_species.R);
#   the coverage cells of species outside the shared pool are left empty.
#
# Inputs : data_clean/species_pool_candidates_<stratum>.csv          (01)
#          data_clean/<stratum>_trait_climatic_coverage.csv          (02)
#          outputs/<stratum>_species_pool_broad_full_clim_cov.csv    (02)
#          outputs/<stratum>_species_pool_strict_full_func_clim_cov.csv (02)
#          04_traits/outputs/<stratum>_traits.csv  (trait sources)
#          02_analysis/02_hmsc/outputs/<pool>_species_flags.csv (01c_prelaunch_checks.R)
#          04_robustness/outputs/occurrence/area_sensitivity_species.csv (03)
#          04_robustness/outputs/accumulation_curve/accumulation_curve_summary.csv (01b)
#          04_robustness/outputs/census_only_cover/census_only_species.csv (07)
# Outputs (outputs/):
#   species_pool_summary_si.csv, species_table_si.csv, species_table_si.docx
# ------------------------------------------------------------------------------

library(here)
library(tidyverse)
library(flextable)
library(officer)
# flextable functions are called with the flextable:: prefix because other
# packages (e.g. terra, ggpubr) mask some of them when loaded in the same
# session

this_module_dir <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)
traits_dir   <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)

# ------------------------------------------------------------------------------
# 0. Settings
# ------------------------------------------------------------------------------
strata <- c("overstory", "understory")

# Functional traits shown per layer (as in 02_full_coverage_species.R)
functional_traits <- list(
  overstory  = c("seed_mass", "wood_density"),
  understory = c("SLA", "seed_mass", "woodiness")
)
trait_labels <- c(seed_mass = "Seed mass", wood_density = "Wood density",
                  SLA = "SLA", woodiness = "Woodiness")

# Species flags written by 01c_prelaunch_checks.R
flags_file <- here("pipeline", "02_analysis", "02_hmsc", "outputs", "strict_species_flags.csv")
# Flags kept from 01c: documented data issues and, for the overstory, few
# abundance records. The automatic understory flags of 01c are replaced by the
# design diagnostics of 04_robustness (see Highlights).
flag_labels <- c(identification        = "Identification inconsistent between times",
                 non_climatic          = "Change driven by a non-climatic factor",
                 stand_dynamics        = "Recruitment pulse",
                 complex               = "Species complex (congeners merged)",
                 few_abundance_records = "Few abundance records")
robustness_dir <- function(...) here("pipeline", "04_robustness", "outputs", ...)
rare_mean_plots <- 10   # rare: fewer plots than this on average across times

# Names in the trait tables -> names in the pools (as in 02_full_coverage_species.R)
name_to_pool <- c("Mahonia nervosa"            = "Berberis nervosa",
                  "Mahonia aquifolium"         = "Berberis aquifolium",
                  "Trautvetteria carolinensis" = "Trautvetteria caroliniensis",
                  "Hypopitys monotropa"        = "Monotropa hypopitys",
                  "Mycelis muralis"            = "Lactuca muralis",
                  "Ozomelis trifida"           = "Mitella trifida",
                  "Xanthocyparis nootkatensis" = "Callitropsis nootkatensis")

font_name <- "Arial Narrow"
font_size <- 8

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
read_module <- function(...) {
  f <- this_module_dir(...)
  if (!file.exists(f)) stop("File not found (run 01 and 02 of this module first): ", f)
  read_csv(f, show_col_types = FALSE)
}
read_species_col <- function(...) {
  df <- read_module(...)
  if ("species" %in% names(df)) df$species else df[[1]]
}

# Trait source code from the <trait>_src provenance tag
source_code <- function(value, src) {
  case_when(is.na(value)                                             ~ "\u2013",
            is.na(src)                                               ~ "O",
            str_detect(src, regex("unresolved", ignore_case = TRUE)) ~ "\u2013",
            str_detect(src, regex("genus", ignore_case = TRUE))      ~ "G",
            str_detect(src, regex("lit", ignore_case = TRUE))        ~ "L",
            TRUE                                                     ~ "O")
}
tick <- function(x) case_when(is.na(x) ~ "", x ~ "\u2713", TRUE ~ "\u2013")

# ------------------------------------------------------------------------------
# 1. Species table per layer
# ------------------------------------------------------------------------------
build_stratum_table <- function(s) {
  cand <- read_module("data_clean", paste0("species_pool_candidates_", s, ".csv"))
  cov  <- read_module("data_clean", paste0(s, "_trait_climatic_coverage.csv"))
  thermo_sp <- read_species_col("outputs", paste0(s, "_species_pool_broad_full_clim_cov.csv"))
  hmsc_sp   <- read_species_col("outputs", paste0(s, "_species_pool_strict_full_func_clim_cov.csv"))
  
  req <- functional_traits[[s]]
  tr <- read_csv(traits_dir("outputs", paste0(s, "_traits.csv")), show_col_types = FALSE) %>%
    mutate(species = recode(species, !!!name_to_pool)) %>%
    distinct(species, .keep_all = TRUE)
  codes <- tibble(species = cov$species) %>% left_join(tr, by = "species")
  spore_genera <- c("Adiantum", "Asplenium", "Athyrium", "Blechnum", "Cystopteris",
                    "Dryopteris", "Equisetum", "Gymnocarpium", "Lycopodium",
                    "Polypodium", "Polystichum", "Pteridium")
  if ("seed_mass" %in% names(codes)) {
    codes$seed_mass[word(codes$species, 1) %in% spore_genera] <- NA
  }
  codes <- bind_cols(
    tibble(species = codes$species),
    map_dfc(set_names(req), function(t) {
      src <- if (paste0(t, "_src") %in% names(codes)) codes[[paste0(t, "_src")]] else NA_character_
      source_code(codes[[t]], src)
    })
  )
  
  out <- cand %>%
    transmute(stratum = s, species, n_plots_past, n_plots_present,
              shared = included_broad, frequent = included_strict) %>%
    left_join(cov %>% select(species, climatic = climatic_ok, functional = functional_ok),
              by = "species") %>%
    left_join(codes, by = "species") %>%
    mutate(thermophilization = species %in% thermo_sp,
           hmsc              = species %in% hmsc_sp)
  
  # Consistency between the coverage table and the lists read by the analyses
  check <- out %>%
    filter(thermophilization != (shared & coalesce(climatic, FALSE)) |
             hmsc != (frequent & coalesce(climatic & functional, FALSE)))
  if (nrow(check) > 0) {
    print(check)
    stop(s, ": analysis lists and coverage table disagree (re-run 02_full_coverage_species.R).")
  }
  out
}

species_tables <- map(set_names(strata), build_stratum_table)

# ------------------------------------------------------------------------------
# 2. Highlights
# ------------------------------------------------------------------------------
read_robustness <- function(...) {
  f <- robustness_dir(...)
  if (!file.exists(f)) stop("File not found (run 04_robustness first): ", f)
  read_csv(f, show_col_types = FALSE) %>% mutate(species = str_replace_all(species, "_", " "))
}

# Documented data issues (01c), overstory few abundance records
flags_01c <- if (file.exists(flags_file)) {
  read_csv(flags_file, show_col_types = FALSE) %>%
    mutate(species = str_replace_all(species, "_", " ")) %>%
    filter(flag %in% names(flag_labels),
           flag != "few_abundance_records" | stratum == "overstory") %>%
    transmute(stratum, species, label = flag_labels[flag], order = match(flag, names(flag_labels)))
} else {
  warning("Species flags not found (run 01c_prelaunch_checks.R): ", flags_file)
  tibble(stratum = character(), species = character(), label = character(), order = integer())
}

# Understory design diagnostics, for the HMSC species
hmsc_under <- species_tables$understory %>% filter(hmsc)
paired <- read_robustness("census_only_cover", "census_only_species.csv") %>%
  select(species, n_plots_paired)
flags_design <- bind_rows(
  hmsc_under %>% filter((n_plots_past + n_plots_present) / 2 < rare_mean_plots) %>%
    transmute(species, label = "Rare (occurrence not interpreted individually)", order = 10),
  hmsc_under %>% inner_join(paired, by = "species") %>%
    filter(n_plots_paired %in% 1:2) %>%
    transmute(species, label = "1\u20132 plots with paired cover (cover not interpreted individually)", order = 11),
  hmsc_under %>% inner_join(paired, by = "species") %>%
    filter(n_plots_paired == 0) %>%
    transmute(species, label = "No paired cover (occurrence only)", order = 12),
  read_robustness("occurrence", "area_sensitivity_species.csv") %>% filter(area_sensitive) %>%
    semi_join(hmsc_under, by = "species") %>%
    transmute(species, label = "Sparse", order = 13),
  read_robustness("accumulation_curve", "accumulation_curve_summary.csv") %>% filter(patchy) %>%
    semi_join(hmsc_under, by = "species") %>%
    transmute(species, label = "Patchy in present quadrats", order = 14)
) %>% mutate(stratum = "understory")

highlights <- bind_rows(flags_01c, flags_design) %>%
  arrange(stratum, species, order) %>%
  group_by(stratum, species) %>%
  summarise(highlights = paste(unique(label), collapse = "; "), .groups = "drop")
cat("Species with highlights:\n")
print(bind_rows(flags_01c, flags_design) %>% count(stratum, label), n = Inf)

species_tables <- map(species_tables, ~ .x %>%
                        left_join(highlights, by = c("stratum", "species")) %>%
                        mutate(highlights = replace_na(highlights, "")) %>%
                        arrange(desc(shared), desc(frequent), species))

flags_outside_hmsc <- bind_rows(species_tables) %>% filter(highlights != "", !hmsc)
if (nrow(flags_outside_hmsc) > 0) {
  message("Flagged species not in the HMSC set: ",
          paste(flags_outside_hmsc$species, collapse = ", "))
}

write_excel_csv(bind_rows(species_tables), this_module_dir("outputs", "species_table_si.csv"))

# ------------------------------------------------------------------------------
# 3. Summary table
# ------------------------------------------------------------------------------
pool_summary <- imap_dfr(species_tables, ~ {
  sh <- filter(.x, shared); fr <- filter(.x, frequent)
  tibble(
    Layer = str_to_sentence(.y),
    `After cleaning` = nrow(.x),
    `Shared: all` = nrow(sh),
    `Shared: climatic affinity` = sum(sh$climatic),
    `Shared: functional traits` = sum(sh$functional),
    `Shared: climatic + functional` = sum(sh$climatic & sh$functional),
    `Frequent: all` = nrow(fr),
    `Frequent: climatic affinity` = sum(fr$climatic),
    `Frequent: functional traits` = sum(fr$functional),
    `Frequent: climatic + functional` = sum(fr$climatic & fr$functional)
  )
})
write_csv(pool_summary, this_module_dir("outputs", "species_pool_summary_si.csv"))
print(pool_summary)

style_ft <- function(ft) {
  ft %>%
    flextable::bold(part = "header") %>%
    flextable::font(fontname = font_name, part = "all") %>%
    flextable::fontsize(size = font_size, part = "all") %>%
    flextable::padding(padding.top = 1, padding.bottom = 1, part = "all") %>%
    flextable::border_remove() %>%
    flextable::hline_top(border = fp_border(width = 1), part = "header") %>%
    flextable::hline_bottom(border = fp_border(width = 0.75), part = "header") %>%
    flextable::hline_bottom(border = fp_border(width = 1), part = "body") %>%
    flextable::fontsize(size = font_size - 1, part = "footer")
}

ft_summary <- flextable::flextable(pool_summary) %>%
  flextable::set_header_labels(values = setNames(as.list(str_to_sentence(str_remove(names(pool_summary), "^(Shared|Frequent): "))),
                                                 names(pool_summary))) %>%
  flextable::add_header_row(values = c("", "", "Shared pool", "Frequent pool"), colwidths = c(1, 1, 4, 4)) %>%
  flextable::align(j = -1, align = "center", part = "all") %>%
  flextable::width(j = 1, width = 0.9) %>%
  flextable::width(j = -1, width = 0.8) %>%
  flextable::add_footer_lines(paste0(
    "After cleaning: species recorded in plots sampled at both times. Shared pool: ",
    "species recorded at both times. Frequent pool: shared species recorded in more than ",
    "five plots at each time. Climatic affinity: mean annual temperature and precipitation ",
    "across the species' occurrences. Functional traits: traits in the HMSC trait matrix ",
    "(overstory: seed mass, wood density; understory: SLA, seed mass, woodiness; seed mass is ",
    "not defined for spore-bearing plants). Climatic-affinity shift ",
    "analyses use shared species with climatic affinity; HMSC models use the restricted pool, i.e. frequent species ",
    "with climatic affinity and functional traits (12 overstory and 58 understory species; understory: 21 woody and 37 non-woody).")) %>%
  style_ft()

# ------------------------------------------------------------------------------
# 4. Species tables (one per stratum)
# ------------------------------------------------------------------------------
make_species_ft <- function(tab, s) {
  req <- functional_traits[[s]]
  display <- tab %>%
    mutate(climatic = if_else(shared, tick(climatic), ""),
           across(c(shared, frequent, thermophilization, hmsc), tick),
           across(all_of(req), ~ replace_na(.x, ""))) %>%
    select(species, n_plots_past, n_plots_present, shared, frequent,
           climatic, all_of(req), thermophilization, hmsc, highlights)
  
  header_top <- c("", "Plots", "Plots", "Species pool", "Species pool", "Coverage",
                  rep("Coverage", length(req)), "Analysis", "Analysis", "")
  header_bottom <- c(species = "Species", n_plots_past = "~1975", n_plots_present = "~2015",
                     shared = "Shared", frequent = "Frequent", climatic = "Climatic affinity",
                     trait_labels[req], thermophilization = "Thermo-philization",
                     hmsc = "HMSC", highlights = "Highlights")
  
  center_cols <- setdiff(names(display), c("species", "highlights"))
  flextable::flextable(display) %>%
    flextable::set_header_labels(values = as.list(header_bottom)) %>%
    flextable::add_header_row(values = header_top, top = TRUE) %>%
    flextable::merge_h(part = "header") %>%
    flextable::italic(j = "species", part = "body") %>%
    flextable::align(j = center_cols, align = "center", part = "all") %>%
    flextable::width(j = "species", width = 1.9) %>%
    flextable::width(j = center_cols, width = 0.65) %>%
    flextable::width(j = "highlights", width = 2.4) %>%
    flextable::add_footer_lines(paste0(
      "Species after cleaning, ordered by pool. Coverage was assessed for shared-pool species. ",
      "Functional-trait sources: O, species-level record (TRY/GIFT/BIEN); G, imputed from congeners; ",
      "L, literature value; \u2013, missing or not applicable (seed mass of spore-bearing plants). Highlights: species retained in all models whose ",
      "estimates are interpreted with caution (Supporting Information S2, S3).")) %>%
    style_ft()
}

doc <- read_docx() %>%
  body_add_par("Species pools: number of species", style = "heading 2") %>%
  body_add_flextable(ft_summary)
for (s in strata) {
  doc <- doc %>%
    body_add_break() %>%
    body_add_par(paste0("Species pools: ", s), style = "heading 2") %>%
    body_add_flextable(make_species_ft(species_tables[[s]], s))
}
doc <- doc %>% body_set_default_section(prop_section(page_size(orient = "landscape")))
print(doc, target = this_module_dir("outputs", "species_table_si.docx"))
