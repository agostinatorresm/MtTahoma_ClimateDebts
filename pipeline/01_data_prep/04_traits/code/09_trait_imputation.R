# ------------------------------------------------------------------------------
# 09_trait_imputation.R
# Module: pipeline/01_data_prep/04_traits
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Trait tables, coverage and validation of genus-level imputation. Runs after
# 08_integration.R; 10_plot_traits.R (phylogeny figures) is independent.
#
# Resolution tiers for a species x trait value:
#   1. species level: TRY, GIFT, BIEN (vegetation name or synonym), or a
#      literature value (08);
#   2. genus level: mean/mode of all records of the genus of the vegetation
#      name (genera of synonyms only when that genus has no records);
#   3. otherwise "unresolved"; such species are dropped from the trait-based
#      pools in 02_species_pool.
#
# The leave-one-out validation (section 3) uses the same evidence pool as
# 05_impute_trait_by_genus.R: TRY + GIFT + BIEN records for seed mass and SLA,
# TRY for woodiness, and the same genus rule. Each species' own records,
# under its vegetation name or any synonym, are excluded before predicting it.
#
# Sections:
#   1. Species x trait long table and species count per source (all traits).
#   2. Trait coverage by source, layer and species pool (figure + CSV), for
#      the traits of the HMSC models of each layer.
#   3. Genus-level imputation leave-one-out validation (seed mass, SLA,
#      woodiness): RMSE (% of the observed mean) and R2 on log scale for the
#      numeric traits, accuracy for woodiness.
#   4. Combined figure: (A) coverage, (B) numeric validation.
#
# Tables -> outputs/ ; figures -> figures/.
# Figures: ggsave PNG, white background, 300 dpi, theme_ipsum without grid and
# with black axis lines (mora_theme), viridis scales.
# ------------------------------------------------------------------------------

library(tidyverse)
library(hrbrthemes)
library(ragg)
library(here)

this_module_dir <- function(...) here("pipeline", "01_data_prep", "04_traits", ...)
pool_dir         <- function(...) here("pipeline", "01_data_prep", "02_species_pool", ...)
outputs_dir      <- function(...) this_module_dir("outputs", ...)
figures_dir      <- function(...) this_module_dir("figures", ...)

dir.create(outputs_dir(), recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir(), recursive = TRUE, showWarnings = FALSE)

VIRIDIS_OPTION <- "D"
QUALITY_DPI    <- 300
GGSAVE_SCALE   <- 1.8

BASE_FONT_SIZE <- 20

mora_theme <- function(base_size = BASE_FONT_SIZE) {
  theme_ipsum(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title    = element_text(size = rel(1.3), face = "bold"),
      plot.subtitle = element_text(size = rel(1.05)),
      plot.caption  = element_text(size = rel(0.85), hjust = 0),
      # Axis titles smaller than the plot title; margins keep them clear of
      # the tick labels.
      axis.title    = element_text(size = rel(1.0), face = "bold"),
      axis.title.x  = element_text(size = rel(1.0), face = "bold", margin = margin(t = 10)),
      axis.title.y  = element_text(size = rel(1.0), face = "bold", margin = margin(r = 10)),
      axis.text     = element_text(size = rel(1.1)),
      strip.text    = element_text(size = rel(1.1), face = "bold"),
      legend.text   = element_text(size = rel(1.2)),
      legend.title  = element_text(size = rel(1.3), face = "bold")
    )
}

with_retry <- function(expr_fn, label, tries = 5, wait_sec = 2) {
  for (i in seq_len(tries)) {
    result <- try(expr_fn(), silent = TRUE)
    if (!inherits(result, "try-error")) return(invisible(result))
    if (i < tries) {
      cat("  (", label, ": write failed, attempt", i, "of", tries, "-- retrying in", wait_sec, "s)\n", sep = "")
      Sys.sleep(wait_sec)
    }
  }
  stop(label, ": could not write after ", tries, " attempts -- check the file isn't open ",
       "elsewhere and that Polybox isn't stuck syncing that folder.\n",
       "Underlying error: ", attr(result, "condition")$message)
}

save_png <- function(plot, path, width, height, bg = "white", scale = GGSAVE_SCALE) {
  with_retry(function() {
    ggsave(filename = path, plot = plot, device = ragg::agg_png, bg = bg,
           width = width, height = height, dpi = QUALITY_DPI, scale = scale,
           units = "in", limitsize = FALSE)
  }, label = basename(path))
}

write_csv_safe <- function(x, path, ...) {
  with_retry(function() write_csv(x, path, ...), label = basename(path))
}

ALL_NUMERIC_TRAITS     <- c("seed_mass", "longevity", "SLA", "wood_density")
ALL_CATEGORICAL_TRAITS <- c("growth_form", "woodiness")
ALL_TRAITS             <- c(ALL_NUMERIC_TRAITS, ALL_CATEGORICAL_TRAITS)

# Display set for every figure in this script: the traits used in the HMSC
# models. Wood density is shown for the overstory only and woodiness for the
# understory only (all overstory species are woody), as in the models.
PLOT_NUMERIC_TRAITS     <- c("seed_mass", "SLA")
PLOT_CATEGORICAL_TRAITS <- c("woodiness")
COVERAGE_TRAITS <- list(overstory  = c("seed_mass", "SLA", "wood_density"),
                        understory = c("seed_mass", "SLA", "woodiness"))

overstory_traits_raw  <- read_csv(outputs_dir("overstory_traits.csv"),  show_col_types = FALSE)
understory_traits_raw <- read_csv(outputs_dir("understory_traits.csv"), show_col_types = FALSE)

# Source class used to select observed (species-level) values for the
# leave-one-out validation. Literature values are kept out of the validation.
classify_src <- function(src) {
  case_when(
    is.na(src) ~ NA_character_,
    src == "unresolved" ~ "unresolved",
    str_starts(src, "genus_mean") | str_starts(src, "genus_mode") ~ "genus",
    str_starts(src, "phylogenetic") ~ "unresolved",
    str_starts(src, "literature:") ~ "literature",
    TRUE ~ "observed"
  )
}

# ==============================================================================
# 1. Species x trait long table + source summary table (FINAL deliverables,
#    cover every trait in ALL_TRAITS -- not display-restricted)
# ==============================================================================
cat("\n--- Section 1: species x trait table + summary (final_ prefix) ---\n")

TRAIT_UNITS <- c(
  seed_mass = "mg", longevity = "check_unit",  # longevity in the TRY StdValue unit
  SLA = "mm2/mg", wood_density = "g/cm3",
  growth_form = "categorical", woodiness = "categorical"
)

species_pool_broad  <- read_csv(pool_dir("data_clean", "species_pool_broad.csv"),  show_col_types = FALSE)
species_pool_strict <- read_csv(pool_dir("data_clean", "species_pool_strict.csv"), show_col_types = FALSE)

pool_flags <- bind_rows(
  species_pool_broad  %>% transmute(species, stratum, pool = "broad"),
  species_pool_strict %>% transmute(species, stratum, pool = "strict")
) %>%
  mutate(pool_label = paste0(stratum, "_", pool)) %>%
  # .keep_all keeps stratum and pool for section 2
  distinct(species, pool_label, .keep_all = TRUE)

pool_summary_per_species <- pool_flags %>%
  group_by(species) %>%
  summarise(pools = paste(sort(unique(pool_label)), collapse = ", "), .groups = "drop")

# Source label for the tables: the species-level databases that contributed
# to the value (e.g. "try+gift+bien"), "genus_imputation", "paper"
# (literature) or "unresolved". Conflicting sources that were not averaged
# are not counted ("+gift:conflict", "+bien:conflict").
classify_source <- function(src) {
  s <- str_remove_all(coalesce(src, ""), "\\+(gift|bien):conflict")
  has_try  <- !(str_starts(s, "gift_only") | str_starts(s, "bien:"))
  has_gift <- str_detect(s, "gift")
  has_bien <- str_detect(s, "bien")
  species_level <- pmap_chr(list(has_try, has_gift, has_bien), function(t, g, b) {
    paste(c("try", "gift", "bien")[c(t, g, b)], collapse = "+")
  })
  case_when(
    is.na(src) ~ NA_character_,
    src == "unresolved" ~ "unresolved",
    str_starts(src, "phylogenetic") ~ "unresolved",
    str_starts(src, "literature:") ~ "paper",
    str_starts(src, "genus_mean") | str_starts(src, "genus_mode") ~ "genus_imputation",
    TRUE ~ species_level
  )
}

build_long_table <- function(traits, stratum_label) {
  present_traits <- intersect(ALL_TRAITS, names(traits))
  values_long <- traits %>% select(species, all_of(present_traits)) %>%
    pivot_longer(cols = -species, names_to = "trait", values_to = "value", values_transform = list(value = as.character))
  src_long <- traits %>% select(species, ends_with("_src")) %>%
    pivot_longer(cols = -species, names_to = "trait_src", values_to = "src") %>%
    mutate(trait = str_remove(trait_src, "_src$")) %>% select(species, trait, src)
  
  values_long %>%
    left_join(src_long, by = c("species", "trait")) %>%
    mutate(stratum = stratum_label, source = classify_source(src), unit = TRAIT_UNITS[trait]) %>%
    left_join(pool_summary_per_species, by = "species") %>%
    mutate(pools = replace_na(pools, "none")) %>%
    select(species, pools, stratum, trait, source, value, unit)
}

species_trait_long <- bind_rows(
  build_long_table(overstory_traits_raw,  "overstory"),
  build_long_table(understory_traits_raw, "understory")
) %>% arrange(species, trait)

write_csv_safe(species_trait_long, outputs_dir("final_species_trait_table_long.csv"))
cat("Species x trait long table:", nrow(species_trait_long), "rows (",
    n_distinct(species_trait_long$species), "species)\n")

source_summary <- species_trait_long %>%
  count(trait, stratum, source) %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
  mutate(n_species = rowSums(across(where(is.numeric)))) %>%
  relocate(n_species, .after = stratum) %>%
  arrange(trait, stratum)

write_csv_safe(source_summary, outputs_dir("final_species_trait_source_summary.csv"))
cat("\nSpecies count per trait x stratum x source:\n")
print(source_summary, n = Inf)
cat("Saved: outputs/final_species_trait_table_long.csv, outputs/final_species_trait_source_summary.csv\n")

# ==============================================================================
# 2. Trait coverage by source (species-level / genus-imputed / unresolved),
#    by layer and species pool (broad vs. strict), for seed mass, SLA and
#    woodiness.
# ==============================================================================
cat("\n--- Section 2: trait coverage by source x pool x stratum ---\n")

source_simple_levels <- c("species-level", "genus-imputed", "unresolved")

pool_coverage_data <- pool_flags %>%
  inner_join(
    species_trait_long %>% filter(trait %in% unique(unlist(COVERAGE_TRAITS))) %>%
      select(species, trait, source),
    by = "species", relationship = "many-to-many"
  ) %>%
  filter(map2_lgl(stratum, trait, ~ .y %in% COVERAGE_TRAITS[[.x]])) %>%
  mutate(source_simple = case_when(
    is.na(source) ~ NA_character_,
    source == "genus_imputation" ~ "genus-imputed",
    source == "unresolved" ~ "unresolved",
    TRUE ~ "species-level"
  )) %>%
  count(pool, stratum, trait, source_simple) %>%
  mutate(source_simple = factor(source_simple, levels = source_simple_levels))

# CSV audit -- wide, one row per pool x layer x trait, same numbers as the
# figure below in table form.
pool_trait_coverage <- pool_coverage_data %>%
  pivot_wider(names_from = source_simple, values_from = n, values_fill = 0) %>%
  relocate(any_of(source_simple_levels), .after = trait) %>%
  arrange(pool, stratum, trait)

write_csv_safe(pool_trait_coverage, outputs_dir("species_pool_trait_coverage_summary.csv"))
cat("Species-pool trait coverage (broad + strict, display trait set):\n")
print(pool_trait_coverage, n = Inf)
cat("Saved: outputs/species_pool_trait_coverage_summary.csv\n")

# Proportions per pool; broad and strict pools drawn as side-by-side bars
# within each layer panel (manual x offset, since ggplot has no combined
# dodge + stack), broad transparent and strict opaque.
pool_coverage_prop <- pool_coverage_data %>%
  group_by(pool, stratum, trait) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup() %>%
  mutate(
    pool_label  = factor(pool, levels = c("broad", "strict"), labels = c("Broad", "Strict")),
    # position of the trait within its layer (three traits per layer)
    trait_pos = map2_int(stratum, trait, ~ match(.y, COVERAGE_TRAITS[[.x]])),
    trait_x   = trait_pos + if_else(pool == "broad", -0.21, 0.21)
  )
trait_names <- c(seed_mass = "Seed mass", SLA = "SLA", wood_density = "Wood density", woodiness = "Woodiness")
# x-axis labels per layer panel (ggplot cannot set different breaks per facet,
# so the labels are drawn as text under each panel)
coverage_axis_labels <- imap_dfr(COVERAGE_TRAITS, ~ tibble(stratum = .y, trait_x = seq_along(.x),
                                                           label = trait_names[.x]))

p_coverage <- ggplot(pool_coverage_prop, aes(x = trait_x, y = prop, fill = source_simple, alpha = pool_label)) +
  geom_col(position = "stack", width = 0.38) +
  facet_wrap(~ stratum, labeller = as_labeller(c(overstory = "Overstory", understory = "Understory"))) +
  geom_text(data = coverage_axis_labels, aes(x = trait_x, y = -0.05, label = label), inherit.aes = FALSE,
            vjust = 1, size = BASE_FONT_SIZE / 3.2) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_viridis_d(option = VIRIDIS_OPTION, name = "Trait source") +
  scale_alpha_manual(values = c(Broad = 0.35, Strict = 1), name = "Species pool",
                     labels = c(Broad = "Shared", Strict = "Frequent")) +
  # Pool legend keys drawn in the species-level colour, so transparency is visible
  guides(alpha = guide_legend(override.aes = list(fill = viridisLite::viridis(1, option = VIRIDIS_OPTION)))) +
  labs(x = NULL, y = "Species (%)") +
  mora_theme() +
  theme(axis.title.y = element_text(hjust = 0.5, margin = margin(r = 12)),
        panel.spacing.y = unit(1, "lines"), plot.margin = margin(10, 10, 30, 10),
        strip.text   = element_text(face = "bold", hjust = 0.5),
        panel.spacing.x = unit(2, "lines"),
        # the two legends on separate rows, so that they fit the figure width
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.text = element_text(size = BASE_FONT_SIZE * 0.85),
        legend.title = element_text(size = BASE_FONT_SIZE * 0.9, face = "bold"),
        legend.box.spacing = unit(1.6, "cm"))   # room for the trait labels under the bars

# ==============================================================================
# 3. Genus-level imputation: leave-one-out validation (seed mass, SLA,
#    woodiness). Each observed species is predicted from the records of its
#    genus set after removing its own records.
# ==============================================================================
cat("\n--- Section 3: genus-level imputation LOO validation ---\n")

extract_genus <- function(species_binomial) word(species_binomial, 1)

load_observed <- function(trait) {
  src_col <- paste0(trait, "_src")
  tables <- purrr::map(list(overstory = overstory_traits_raw, understory = understory_traits_raw), function(df) {
    if (!trait %in% names(df)) return(NULL)
    df %>% select(species, value = all_of(trait), src = all_of(src_col))
  })
  bind_rows(tables) %>%
    distinct(species, .keep_all = TRUE) %>%
    mutate(src_class = classify_src(src)) %>%
    filter(src_class == "observed", !is.na(value))
}

load_observed_for <- function(trait, numeric) {
  o <- load_observed(trait) %>% mutate(trait = trait)
  if (numeric) o <- o %>% mutate(observed = as.numeric(value)) else o <- o %>% mutate(observed = value)
  o %>% select(species, trait, observed)
}

observed_numeric_tbl     <- map_dfr(PLOT_NUMERIC_TRAITS,     load_observed_for, numeric = TRUE)
observed_categorical_tbl <- map_dfr(PLOT_CATEGORICAL_TRAITS, load_observed_for, numeric = FALSE)

raw_obj <- readRDS(this_module_dir("tmp_files", "trait_matching_intermediate.rds"))
list2env(raw_obj, envir = environment())

# Evidence pool identical to 05_impute_trait_by_genus.R
add_raw <- function(base, label, tr) {
  if (!exists(label)) return(base)
  extra <- get(label)[[tr]]
  if (is.null(extra)) return(base)
  bind_rows(base %>% select(species, value), extra %>% select(species, value))
}
numeric_raw_tables_all <- list(seed_mass = seed_mass_raw, SLA = sla_raw)
for (tr in names(numeric_raw_tables_all)) {
  for (label in c("gift_raw_by_trait", "bien_raw_by_trait")) {
    numeric_raw_tables_all[[tr]] <- add_raw(numeric_raw_tables_all[[tr]], label, tr)
  }
}
cat("LOO evidence pool: seed_mass/SLA from TRY",
    if (exists("gift_raw_by_trait")) "+ GIFT", if (exists("bien_raw_by_trait")) "+ BIEN",
    "; woodiness from TRY\n")
categorical_raw_tables <- list(woodiness = woodiness_raw)

# Names under which a species' own records can appear: its vegetation name,
# its synonyms, the name GIFT matched it through, and names in BIEN (BIEN
# records are matched through the same names). Own records are removed by
# this alias table before predicting a species.
species_aliases <- bind_rows(
  tibble(canonical = project_species$species, alias = project_species$species),  # a species is its own alias
  spp_synonyms %>% transmute(canonical = final_species, alias = synonym)
)

gift_species_level_path <- this_module_dir("tmp_files", "gift_species_level.csv")
if (file.exists(gift_species_level_path)) {
  gift_aliases <- read_csv(gift_species_level_path, show_col_types = FALSE) %>%
    distinct(match_species, species) %>%
    transmute(canonical = match_species, alias = species)
  species_aliases <- bind_rows(species_aliases, gift_aliases)
} else {
  cat("NOTE: tmp_files/gift_species_level.csv not found; GIFT match names not",
      "added to the alias table.\n")
}
species_aliases <- species_aliases %>% distinct(canonical, alias)

# Aliases shared by more than one project species are listed; exclusion
# applies to every species sharing the alias.
ambiguous_aliases <- species_aliases %>% count(alias, name = "n_canonical") %>% filter(n_canonical > 1)
if (nrow(ambiguous_aliases) > 0) {
  cat("NOTE:", nrow(ambiguous_aliases), "alias name(s) are shared by more than",
      "one canonical species (ambiguous synonym) -- LOO exclusion still",
      "applies to all matching canonical species for that alias:\n")
  print(species_aliases %>% semi_join(ambiguous_aliases, by = "alias") %>% arrange(alias))
}

# Candidate genera of each project species, as in 05: genus of the vegetation
# name (priority 1); genera of its synonyms only as a fallback (priority 2)
species_genera <- bind_rows(
  project_species %>% transmute(species, genus = extract_genus(species), priority = 1L),
  spp_synonyms %>% filter(!is.na(synonym)) %>%
    transmute(species = final_species, genus = extract_genus(synonym), priority = 2L)
) %>%
  semi_join(project_species, by = "species") %>%
  group_by(species, genus) %>% summarise(priority = min(priority), .groups = "drop")

# Records available to predict each observed species: records of its genus
# (or, if none remain, of its synonyms' genera) except those filed under any
# of its own names. A species without own records is predicted from the full
# genus, as in production. Species with no other records are not scored.
loo_pool <- function(raw_df, trait_label, observed_tbl) {
  obs <- observed_tbl %>% filter(trait == trait_label) %>% select(species, observed)
  own <- species_aliases %>% semi_join(obs, by = c("canonical" = "species"))
  obs %>%
    inner_join(species_genera, by = "species") %>%
    inner_join(raw_df %>% transmute(rec_species = species, value, genus = extract_genus(species)),
               by = "genus", relationship = "many-to-many") %>%
    anti_join(own, by = c("species" = "canonical", "rec_species" = "alias")) %>%
    group_by(species) %>%
    slice_min(priority, with_ties = TRUE) %>%
    ungroup()
}

genus_loo_numeric_points <- function(raw_df, trait_label, observed_tbl) {
  loo_pool(raw_df, trait_label, observed_tbl) %>%
    group_by(species, observed) %>%
    summarise(predicted = mean(value, na.rm = TRUE), .groups = "drop") %>%
    transmute(species, trait = trait_label, observed, predicted)
}

genus_loo_categorical_points <- function(raw_df, trait_label, observed_tbl) {
  loo_pool(raw_df, trait_label, observed_tbl) %>%
    count(species, observed, value, name = "n") %>%
    group_by(species, observed) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(species, trait = trait_label, observed, predicted = value)
}

genus_points_numeric <- map_dfr(PLOT_NUMERIC_TRAITS, function(tr) {
  genus_loo_numeric_points(numeric_raw_tables_all[[tr]], tr, observed_numeric_tbl)
})
genus_points_categorical <- map_dfr(PLOT_CATEGORICAL_TRAITS, function(tr) {
  genus_loo_categorical_points(categorical_raw_tables[[tr]], tr, observed_categorical_tbl)
})

# RMSE as % of the observed mean, so seed mass and SLA are comparable.
# r_squared: squared Pearson correlation of observed vs predicted.
# r_squared_log: the same on log values; seed mass is strongly right-skewed and
# the log scale weights each species by its relative deviation.
numeric_metrics <- genus_points_numeric %>%
  group_by(trait) %>%
  summarise(n = n(),
            rmse = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
            mean_observed = mean(observed, na.rm = TRUE),
            rmse_pct = round(100 * rmse / mean_observed, 1),
            r_squared = round(cor(observed, predicted, use = "complete.obs")^2, 3),
            r_squared_log = round(cor(log(observed[observed > 0 & predicted > 0]),
                                      log(predicted[observed > 0 & predicted > 0]),
                                      use = "complete.obs")^2, 3),
            # Label position in data units (top-left of each facet)
            label_x = min(observed, na.rm = TRUE),
            label_y = max(predicted, na.rm = TRUE),
            .groups = "drop")

write_csv_safe(numeric_metrics, outputs_dir("genus_imputation_validation_numeric.csv"))
cat("\nGenus-level imputation LOO -- numeric traits:\n"); print(numeric_metrics)

categorical_metrics <- genus_points_categorical %>%
  group_by(trait) %>%
  summarise(n = n(), accuracy_pct = round(100 * mean(predicted == observed), 1), .groups = "drop")

write_csv_safe(categorical_metrics, outputs_dir("genus_imputation_validation_categorical.csv"))
cat("\nGenus-level imputation LOO -- categorical traits:\n"); print(categorical_metrics)

# Log-log axes show the full range of seed mass without a display cap and
# match r_squared_log; the dashed line is the 1:1 line.
trait_facet_labels <- c(seed_mass = "Seed mass", SLA = "SLA")

p_numeric <- ggplot(genus_points_numeric, aes(x = observed, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.7, size = 2, color = "black") +
  scale_x_log10() +
  scale_y_log10() +
  geom_text(data = numeric_metrics,
            aes(x = label_x, y = label_y, label = paste0("R² (log) = ", r_squared_log, "\n(n = ", n, ")")),
            hjust = 0, vjust = 1, inherit.aes = FALSE, size = 5.5) +
  facet_wrap(~ trait, scales = "free", labeller = as_labeller(trait_facet_labels)) +
  labs(x = "Observed value (log scale)", y = "Predicted value (log scale)") +
  mora_theme() +
  theme(axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        strip.text   = element_text(face = "bold", hjust = 0.5),
        panel.spacing.x = unit(2, "lines"))

p_categorical <- ggplot(categorical_metrics, aes(x = trait, y = accuracy_pct, fill = trait)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = paste0(accuracy_pct, "%\n(n=", n, ")")), vjust = -0.3, size = 5.5) +
  scale_fill_viridis_d(option = VIRIDIS_OPTION, guide = "none") +
  labs(title = "Genus-level imputation validation -- woodiness", x = NULL, y = "Accuracy (%)") +
  ylim(0, 105) + mora_theme() +
  theme(axis.text.x = element_text(size = rel(1.2), face = "bold"))

# ==============================================================================
# 4. Combined figure: (A) coverage, (B) numeric LOO validation. Woodiness
#    accuracy is reported in its CSV (p_categorical is built but not combined).
# ==============================================================================
cat("\n--- Section 4: combined coverage + validation figure ---\n")

if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
library(ggpubr)

p_coverage_and_validation <- ggarrange(
  p_coverage + theme(plot.margin = margin(10, 10, 35, 10)),   # bottom margin for the trait labels
  NULL,
  p_numeric  + theme(plot.margin = margin(10, 10, 10, 10)),
  nrow = 3, heights = c(1.15, 0.06, 1.1),
  labels = c("A", "", "B"), font.label = list(size = BASE_FONT_SIZE * 1.2, face = "bold")
)

# scale = 1: larger ggsave scale values shrink text.
out_path <- figures_dir("trait_coverage_and_imputation_validation.png")
save_png(p_coverage_and_validation, out_path, width = 13, height = 14, scale = 1)
cat("Saved:", out_path, "\n")

open_fig <- function(f) if (.Platform$OS.type == "windows") shell.exec(f) else system2("open", shQuote(f))
open_fig(out_path)

cat("\n========================================\n")
cat("09_trait_imputation.R COMPLETE\n")
cat("Tables  -> outputs/final_species_trait_table_long.csv, outputs/final_species_trait_source_summary.csv,\n")
cat("           outputs/species_pool_trait_coverage_summary.csv,\n")
cat("           outputs/genus_imputation_validation_numeric.csv, outputs/genus_imputation_validation_categorical.csv\n")
cat("Figure  -> figures/trait_coverage_and_imputation_validation.png\n")
cat("========================================\n")