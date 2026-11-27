
# MORAclimshift

Code for: Torres, A., Bektas, B., Pradhan, K., Angert, A. & Hille Ris Lambers, J. *Climate debts in old-growth forests: contrasting community reassembly trajectories across forest layers and elevations*. Global Change Biology (submitted).

Authors: Agostina Torres<sup>1</sup>, Billur Bektas<sup>1</sup>, Kavya Pradhan<sup>2,3</sup>, Amy Angert<sup>4</sup>, Janneke Hille Ris Lambers<sup>1,2</sup>

<sup>1</sup> Institute of Integrative Biology, ETH Zürich, Zürich, Switzerland
<sup>2</sup> Department of Biology, University of Washington, Seattle, WA, USA
<sup>3</sup> Conservation International
<sup>4</sup> Biodiversity Research Centre and Department of Botany, University of British Columbia, Vancouver, British Columbia, Canada

Author of the code: Agostina Torres (ETH Zürich).

The pipeline analyses changes in overstory and understory plant communities in 324 permanent forest plots (316 for the understory) in Mount Rainier National Park, on Mt. Tahoma (WA, USA), first surveyed in about 1975 and resurveyed in 2014–2019. It computes climatic-affinity shifts (community temperature and moisture indices, CTI and CMI, from continent-wide species affinities) and their lags relative to climate change, and fits hierarchical joint species distribution models (HMSC) of species change with traits.

## Data

The data are archived in Dryad: https://doi.org/10.5061/dryad.p2ngf1w82. The archive contains:

- the original survey files (`raw_surveys.zip`), from which the cleaning scripts rebuild the vegetation tables;
- the processed tables used in the analyses (vegetation, plots and climate, species climatic affinities, traits, species pools and phylogeny), described column by column in the archive README (copy in `pipeline/07_data_archive/archive/README.md`).

To rerun the pipeline from the original survey files, unzip `raw_surveys.zip` and place the `_Org.csv` files in `pipeline/01_data_prep/01_cleaning/data/data_raw/01_raw_surveys/` and `MORA_PlotNames.csv` in `pipeline/01_data_prep/01_cleaning/data/data_raw/02_earlier_cleaning/`.

The other inputs come from public sources:

- Species occurrences: two GBIF downloads (GBIF.org, https://doi.org/10.15468/dl.cars7w and https://doi.org/10.15468/dl.ndtn3b), retrieved by the scripts in `pipeline/01_data_prep/05_species_distribution/01_gbif_occurrences/`. GBIF downloads require a GBIF account; the credentials are read from `.Renviron` (`GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL`).
- Climate: CHELSA v2.1 monthly time series and 1981–2010 bioclimatic variables, downloaded by the scripts in `pipeline/01_data_prep/03_climate/01_download_chelsa/` and `pipeline/01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/`.
- Traits: public TRY data (requests 38172, 40000, 52607 and 52732; https://www.try-db.org), placed in `pipeline/01_data_prep/04_traits/data_raw/`; GIFT and BIEN data are retrieved by the trait scripts through their R packages.

### Data cleaning

The cleaning code in `pipeline/01_data_prep/01_cleaning/` was written independently and starts from the original survey files. As a guide, I used the earlier cleaning of the same surveys by Kavya Pradhan (https://github.com/kavs-P/dissimilarityPNW), which drew on knowledge of the field surveys that is not recorded in the data files (e.g. how plots were named, which census to keep for plots surveyed twice, and known data-entry errors). Its plot-name table and documented corrections were used to check plot identifiers, coordinates and survey years.

## Requirements

R ≥ 4.3. Main packages: tidyverse, here, Hmsc, brms, lme4, terra, sf, rgbif, TNRS, WorldFlora, GIFT, BIEN, V.PhyloMaker2, ape, phytools, castor, caper, ggtree, hrbrthemes, patchwork, cowplot, circlize, viridisLite, data.table, officer and flextable. Versions of R and of all packages used in the scripts: `package_versions.txt`. All paths are relative to the project root through `here()`; open `MORAclimshift.Rproj` before running any script.

HMSC models were fitted on a 48-core Windows workstation; the full understory models with 100 reduced resamplings take several days.

## Structure and run order

Each module has a `code/` folder and writes to `data_clean/` (intermediate files), `outputs/` (tables) and `figures/`. Scripts are numbered in run order within each module; a script with a letter suffix (e.g. `01b`) runs after the script with the same number.

| Folder | What it does |
|---|---|
| `pipeline/01_data_prep/01_cleaning` | Plot identifiers, species names (TNRS), selection of the present-day census, overstory basal area and counts, understory cover and occurrence |
| `pipeline/01_data_prep/02_species_pool` | Shared and frequent species pools; trait and affinity coverage; species tables |
| `pipeline/01_data_prep/03_climate` | CHELSA monthly climate, annual and seasonal variables, past (1980–1983; 1979–1983 for precipitation) and present (2013–2017) windows, climate change along elevation |
| `pipeline/01_data_prep/04_traits` | Synonyms; traits from TRY, GIFT and BIEN; genus imputation and its validation; literature values; phylogeny; phylogenetic signal |
| `pipeline/01_data_prep/05_species_distribution` | GBIF download and cleaning (`01_gbif_occurrences`), then CHELSA bioclimatic values at the occurrences, species climatic affinities and their resamples, regional and ecoregion sensitivity (`02_integrate_gbif_chelsa`) |
| `pipeline/02_analysis/01_thermophilization/01_overstory`, `02_understory` | CTI and CMI per plot and survey, rates per decade, expected rates and lags (space-for-time), community tails, species and life-form contributions, replacement flows |
| `pipeline/02_analysis/01_thermophilization/03_integration_over_under` | Both layers together: Figure 2, indices along elevation, contributions, Bayesian rate models, tables, space-for-time relationships |
| `pipeline/02_analysis/02_hmsc` | Model inputs and pre-launch checks (01–01c); fitting and convergence (02–03); results and figures (04–05g); reduced resamplings of the past understory, pooling, robustness comparisons and variance partitioning (11–35) |
| `pipeline/04_robustness` | Robustness of the results to the difference in understory sampling between surveys (see its README) |
| `pipeline/05_manuscript_figures` | Figure 1 (conceptual panels) and Figure S1.1 (study area) |
| `pipeline/07_data_archive` | Export of the files archived in Dryad |

Exploratory scripts, not needed to reproduce the reported results: `pipeline/01_data_prep/06_exploring/`, `04_traits/code/07_impute_phylogenetic.R`, `02_hmsc/code/06_run_preliminary_medium.R` and `02_hmsc/code/08_plot_only_build_and_launch.R`.

## Reproducing a figure or table

Paths are relative to `pipeline/`.

| Item | Script |
|---|---|
| Figure 1 | `05_manuscript_figures/figure1/code/01_figure1_icons.R`, `02_figure1_panels.R` (panels assembled in PowerPoint) |
| Figure 2 | `02_analysis/01_thermophilization/03_integration_over_under/code/01_thermo_xero_lag_integration.R` |
| Figure 3 | `02_analysis/02_hmsc/code/05b_species_change_centered_figure.R` |
| Figure 4 | `02_analysis/02_hmsc/code/05_hmsc_figures_combined.R` (trait effects, Γ) |
| Figure S1.1 | `05_manuscript_figures/figureS1/code/figureS1_study_area.R` |
| Figures S1.2–S1.4 | `01_data_prep/03_climate/02_explore_climate/code/01_climate_trends.R`, `02_climate_change_spatial_structure.R`, `03_climate_pca_layout.R` |
| Table S2.1 | `01_data_prep/02_species_pool/code/04_species_pools_table.R` |
| Figure S2.1 | `01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/09_species_affinity_figure.R` |
| Figure S2.2 | `01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/08_affinity_correlation_figure.R` |
| Figures S2.3, S4.7 | `01_data_prep/05_species_distribution/02_integrate_gbif_chelsa/code/07_ecoregion_affinity_sensitivity.R` |
| Figure S2.4 | `01_data_prep/04_traits/code/09_trait_imputation.R` |
| Figure S2.5 | `01_data_prep/04_traits/code/12_phylogeny_traits_figure.R` |
| Table S2.2 | `01_data_prep/04_traits/code/11_phylogenetic_signal_traits.R` |
| Table S3.1 | `04_robustness/code/13_table_S1_understory_design_bias.R` |
| Figures S3.1, S3.3 | `04_robustness/code/08_diagnosis_figures.R` |
| Figure S3.2 | `04_robustness/code/05_rarefaction_occurrence.R` |
| Figure S4.1 | `04_robustness/code/12_rarefied_rate_models_bayes.R` |
| Figure S4.2 | `02_analysis/01_thermophilization/03_integration_over_under/code/03_cti_cmi_past_present_elevation_integration.R` |
| Figure S4.3 | `02_analysis/01_thermophilization/03_integration_over_under/code/06_rate_models.R` |
| Figure S4.4 | `02_analysis/01_thermophilization/03_integration_over_under/code/02_affinity_edges_integration.R` |
| Table S4.1 | `02_analysis/01_thermophilization/03_integration_over_under/code/14_table_climatic_affinity_shifts.R` |
| Figure S4.5 | `02_analysis/01_thermophilization/03_integration_over_under/code/04_species_lifeform_contribution_integration.R` |
| Table S4.2 | `02_analysis/01_thermophilization/03_integration_over_under/code/18_table_lag_linear_vs_loess.R` |
| Figure S4.6 | `02_analysis/01_thermophilization/03_integration_over_under/code/17_species_replacement_chord_figure.R` |
| Figure S4.8 | `02_analysis/01_thermophilization/03_integration_over_under/code/16_baseline_affinity_mismatch.R` |
| Figure S4.9 | `02_analysis/01_thermophilization/03_integration_over_under/code/15_space_for_time_and_affinity_correlations.R` |
| Figure S5.1 | `02_analysis/02_hmsc/code/18_convergence_figure.R` |
| Figure S5.2 | `02_analysis/02_hmsc/code/05_hmsc_figures_combined.R` |
| Figures S5.3, S5.4 | `02_analysis/02_hmsc/code/05f_hmsc_checks_figures.R` |
| Figure S5.5 | `02_analysis/02_hmsc/code/17_compare_spatial_oct09.R` |
| Figure S6.1 | `02_analysis/02_hmsc/code/35_variance_partitioning_figure.R` |
| Figure S6.2 | `02_analysis/02_hmsc/code/05d_species_change_elevation_heatmap.R` |
| Figure S6.3 | `02_analysis/02_hmsc/code/05c_trait_effects_elevation_figure.R` |
| Figure S6.4 | `02_analysis/02_hmsc/code/05g_species_change_range_margins.R` |
| Figure S6.5 | `02_analysis/02_hmsc/code/05e_residual_associations_figure.R` |

## Licence

Code: MIT licence (see `LICENSE`). Data: CC0 1.0, archived in Dryad (https://doi.org/10.5061/dryad.p2ngf1w82).
