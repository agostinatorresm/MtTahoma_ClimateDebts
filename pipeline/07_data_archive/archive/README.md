
# Data for: *Climate debts in old-growth forests: contrasting community reassembly trajectories across forest layers and elevations*

Agostina Torres (1), Billur Bektas (1), Kavya Pradhan (2, 3), Amy Angert (4), Janneke Hille Ris Lambers (1, 2)

(1) Institute of Integrative Biology, ETH Zürich, Zürich, Switzerland

(2) Department of Biology, University of Washington, Seattle, WA, USA

(3) Conservation International

(4) Biodiversity Research Centre and Department of Botany, University of British Columbia, Vancouver, British Columbia, Canada

Global Change Biology (submitted). Contact: Agostina Torres (ETH Zürich).

Code: the analysis code is archived with this dataset (Software, Zenodo) and is also available at https://github.com/agostinatorresm/MtTahoma_ClimateDebts.

Vegetation of 324 forest plots (316 for the understory) in Mount Rainier National Park, on Mt. Tahoma (WA, USA), first surveyed in about 1975 and resurveyed in 2014–2019, with the species climatic affinities, traits and phylogeny used in the analyses. Except for `raw_surveys.zip`, all files are the data after cleaning, as used in the analyses; the cleaning scripts are part of the code. Plot identifiers (`plot_id`) and species names (`species`) link all files.

## vegetation_overstory.csv

One row per plot × survey × tree species, for plots surveyed in both periods and species recorded at least once.

| Column | Description | Units |
|---|---|---|
| plot_id | Plot identifier | |
| period | Survey: `past` (about 1975) or `present` (2014–2019) | |
| survey_year | Year of the survey of that plot; for plots with two present-day censuses averaged, the mean of the two years | year |
| species | Accepted species name | |
| occurrence | 1 if the species was recorded in the plot, 0 otherwise | |
| stem_count | Number of stems (mean of the two censuses where averaged) | stems |
| basal_area_m2 | Basal area of the species in the plot, from stem diameters grouped in 10-cm classes (class mid-points), in both surveys | m² |
| plot_area_m2 | Plot area (past: as recorded; present: 500 m²) | m² |

## vegetation_understory.csv

One row per plot × survey × understory species.

| Column | Description | Units |
|---|---|---|
| plot_id, period, survey_year, species | As above | |
| occurrence | 1 if recorded in the plot (present: in any of the four 1-m² quadrats or the 50-m² subplot), 0 otherwise | |
| cover_percent | Cover (past: plot-level estimate; present: mean over the four 1-m² quadrats, with zeros where absent). NA where the species was recorded only in the 50-m² subplot | % |

## plots.csv

One row per plot.

| Column | Description | Units |
|---|---|---|
| plot_id | Plot identifier | |
| latitude, longitude | Plot location (WGS84) | decimal degrees |
| elevation_m | Elevation | m |
| aspect_deg, slope_deg | Aspect and slope | degrees |
| plot_area_past_m2 | Plot area in the past survey | m² |
| year_past | Year of the past survey | year |
| year_present_overstory, year_present_understory | Year of the present-day survey used for each layer | year |
| selection_present_overstory, selection_present_understory | How the present-day census was chosen: single census, documented year, complete census, or two censuses averaged | |
| MAT_past, MAT_present | Mean annual temperature, CHELSA, mean of 1980–1983 (past) and 2013–2017 (present) | °C |
| MAT_min_past, MAT_min_present | Mean annual minimum temperature, same windows | °C |
| MAP_past, MAP_present | Mean annual precipitation, CHELSA, mean of 1979–1983 (past) and 2013–2017 (present) | mm |
| PAS_past, PAS_present | Precipitation as snow (precipitation in months with mean temperature ≤ 0 °C), 1980–1983 (past) and 2013–2017 (present) | mm |
| n_years_* | Number of years with data in each window | years |

## species_climatic_affinity.csv

One row per species. Affinities from GBIF occurrence records in North America up to 1985 (coordinate uncertainty < 100 km, one record per species per ~1-km cell) and CHELSA 1981–2010 climate.

| Column | Description | Units |
|---|---|---|
| species | Accepted species name | |
| n_occurrences | Records after cleaning | |
| n_truncated_MAT, n_truncated_MAP | Records within the 5–95 % range of each variable | |
| MAT_mean | Thermal affinity: mean annual temperature of records within the 5–95 % range | °C |
| MAP_mean | Moisture affinity: mean annual precipitation of records within the 5–95 % range | mm |

## species_climatic_affinity_resamples.csv

100 resamples per species of the occurrence climates within the 5–95 % range, used to compute community indices.

| Column | Description | Units |
|---|---|---|
| species | Accepted species name | |
| rep | Resample number (1–100) | |
| MAT_resample, MAP_resample | Mean annual temperature and precipitation of a resampled occurrence | °C, mm |

## species_traits.csv

One row per layer × species.

| Column | Description | Units |
|---|---|---|
| layer | `overstory` or `understory` | |
| species | Accepted species name | |
| seed_mass | Seed mass | mg |
| SLA | Specific leaf area | mm² mg⁻¹ |
| wood_density | Wood density (overstory only) | g cm⁻³ |
| woodiness | `woody` or `non-woody` | |
| growth_form | tree, shrub, vine, herb or grass | |
| longevity | Plant longevity | years |
| *_src | Source of the value: species-level records (TRY, GIFT, BIEN), genus mean or mode (`genus_mean_…:Genus`), literature, or `unresolved` | |
| *_n_records | Number of records the value is based on | |

Seed mass is undefined for spore-bearing plants and is missing for them.

## try_references.csv

One row per TRY dataset that contributed records to the trait values in `species_traits.csv` (species-level values, and genus means or modes through the records of the genus). Public TRY data (Kattge et al. 2020) are released under CC BY 4.0, which requires citing TRY and the original datasets.

| Column | Description |
|---|---|
| DatasetID | TRY dataset identifier |
| Dataset | Dataset name |
| reference | Reference(s) of the dataset, separated by " \| " |
| custodian | Dataset custodian(s) |
| traits | Traits to which the dataset contributed records |
| n_records | Number of records used |

## species_pools.csv

One row per layer × species recorded after cleaning.

| Column | Description |
|---|---|
| layer, species | As above |
| n_plots_past, n_plots_present | Number of plots where the species was recorded in each survey |
| shared_pool | Recorded in both surveys (used for community indices) |
| frequent_pool | Recorded in more than five plots in each survey (used for HMSC) |

## phylogeny.tre

Newick tree of all species, from V.PhyloMaker2 (GBOTB.extended.WP backbone, scenario S3).

## raw_surveys.zip

Original survey files (about 1975 and 2014–2019) read by the cleaning scripts in `pipeline/01_data_prep/01_cleaning/` of the code, with observer names, directions to the plots, photo file names and free-text notes removed. All other values are as recorded, including data-entry errors that the cleaning corrects. With these files, the cleaning can be reproduced from the original data: place the `_Org.csv` files in `pipeline/01_data_prep/01_cleaning/data/data_raw/01_raw_surveys/` and `MORA_PlotNames.csv` in `pipeline/01_data_prep/01_cleaning/data/data_raw/02_earlier_cleaning/`.

| File | Content |
|---|---|
| MORA_TreeData_Legacy_Org.csv, MORA_TreeData_Resurvey_Org.csv | Tree counts by 10-cm DBH class (past) and individual tree diameters (present) |
| MORA_VegData_Legacy_Org.csv, MORA_VegData_Resurvey_Org.csv | Understory cover per plot (past, one column per species code) and per 1-m² quadrat and 50-m² subplot (present) |
| MORA_Seedlings_Legacy_Org.csv, MORA_Seedlings_Resurvey_Org.csv | Tree seedling counts (not used in the analyses) |
| MORA_Plot_Legacy_Org.csv, MORA_Plot_Resurvey_Org.csv | Plot location, topography and survey information |
| MORA_PlantList_AllSurvey_Org.csv, MORA_PlantList_Legacy_Org.csv, MORA_PlantList_Resurvey_Org.csv | Species codes and names used in each survey |
| MORA_PlotNames.csv | Plot-name table of the earlier cleaning of these surveys, used to match plot identifiers between surveys |

## Data derived from other sources

- Species occurrences: GBIF.org, https://doi.org/10.15468/dl.cars7w and https://doi.org/10.15468/dl.ndtn3b.
- Climate: CHELSA v2.1 (Karger et al. 2017), monthly time series and 1981–2010 bioclimatic variables.
- Traits: TRY (Kattge et al. 2020; datasets in `try_references.csv`), GIFT (Weigelt et al. 2020) and BIEN (Maitner et al. 2018).
- Phylogeny: V.PhyloMaker2 (Jin & Qian 2022).

## Licence

CC0 1.0 (public domain dedication), as for all Dryad datasets.
