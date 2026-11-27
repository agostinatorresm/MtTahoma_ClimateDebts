
# 04_robustness

Robustness of the results to the difference in understory sampling between surveys. The past survey recorded occurrence and cover over the whole plot (250, 500 or 1000 m²); the present survey recorded cover in four 1-m² quadrats and occurrence in those quadrats plus a nested 50-m² subplot. These scripts quantify the consequences of that difference and test whether the main results hold when it is accounted for (Supporting Information S3).

## Scripts (run order)

| Script | What it does |
|---|---|
| `00_robustness_inputs.R` | Understory tables for the diagnostics, from the outputs of `01_cleaning` (same plots, species, present-day census and species complexes as the analyses) |
| `01a_accumulation_curves.R` | Sampling completeness of the present-day design: species and cover accumulation over 1–4 quadrats and the 50-m² subplot |
| `01b_patchy_species.R` | Species-level quadrat occupancy; flags patchy species |
| `02_bootstrapping.R` | Whether the change in cover between surveys exceeds the sampling variation of the present-day quadrats |
| `03_occurrence_sampling.R` | Sensitivity of occurrence to the difference in sampled area |
| `04_reliability_diagnostic.R` | One table per species combining the diagnostics, and the species lists for the HMSC robustness models |
| `05_rarefaction_occurrence.R` | Rarefaction of past occurrence to the present sampling effort (Wilson et al. 2024; Figure S3.2) |
| `06_rarefied_glm_occurrence.R` | Single-species GLMs of occurrence change on observed and rarefied data |
| `07_census_only_cover.R` | Cover change with abundance paired by plot |
| `08_diagnosis_figures.R` | Combined figures of the sampling diagnosis (Figures S3.1 and S3.3) |
| `09_rarefied_thermophilization.R` | Understory CTI and CMI with the past data reduced to the present sampling effort |
| `10_raw_cover_diagnosis.R` | Distribution of raw cover by type of record for the HMSC species |
| `11_glm_reduced_10_vs_100.R` | Check that 10 reduced past datasets give the same species changes as 100 |
| `12_rarefied_rate_models_bayes.R` | Bayesian mean rates of CTI and CMI with the reduced past data (Figure S4.1) |
| `13_table_S1_understory_design_bias.R` | Table S3.1: sources of bias, corrections, sensitivity analyses and cautions |

`00_run_all_robustness.R` runs scripts 00–10, with 10 before 08, whose figures use its output. Script 10 reads the HMSC inputs (`02_analysis/02_hmsc/code/01_prepare_hmsc_inputs.R`).

## Inputs and outputs

Inputs are read from `pipeline/01_data_prep/01_cleaning/outputs/`. Final tables are written to `outputs/`, intermediate files to `tmp_files/`, and the supplementary figures to `figures/`.
