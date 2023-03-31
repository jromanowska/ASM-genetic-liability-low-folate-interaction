## SCRIPTS

### General

- `00_params.R`

common variables and parameters - check here for some definitions!

- `01_data_exploration.R`

exploring the entire dataset; filtering the individuals and merging
with PRS dataset; exploring the case-control status; plotting some
histograms;
creating table 1 (descriptive statistics)

### Part 1: Risk of language impairment

- `02_prs_folate_risk_lang_impair_check_data.R`

checking dataset and preparing for analyses

- `03_function_tidy_res.R`

functions defined here:

    - `tidy_res_robust`
    - `compare_all_adjustments`
    - `extract_relevant_adjustments`

These will be used during the preparation of tables, in the subsequent .Rmd files

- `04_1_models_run.R` - functions to run the main models

- `04_prs_folate_risk_lang_impair_models_run.R`    
main analyses:

    - running the adjusted models for each age;
    - stratifying by use of AED;
    - no stratification, using entire data, but including interaction;
    - using `04_1_models_run.R`

- `05_prs_folate_risk_lang_impair_report.Rmd` - reporting final results

- `06_check_women_group_vs_child_age.R` - checking if the language impairment
reporting frequency depend on the child's age (i.e., if there is a significant
difference between the questionnaires)

### Part 2: Risk of autism

- `07_data_exploration_tabl3.R` - exploration of data for analyses with autism
as an outcome; similar to script `02_prs_folate_risk_lang_impair_check_data.R`
  
- `10_prs_folate_risk_autist_models_run.R`    
main analyses:

    - running the adjusted models for each age;
    - stratifying by use of AED;
    - no stratification, using entire data, but including interaction

- `11_prs_folate_risk_autism_report.Rmd` - reporting final results


### Part 3: OR of autism/language impairment in various subgroups

- `12_check_AED_stratif.R`

### Part 4: predicted outcome vs (continuous) PRS or MTHFR genotype strata,
  stratified by AED use

- `13_interaction_plot_AED_use-PRS.R` - 
  preparing data for plotting: running the crude models stratified by AED use
  and extracting prediction values for the entire span of each PRS (three types:
  PRS automatically chosen by PRSice2, PRS constructed manually, and MTHFR)
  
- `14_interaction_plots_report.Rmd` - report with final plots
  
### Part 5: Examine the relationship between the PRS score and folate concentration
(“fol_hydro”) stratified for types of ASM treatment with non-parametric
correlation analysis

- `15_corr_PRS-fol_hydro_AED_strat.Rmd`
