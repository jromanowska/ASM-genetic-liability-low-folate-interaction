# DESCRIPTION: Running the models for table 2
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-06-30
# DATE MODIFIED: 2022-12-07

# SETUP ----
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(here)

# source parameters and constants:
source("00_params.R")
# source functions:
#   - run_main_model_strat
#   - run_main_model
source("04_1_models_run.R")

# READ DATA ----
# result types
all_res_types # read from '00_params.R'

input_data <- read_delim(
  here("DATA", "data_ready_analysis_table2.txt"),
  delim = "\t"
  ) %>%
  mutate(
    across(
      c(minor_al_load, a1, a2, CC, IK_EPI, EPILEPSY2_c, FOL_PRE_or_TRIM1,
        starts_with("lang_imp"), starts_with("autistic"), any_AED,
        AED_PREG_c, NOAED_PREG_c, CC_STATUS, MTHFR_2strat),
      ~ factor(.x)
    )
  )
all_ages <- levels(input_data$lang_imp_age)
all_ages

# RUN MODELS WITH STRATIFICATION ----

# the analyses 'mthfr' and 'fil_int' have additional stratification
strat_mthfr_2strat <- levels(input_data$MTHFR_2strat)
strat_fol_int <- levels(input_data$FOL_PRE_or_TRIM1)
add_strata <- c(list(NA), list(NA), list(NA),
                list(strat_mthfr_2strat), list(strat_mthfr_2strat),
                list(strat_fol_int), list(strat_fol_int))
names(add_strata) <- all_res_types
add_strata

# run all models
tab2_all_res <- run_main_model_strat( #function defined in "04_1_models_run.R"
  data = input_data,
  # choosing only the models without interaction:
  res_types = stringr::str_subset(all_res_types, "_interact", negate = TRUE),
  explanatory_vars = explanatory_vars,
  all_prs_types = all_prs_types,
  all_ages = all_ages,
  add_strata = add_strata,
  dependent_var = "lang_imp",
  dependent_age_var = "lang_imp_age"
)

saveRDS(
  object = tab2_all_res,
  file = here("RESULTS", "analysis_1_fit_adj_model_list.rds")
)

# RUN MODELS WITH INTERACTION ----
# I will take entire dataset, but now case/control status will be only
# epilepsy/no epilepsy (EPILEPSY2_c); then, PRS will interact with the var
# that indicates taking AED (independent from CC_STATUS)
explanatory_vars <- c("any_AED", explanatory_vars)
explanatory_names <- c("any_AED", explanatory_names)
names(explanatory_vars) <- explanatory_names

tab2_interact_all_res <- run_main_model( #function from 04_1_models_run.R
  data = input_data,
  res_type = "res_prelim_interact",
  explanatory_vars = explanatory_vars,
  all_prs_types = all_prs_types,
  all_ages = all_ages,
  dependent_var = "lang_imp",
  dependent_age_var = "lang_imp_age"
)

# trying also without 'EPILEPSY2_c' variable
saveRDS(
  object = tab2_interact_all_res,
  file = here("RESULTS", "analysis_1_fit_adj_model_interact_list.rds")
)

# # trying with additional stratification by use of folate
# explanatory_vars <- stringr::str_subset(
#   explanatory_vars,
#   "FOL_PRE_or_TRIM1",
#   negate = TRUE
# )
# tab2_interact_fol_users_res <- run_main_model( #function from 04_1_models_run.R
#   data = input_data %>%
#     filter(FOL_PRE_or_TRIM1 == 1),
#   res_type = "res_prelim_interact",
#   explanatory_vars = explanatory_vars,
#   all_prs_types = all_prs_types,
#   all_ages = all_ages,
#   dependent_var = "lang_imp",
#   dependent_age_var = "lang_imp_age"
# )
# 
# tab2_interact_no_fol_res <- run_main_model( #function from 04_1_models_run.R
#   data = input_data %>%
#     filter(FOL_PRE_or_TRIM1 == 0),
#   res_type = "res_prelim_interact",
#   explanatory_vars = explanatory_vars,
#   all_prs_types = all_prs_types,
#   all_ages = all_ages,
#   dependent_var = "lang_imp",
#   dependent_age_var = "lang_imp_age"
# )
