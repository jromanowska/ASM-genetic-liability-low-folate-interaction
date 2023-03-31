# DESCRIPTION: Running the models for table 3
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-08-25
# DATE MODIFIED: 2022-12-22

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
all_res_types # read from '00_params.R'
all_prs_types

input_data <- read_delim(
  file = here("DATA", "data_ready_analysis_table3.txt"),
  delim = "\t"
  ) %>%
  select(-starts_with("lang_imp")) %>%
  select(iid:NOAED_PREG_c, PLURAL, DUPLIKAT, FOL_PRE_or_TRIM1,
         FOLAT, s_folat, MTHFR_2strat:autistic_trait) %>%
  mutate(
    across(
      c(minor_al_load, a1, a2, CC, IK_EPI, EPILEPSY2_c, FOL_PRE_or_TRIM1,
        autistic_trait, AED_PREG_c, NOAED_PREG_c, CC_STATUS, any_AED, MTHFR_2strat),
      ~ factor(.x)
    )
  )
all_ages <- unique(input_data$autist_age)

# RUN MODELS WITH STRATIFICATION ----
# the analyses 'mthfr' and 'fil_int' have additional stratification
strat_mthfr_2strat <- levels(input_data$MTHFR_2strat)
strat_fol_int <- levels(input_data$FOL_PRE_or_TRIM1)
add_strata <- c(list(NA), list(NA), list(NA),
                list(strat_mthfr_2strat), list(strat_mthfr_2strat),
                list(strat_fol_int), list(strat_fol_int))
names(add_strata) <- all_res_types
add_strata


tab3_all_res <- run_main_model_strat( #function defined in "04_1_models_run.R"
  data = input_data,
  # choosing only the models without interaction:
  res_types = stringr::str_subset(all_res_types, "_interact", negate = TRUE),
  explanatory_vars = explanatory_vars,
  all_prs_types = all_prs_types,
  all_ages = all_ages,
  add_strata = add_strata,
  dependent_var = "autistic_trait",
  dependent_age_var = "autist_age"
)
  
saveRDS(
  object = tab3_all_res,
  file = here("RESULTS", "analysis_2_fit_adj_model_list.rds")
)

# RUN MODELS WITH INTERACTION ----
# I will take entire dataset, but now case/control status will be only
# use and no use of AED, PRS will interact with the var
# that indicates taking AED (independent from CC_STATUS)
explanatory_vars <- c("any_AED", explanatory_vars)
explanatory_names <- c("any_AED", explanatory_names)
names(explanatory_vars) <- explanatory_names

tab3_interact_all_res <- run_main_model( #function from 04_1_models_run.R
  data = input_data,
  res_type = "res_prelim_interact",
  explanatory_vars = explanatory_vars,
  all_prs_types = all_prs_types,
  all_ages = all_ages,
  dependent_var = "autistic_trait",
  dependent_age_var = "autist_age"
)

# trying also without 'EPILEPSY2_c' variable
saveRDS(
  object = tab3_interact_all_res,
  file = here("RESULTS", "analysis_2_fit_adj_model_interact_list.rds")
)
