# DESCRIPTION: Check the number of users in each AED type group, for different
#   age categories
# AUTHOR: Julia ROmanowska
# DATE CREATED: 2022-09-19
# DATE MODIFIED: 

# SETUP ----
library(here)
library(tidyr)
library(readr)
library(dplyr)
library(purrr)

# READ DATA ----
input_data <- read_delim(
  file = here("DATA", "data_ready_analysis_table2.txt"),
  delim = "\t"
) %>%
  mutate(
    across(
      c(minor_al_load, a1, a2, CC, IK_EPI, EPILEPSY2_c, PARITET_MFR:ALC_PREG, 
        SCL1_depr22, AED_POLY, ends_with("_MONO"), FOL_PRE_or_TRIM1,
        VIT_FOL_PRETRIM1TRIM2, starts_with("lang_imp"), starts_with("autistic"),
        AED_PREG_c, NOAED_PREG_c, CC_STATUS, SES_sum, BMIbefore_cat),
      ~ factor(.x)
    )
  )

# CHECK GROUPING ----
all_aed_groups <- c(
  stringr::str_subset(
    names(input_data), pattern = "_MONO"
  ),
  "AED_POLY"
)

## first - language impairment ----
walk(
  all_aed_groups,
  function(cur_aed){
    print(knitr::kable(
      input_data %>%
        filter(!is.na(lang_imp)) %>%
        janitor::tabyl(!!sym(cur_aed), lang_imp_age)
    ))
  }
)

## then, autistic traits ----
# I need to re-pivot the data first!
input_data <- input_data %>%
  pivot_wider(
    names_from = lang_imp_age,
    values_from = lang_imp,
    names_prefix = "lang_imp_"
  ) %>%
  pivot_longer(
    cols = starts_with("autistic_trait"),
    names_to = "autist_age",
    values_to = "autistic_trait"
  ) %>%
  mutate(autist_age = stringr::str_sub(autist_age, -3, -1))

walk(
  all_aed_groups,
  function(cur_aed){
    print(knitr::kable(
      input_data %>%
        filter(!is.na(autistic_trait)) %>%
        janitor::tabyl(!!sym(cur_aed), autist_age)
    ))
  }
)

