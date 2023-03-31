# DESCRIPTION: Check the data from Elisabeth; filter the women that have PRS;
#   check specifically for analysis in table 3
# AUTHOR: Julia ROmanowska
# DATE CREATED: 2022-08-16
# DATE MODIFIED: 2022-08-25

# SETUP ----
library(here)
library(tidyr)
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

# READ DATA ----
input_data <- read_delim(
  file = here("DATA", "data_ready_analysis_table2.txt"),
  delim = "\t"
) %>%
  mutate(
    across(
      c(minor_al_load, a1, a2, CC, IK_EPI, EPILEPSY2_c, FOL_PRE_or_TRIM1,
        starts_with("lang_imp"), starts_with("autistic"), MTHFR_2strat,
        AED_PREG_c, NOAED_PREG_c, CC_STATUS),
      ~ factor(.x)
    )
  )
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
  mutate(autist_age = stringr::str_sub(autist_age, -3, -1)) %>%
  select(-starts_with("lang_imp_"))
all_ages <- unique(input_data$autist_age)

skimr::skim(input_data)

# Check missingness ----
purrr::walk(
  all_ages,
  function(cur_age){
    cur_data <- input_data %>%
      filter(autist_age == cur_age)
    
    cat("\n ---- AGE: ", cur_age, " ----\n")
    summ_miss_case <- naniar::miss_case_summary(cur_data)
    print(knitr::kable(
      summ_miss_case %>%
        summarise(
          max_pct_miss = max(pct_miss),
          min_pct_miss = min(pct_miss),
          mean_pct_miss = mean(pct_miss),
          median_pct_miss = median(pct_miss),
          n_any_miss = nrow(filter(., pct_miss > 0))
        ),
      caption = "Summary of missing values per row"
    ))
    print(summ_miss_case)
  }
)

purrr::walk(
  all_ages,
  function(cur_age){
    cur_data <- input_data %>%
      filter(autist_age == cur_age)
    
    cat("\n ---- AGE: ", cur_age, " ----\n")
    summ_miss_var <- naniar::miss_var_summary(cur_data)
    print(knitr::kable(
      summ_miss_var %>%
        summarise(
          max_pct_miss = max(pct_miss),
          min_pct_miss = min(pct_miss),
          mean_pct_miss = mean(pct_miss),
          median_pct_miss = median(pct_miss),
          n_any_miss = nrow(filter(., pct_miss > 0))
        ),
      caption = "Summary of missing values per variable"
    ))
    print(summ_miss_var)
  }
)

# stratification variables ----
walk(c("FOL_PRE_or_TRIM1", "MTHFR_2strat", "any_AED"), function(strat_var){
  cat("\n", strat_var, " vs autistic_trait \n")
  walk(all_ages, function(cur_age){
    cat("\t cur_age: ", cur_age, "\n")
    print(knitr::kable(
      janitor::tabyl(input_data, !!sym(strat_var), autistic_trait)
    ))
  })
})

# Check the explanatory variables ----
purrr::walk(
  all_ages,
  function(cur_age){
    cur_data <- input_data %>%
      filter(autist_age == cur_age)
    
    cat("\n ---- AGE: ", cur_age, " ----\n")
    print(
      janitor::tabyl(cur_data, autistic_trait)
    )
    cat("\n")
    print(
      janitor::tabyl(cur_data, CC_STATUS, autistic_trait)
    )
    cat("\n")
  }
)


## Save data for later ----
write_delim(
  input_data,
  file = here("DATA", "data_ready_analysis_table3.txt"),
  delim = "\t"
)
