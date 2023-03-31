# DESCRIPTION: Plot the interaction plots (interaction between PRS and AED use)
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-09-19
# DATE MODIFIED: 2022-11-25

# SETUP ----
library(here)
library(tidyr)
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(naniar)
library(emmeans)

# READ DATA ----
input_data <- read_delim(
  here("DATA", "data_ready_analysis_table2.txt"),
  delim = "\t"
) %>%
  mutate(
    across(
      c(minor_al_load, a1, a2, CC, IK_EPI, EPILEPSY2_c, FOL_PRE_or_TRIM1,
        starts_with("lang_imp"), starts_with("autistic"),
        AED_PREG_c, NOAED_PREG_c, CC_STATUS, MTHFR_2strat),
      ~ factor(.x)
    )
  )
source("00_params.R")

# PLOTTING ----

# this function will return a list of lists of interaction plots, one for each
#   prs type and for each age; each element of the inner list contains:
#   - glm_run - glm object with result of run of a model: outcome ~ prs_type * any_AED
#   - trends - result of 'emtrends' function; showing whether the interaction
#              is plausible
#   - plot - ggplot2 with interaction plot
plot_interact_prs <- function(
    outcome, # here: 'lang_imp' or 'autistic_trait'
    data, # input data
    all_ages, # all age categories (character vector)
    all_prs_types, # all prs types (character vector with variable names)
    outcome_age_cat # name of the variable in 'data' that has age categories
){
  all_results <- map(all_ages, function(cur_age){
    cat("--- CURRENT AGE CAT.: ", cur_age, " ---\n")

    plot_data <- data %>%
      filter(!!sym(outcome_age_cat) == cur_age)
    
    # check how the missing patterns look like in this age category
    print(
      knitr::kable(
        plot_data %>%
          janitor::tabyl(any_AED, !!sym(outcome))
      )
    )
    
    cur_age_res <- map(all_prs_types, function(cur_prs_type){
      # cur_prs_type <- all_prs_types[1]
      # cur_age <- all_ages[2]
      cat("--- CURRENT PRS TYPE: ", cur_prs_type, " ---\n")

      if(isa(data[[cur_prs_type]], "factor")){
        print(
          knitr::kable(
            plot_data %>%
              janitor::tabyl(any_AED, !!sym(cur_prs_type))
          )
        )
      } else {
        print(
          ggplot(plot_data, aes(!!sym(cur_prs_type), !!sym(outcome))) +
            geom_miss_point() +
            labs(title = paste0("Age: ", cur_age, ", outcome: ", outcome))
        )
      }
      
      plot_data <- plot_data %>%
        select(all_of(c(cur_prs_type, outcome)), any_AED)
            
      # first, run simple model
      cur_formula <- paste0(outcome, " ~ ", cur_prs_type, "*any_AED")
      aed_prs_interact <- glm(
        formula = cur_formula,
        data = plot_data,
        family = "binomial"
      )
      # this is needed to use 'easystats' functions on the result afterwards
      aed_prs_interact$call$formula <- as.formula(cur_formula)
      
      # check the interaction - this differs for categorical or continuous PRS
      if(isa(data[[cur_prs_type]], "factor")){
        aed_prs_trends <- emmeans(
          aed_prs_interact,
          specs = as.formula(paste0(" ~ any_AED*", cur_prs_type))
        )
        
        int_plot_data <- emmip(
          aed_prs_interact,
          as.formula(paste0("any_AED ~", cur_prs_type)),
          CIs = TRUE, # this will calculate the CI as range (LCL, UCL)
          plotit = FALSE # return data only
        )
        int_plot <- ggplot(int_plot_data) +
          aes(!!sym(cur_prs_type), yvar) +
          geom_line(aes(col = any_AED, group = any_AED)) +
          geom_pointrange(aes(ymax = UCL, ymin = LCL, col = any_AED))
      } else {
        # continuous PRS
        aed_prs_trends <- emtrends(
          aed_prs_interact,
          specs = pairwise ~ any_AED,
          var = cur_prs_type
        )
        
        # plot the interaction
        coords_list <- list(
          cur_prs = seq(
            min(plot_data[[cur_prs_type]]),
            max(plot_data[[cur_prs_type]]),
            length.out = 10
          ),
          any_AED = c(FALSE, TRUE)
        )
        names(coords_list) <- c(cur_prs_type, "any_AED")
      
        int_plot_data <- emmip(
          aed_prs_interact,
          as.formula(paste0("any_AED ~", cur_prs_type)),
          at = coords_list,
          CIs = TRUE, # this will calculate the CI as range (LCL, UCL)
          plotit = FALSE # return data only
        )
        int_plot <- ggplot(int_plot_data) +
          aes(!!sym(cur_prs_type), yvar) +
          geom_line(aes(col = any_AED)) +
          geom_ribbon(aes(ymax = UCL, ymin = LCL, fill = any_AED), alpha = 0.3)
      }
      
      return(
        list(
          glm_run = aed_prs_interact,
          trends = aed_prs_trends,
          plot = int_plot
        )
      )
    })
    names(cur_age_res) <- all_prs_types
    return(cur_age_res)
  })
  names(all_results) <- all_ages
  return(all_results)
}

## first, language impairment ----
outcome <- "lang_imp"
all_ages <- levels(input_data$lang_imp_age)

interact_lang_imp_list <- plot_interact_prs(
  outcome,
  input_data,
  all_ages,
  all_prs_types,
  "lang_imp_age"
)

## next, autistic trait ----
outcome <- "autistic_trait"

# re-pivot data
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

all_ages <- unique(input_data$autist_age)

interact_autist_list <- plot_interact_prs(
  outcome,
  input_data,
  all_ages,
  all_prs_types,
  "autist_age"
)

# SAVE RESULTS ----
saveRDS(
  interact_lang_imp_list,
  here("RESULTS", "interact_lang_imp_list.rds")
)

saveRDS(
  interact_autist_list,
  here("RESULTS", "interact_autist_list.rds")
)
