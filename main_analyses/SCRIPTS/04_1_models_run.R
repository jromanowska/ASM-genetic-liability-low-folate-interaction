#' Run the main models
#' 
#' Take the data and covariates that were chosen to change the crude OR and
#'   construct models that use stratified data to run models on each age
#'  separately. All arguments are necessary.
#' 
#' @param data - dataset to run the models
#' @param res_types - character vector with names of the result groups
#' @param explanatory_vars - named character vector with all the explanatory
#'   variables - this will be used to construct the formula, names will be used
#'   when extracting variables for each model
#' @param all_prs_types - character vector with names of the PRS variables
#' @param all_ages - character vector with levels of age categories
#' @param add_strata - list of lists of additional strata, useful for the 
#'   analyses that are stratified by folate intake and MTHFR_SNP genotype
#' @param dependent_var - name of the dependent variable
#' @param dependent_age_var - name of the variable containing age categories
#' 
#' @return list of lists of tibbles with results (model objects), per age
#'   category and per result type (from res_types)

run_main_model_strat <- function(
  data,
  res_types,
  explanatory_vars,
  all_prs_types,
  all_ages,
  add_strata,
  dependent_var,
  dependent_age_var
){
  crude_form <- as.formula(paste(dependent_var, "~ CC_STATUS"))
  explanatory_names <- names(explanatory_vars)
  out <- map(
    res_types,
    function(res_type){
      cat("\n-- RESULTS TYPE: ", res_type, "--\n")
      # I need to make it explicit here because of filtering of data
      if(stringr::str_detect(string = res_type, pattern = "_no_aed")){
        cat("--- CURRENT STRATUM: women who did not take AED ---\n")
        cur_data <- data %>%
          filter(CC_STATUS != "epi_AED") %>%
          mutate(CC_STATUS = forcats::fct_drop(CC_STATUS)) %>%
          mutate(CC_STATUS = forcats::fct_relevel(CC_STATUS, "no_epi", "epi_no_AED"))
      } else {
        cat("--- CURRENT STRATUM: women who took AED ---\n")
        cur_data <- data %>%
          filter(CC_STATUS != "epi_no_AED") %>%
          mutate(CC_STATUS = forcats::fct_drop(CC_STATUS)) %>%
          mutate(CC_STATUS = forcats::fct_relevel(CC_STATUS, "no_epi", "epi_AED"))
      }
      cur_data <- cur_data %>%
        select(iid, a1, a2, CC_STATUS, m_id,
               all_of(c(explanatory_names, all_prs_types, dependent_var, dependent_age_var)))
      
      cur_res_age_all <- map(
        all_ages,
        function(cur_age){
          cat("--- CURRENT AGE: ", cur_age, " ---\n")

          if(stringr::str_detect(res_type, "mthfr|fol_int")){
            # here, there are additional strata level
            all_strata <- add_strata[[res_type]]
            # I also need the name of the variable used for stratification
            #  hard-coded for now:
            stratum_var <- if_else(
              stringr::str_detect(res_type, "mthfr"),
              "MTHFR_2strat",
              "FOL_PRE_or_TRIM1"
            )
            
            cur_res_age <- map(
              all_strata,
              function(cur_stratum){
                cat("------ current stratum: ", cur_stratum, "------\n")
                cur_data_age <- cur_data %>%
                  filter(!!sym(dependent_age_var) == cur_age &
                           !!sym(stratum_var) == cur_stratum) %>%
                  select(iid, CC_STATUS, m_id, all_of(
                    c(explanatory_vars, dependent_var, all_prs_types))
                  )
                
                # since the MTHFR-stratification includes genotype data,
                # we won't be looking at PRS types...
                if(stringr::str_detect(res_type, "mthfr")){
                  cur_formula <- as.formula(
                    paste(dependent_var, "~ CC_STATUS +",
                          paste(as.character(explanatory_vars),
                                collapse = "+")
                    )
                  )
                  cur_res_strat <- glm(
                    formula = cur_formula,
                    data = cur_data_age,
                    family = "binomial"
                  )
                  cur_res_strat$call$formula <- cur_formula
                  print(cur_res_strat)
                  
                  # in the other stratified analysis we get a list, so
                  #   here, we need to also return a list
                  cur_res_strat <- list(cur_res_strat)
                  names(cur_res_strat) <- cur_stratum
                  
                } else { # however, for stratification by folate intake,
                  # any PRS can be considered in the formula
                  cur_res_strat <- run_per_prs_type(
                    all_prs_types,
                    cur_data_age,
                    dependent_var,
                    stringr::str_subset(explanatory_vars, stratum_var, negate = TRUE)
                  )
                }
                
                # add crude model
                cur_formula <- as.formula(
                  paste(dependent_var, "~ CC_STATUS")
                )
                crude_res <- glm(
                  formula = cur_formula,
                  data = cur_data_age,
                  family = "binomial"
                )
                crude_res$call$formula <- cur_formula
                print(crude_res)
                cur_res_strat <- c(cur_res_strat, list(crude = crude_res))
                return(cur_res_strat)
              }
            )
            names(cur_res_age) <- all_strata
            
          } else { # running for all other models, without additional stratification
            cur_data_age <- cur_data %>%
              filter(!!sym(dependent_age_var) == cur_age) %>%
              select(iid, CC_STATUS, m_id, all_of(c(
                explanatory_vars, dependent_var, all_prs_types
                )
              ))
          
            cur_res_age <- run_per_prs_type(
              all_prs_types,
              cur_data_age,
              dependent_var,
              explanatory_vars
            )
          
            # add crude model
            cur_formula <- as.formula(
              paste(dependent_var, "~ CC_STATUS")
            )
            crude_res <- glm(
              formula = cur_formula,
              data = cur_data_age,
              family = "binomial"
            )
            crude_res$call$formula <- cur_formula
            print(crude_res)
            cur_res_age <- c(cur_res_age, list(crude = crude_res))
          }
          
          return(cur_res_age)
        }
      )
      names(cur_res_age_all) <- all_ages
      return(cur_res_age_all)
    }
  )
  names(out) <- res_types
  return(out)
}

#' Small function running the model per each PRS type
#' 
run_per_prs_type <- function(
  all_prs_types,
  cur_data_age,
  dependent_var,
  explanatory_vars
){
  # since these are only the covariates that were changing the crude
  #   estimate, there might not be any PRS...
  cur_res_age <- map(all_prs_types, function(prs_type){
    cur_formula <- as.formula(
      paste(dependent_var, "~ CC_STATUS +",
          paste(as.character(explanatory_vars), collapse = "+"),
          " + ", prs_type)
    )
    cur_fit <- glm(
      formula = cur_formula,
      data = cur_data_age,
      family = "binomial"
    )
    cur_fit$call$formula <- cur_formula
    print(cur_fit)
    return(cur_fit)
  })
  names(cur_res_age) <- all_prs_types
  
  # running once more, without PRS
  cur_formula <- as.formula(
    paste(dependent_var, "~ CC_STATUS +",
        paste(as.character(explanatory_vars), collapse = "+")
   )
  )
  cur_res_no_prs <- glm(
    formula = cur_formula,
    data = cur_data_age,
    family = "binomial"
  )
  cur_res_no_prs$call$formula <- cur_formula
  print(cur_res_no_prs)
  cur_res_age <- c(
    cur_res_age,
    list(no_prs = cur_res_no_prs)
  )
  
  return(cur_res_age)
}

#' Run the main models - without stratification
#' 
#' Take the data and covariates that were chosen to change the crude OR and
#'   construct models that use stratified data to run models on each age
#'  separately. All arguments are necessary.
#' 
#' @param data - dataset to run the models
#' @param res_type - name of the result group to extract covariates from 
#'   relev_covs_prelim
#' @param explanatory_vars - named character vector with all the explanatory
#'   variables - this will be used to construct the formula, names will be used
#'   when extracting variables for each model
#' @param all_prs_types - character vector with names of the PRS variables
#' @param all_ages - character vector with levels of age categories
#' @param dependent_var - name of the dependent variable
#' @param dependent_age_var - name of the variable containing age categories
#' 
#' @return list of tibbles with results (model objects), one per age category

run_main_model <- function(
  data,
  res_type,
  explanatory_vars,
  all_prs_types,
  all_ages,
  dependent_var,
  dependent_age_var
){
  # crude_form <- as.formula(paste(dependent_var, "~ EPILEPSY2_c"))
  crude_form <- as.formula(paste(dependent_var, "~ any_AED"))
  explanatory_names <- names(explanatory_vars)
  out <- map(
    all_ages,
    function(cur_age){
      cat("--- CURRENT AGE: ", cur_age, " ---\n")
      cur_res_type <- res_type
      cur_data_age <- data %>%
        filter(!!sym(dependent_age_var) == cur_age) %>%
        select(iid, EPILEPSY2_c, any_AED, m_id, all_of(explanatory_vars),
               all_of(c(all_prs_types, dependent_var)))
      
      # since these are only the covariates that were changing the crude
      #   estimate, there might not be any PRS...
      cur_res_age <- map(all_prs_types, function(prs_type){
        cat("PRS TYPE: ", prs_type, "\n")
        # check the model with the PRS
        cur_formula <- as.formula(
          # paste(dependent_var, " ~ EPILEPSY2_c +",
          paste(dependent_var, " ~",
                paste(as.character(explanatory_vars), collapse = "+"),
                " + ", prs_type, " + ", prs_type, "*any_AED")
        )
        cur_fit <- glm(
          formula = cur_formula,
          data = cur_data_age,
          family = "binomial"
        )
        cur_fit$call$formula <- cur_formula
        print(cur_fit)
        return(cur_fit)
      })
      names(cur_res_age) <- all_prs_types
      
      # running once more, without PRS
      cur_formula <- as.formula(
        # paste(dependent_var, " ~ EPILEPSY2_c +",
        paste(dependent_var, " ~ ",
              paste(as.character(explanatory_vars), collapse = "+")
        )
      )
      cur_res_no_prs <- glm(
        formula = cur_formula,
        data = cur_data_age,
        family = "binomial"
      )
      cur_res_no_prs$call$formula <- cur_formula
      cur_res_age <- c(
        cur_res_age,
        list(no_prs = cur_res_no_prs)
      )
      
      return(cur_res_age)
    }
  )
  names(out) <- all_ages
  return(out)
}