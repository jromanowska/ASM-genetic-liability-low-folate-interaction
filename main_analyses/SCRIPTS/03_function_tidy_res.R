# DESCRIPTION: Functions to perform tidying of glm results, calculating robust
#    CI; gathering all in tibbles; checking difference with the crude model.
#AUTHOR: Julia Romanowska
#DATE CREATED: 2022-09-23
#DATE MODIFIED: 2022-12-10


#' Tidy the results and add robust CI
#' 
#' @param cur_model - result of glm run
#' @param cluster_var - variable indicating clusters within the data (numeric vec.)
#' @param dependent - name of the dependent variable
#' @param exp - should the results be given as coef (exp = FALSE) or OR
#'    (exp = TRUE, default)?
#' @param interaction - vector with two elements: names of the variables
#'    involved in the interaction term
#' 
#' @return tibble with tidied results
#' 
tidy_res_robust <- function(
  cur_model,
  cluster_var,
  dependent,
  exp = TRUE,
  interaction = NULL
){
  robust_ci <- lmtest::coefci( # calculate robust CI
    cur_model,
    vcov. = sandwich::vcovCL(cur_model, cluster = cluster_var)
  )
  if(exp){
    robust_ci <- exp(robust_ci)
  }
  robust_ci_tbl <- as_tibble(robust_ci, rownames = "term")
  colnames(robust_ci_tbl) <- c("term", "conf.low", "conf.high")
  
  tidy_res <- broom::tidy(cur_model, exp = exp) %>%
    select(term, estimate) %>%
    left_join(robust_ci_tbl, by = "term")
  
  # check whether the estimates include interaction term
  none_numeric <- FALSE
  interact_term <- ""
  if(!is.null(interaction)){
    # just in case the given vector is larger:
    interaction <- interaction[1:2]
    
    interact_term <- paste0(interaction, collapse = ":")

    none_numeric <- cur_model$model %>%
      select(!!sym(interaction[1]), !!sym(interaction[2])) %>%
      none(is.numeric)
      
    if(none_numeric){
      cur_model$model <- cur_model$model %>%
        mutate(across(c(!!sym(interaction[1]), !!sym(interaction[2])),
                      as.factor)) %>%
        finalfit::ff_interaction(
          !!sym(interaction[1]), !!sym(interaction[2]),
          var_sep = ":"
        )
      # for the variable to be recognized as interaction term by {finalfit},
      #   the two parts should be separated by "_"
      names(cur_model$model)[length(names(cur_model$model))] <- 
        str_replace_all(pattern = "_",
             replacement = ".",
             string = interact_term) %>%
        str_replace(pattern = ":", replacement = "_")
    }
  }
  explanatory <- setdiff(
    names(cur_model$model),
    dependent
  )
  summary_cc <- suppressMessages(
    finalfit::summary_factorlist(
      cur_model$model,
      explanatory = explanatory,
      dependent = dependent,
      fit_id = TRUE,
      column = FALSE # percentages will be computed by row, not by column
    )
  )
  
  if(!is.null(interaction)){
    if(none_numeric){
      # need to change the name of the interaction term in tidy table
      interaction_OR <- tidy_res %>%
        tail(1)
      tidy_res <- tidy_res %>%
        head(-1) %>%
        add_row(
          term = summary_cc %>% tail(1) %>% pull(fit_id),
          estimate = interaction_OR$estimate,
          conf.low = interaction_OR$conf.low,
          conf.high = interaction_OR$conf.high
        )

    } else {
      #there was a problem when computing mean and SD for prs_auto!
      cur_prs <- str_subset(interaction, "prs")
      mean_prs0 <- cur_model$model %>%
        filter(!!sym(dependent) == 0) %>%
        pull(!!sym(cur_prs)) %>%
        mean()
      sd_prs0 <- cur_model$model %>%
        filter(!!sym(dependent) == 0) %>%
        pull(!!sym(cur_prs)) %>%
        sd()
      mean_prs1 <- cur_model$model %>%
        filter(!!sym(dependent) == 1) %>%
        pull(!!sym(cur_prs)) %>%
        mean()
      sd_prs1 <- cur_model$model %>%
        filter(!!sym(dependent) == 1) %>%
        pull(!!sym(cur_prs)) %>%
        sd()
      
      summary_cc <- summary_cc %>%
        filter(label != cur_prs) %>%
        add_row(
          label = cur_prs,
          levels = "Mean (SD)",
          `0` = sprintf("%.2e (%.2e)", mean_prs0, sd_prs0),
          `1` = sprintf("%.2e (%.2e)", mean_prs1, sd_prs1),
          fit_id = cur_prs
        ) %>%
        # need to add name of the interaction to summary_cc
        add_row(
          label = interact_term,
          fit_id = tidy_res %>% tail(1) %>% pull(term)
        )
      
      # tidy_res <- as_tibble(summary_cc) %>%
      #   full_join(
      #     tidy_res,
      #     by = c("fit_id" = "term")
      #   ) %>%
      #   rename(term = fit_id) %>%
      #   filter(term != "(Intercept)")
    }
  } #else {
    tidy_res <- as_tibble(summary_cc) %>%
      left_join(
        tidy_res,
        by = c("fit_id" = "term")
      ) %>%
      rename(term = fit_id)
  # }
  
  return(tidy_res)
}


#' Create a tibble that compares coefficients for the selected terms between
#' all the models
#'
#' @param tidy_list - list of tidied results for each model (named)
#' @param filter_terms - character vector with all the terms to filter from 
#'   each element of the list
#'   
#' @return tibble gathering all the coefficients, with variables:
#'   - adj - name of the model
#'   - estimate - coefficient
#'   - conf.low, conf.high - CI

compare_all_adjustments <- function(tidy_list, filter_terms){
  all_adj_compared <- map(names(tidy_list), function(cur_adj){
    cur_model <- tidy_list[[cur_adj]]
    out <- cur_model %>%
      filter(term %in% filter_terms) %>%
      select(estimate, conf.low, conf.high) %>%
      tibble::add_column(adj = cur_adj, .before = "estimate")
    return(out)
  }) %>%
    bind_rows()
  return(all_adj_compared)
}

#' Compare the adjusted OR with the crude one and extract the adjustments that
#'   change the crude OR by more than a given threshold
#'   
#' @param estimates - tibble with crude and adjusted estimates (obtained from
#'   'compare_all_adjustments()')
#' @param threshold - how large should the change be to accept this adjustment?
#'   (default: 0.01)
#' @param filter_crude - string: name of the crude estimate in the 'estimates'
#'   (default: "crude")
#' 
#' @return tibble with only relevant adjustments and estimates

extract_relevant_adjustments <- function(
    estimates,
    threshold = 0.01,
    filter_crude = "crude"
){
  cur_crude_est <- estimates %>%
    filter(adj == filter_crude) %>%
    pull(estimate)
  out <- estimates %>%
    mutate(est_diff_crude = abs(estimate - !!cur_crude_est)) %>%
    filter(est_diff_crude > threshold)
  out <- bind_rows(
    estimates %>%
      filter(adj == filter_crude),
    out
  )
  return(out)
}


tidy_and_extract <- function(
    cur_data,
    cur_res_type,
    cur_age,
    main_var = "CC_STATUS",
    dept_var = "lang_imp"
){
  map(
    names(cur_data),
    function(cur_adj){
      # cat("\t\t cur_adj: ", cur_adj, "\n")
      cur_fit <- cur_data[[cur_adj]]
      
      if(is.null(cur_fit)){
        return(NULL)
      }
      cur_formula <- as.character(cur_fit$formula)
      cur_formula <- paste(
        cur_formula[2], cur_formula[1], cur_formula[3]
      )
      cur_adj_table <- tidy_res_robust(
        cur_fit,
        cur_fit$data$m_id,
        dept_var
      )
      
      cur_adj_table <- cur_adj_table %>%
        filter(stringr::str_starts(
          term,
          as.character(main_var))
        ) %>%
        tibble::add_column(
          age = cur_age,
          res_type = cur_res_type,
          formula = cur_formula,
          adj = cur_adj
        ) %>%
        select(-index)
      return(cur_adj_table)
    }
  ) %>%
    bind_rows()
}

#' Creating nice table with results, a.k.a. table 2
#' 
#' @param data - current estimates
#' @param subtitle - title for this table
#' @param levels - levels of the dependent variable
#' 
#' @return flextable that looks nice
#' 
make_table2 <- function(
    data,
    subtitle,
    levels = list(
      `0`= "no lang.impair.",
      `1` = "lang.impair."
    )
){
  small_border <- officer::fp_border(color = "gray", width = 1)
  header_border <- officer::fp_border(color = "gray30", width = 2)
  
  as_grouped_data(
    data %>% 
      select(-label) %>%
      select(adj, everything()) %>%
      mutate(
        estimate = ifelse(
          is.na(estimate),
          yes = "",
          no = sprintf(
            "%.2f (%.2f - %.2f)", estimate, conf.low, conf.high
          )
        ),
        adj = ifelse(
          is.na(adj),
          yes = "crude",
          no = adj
        )
      ) %>%
      arrange(age, adj) %>%
      select(-conf.low, -conf.high),
    groups = "age"
  ) %>%
    as_flextable() %>% 
    bold(bold = TRUE, part = "header") %>%
    align(i = ~ !is.na(age), align = "center") %>% 
    bold(i = ~ !is.na(age)) %>%
    set_header_labels(
      adj = "model name",
      estimate = "OR (95% CI)",
      levels = "case-control group",
      `0` = levels[["0"]],
      `1` = levels[["1"]]
    ) %>%
    add_header_row(
      values = c("", "N (% of total)", ""),
      colwidths = c(2, 2, 2)
    ) %>%
    merge_v(
      j = c("adj", "formula")
    ) %>%
    border_remove() %>%
    hline_bottom(
      j = 1:6,
      border = header_border
    ) %>%
    hline(
      border = header_border,
      j = 3:4,
      part = "header"
    ) %>%
    bg(
      i = ~!is.na(age),
      bg = "gray80"
    ) %>%
    # THIS DIDN'T WORK:
    # hline(
    #   border = header_border,
    #   i = ~ !is.na(age)
    # ) %>%
    hline_top(
      j = 1:6,
      border = header_border
    ) %>%
    border_inner_h(
      border = small_border,
      part = "body"
    ) %>%
    fix_border_issues() %>%
    fontsize(
      j = "formula",
      size = 8,
      part = "body"
    ) %>%
    # set_table_properties(layout = "autofit") %>%
    width(
      j = "formula",
      width = 5,
      unit = "cm"
    ) %>%
    width(
      j = "estimate",
      width = 4,
      unit = "cm"
    ) %>%
    set_caption(
      caption = subtitle
    )
}

make_table2_interact <- function(
    data,
    title,
    levels = list(
      `0`= "no lang.impair.",
      `1` = "lang.impair."
    )
){
  small_border <- officer::fp_border(color = "gray", width = 1)
  header_border <- officer::fp_border(color = "gray30", width = 2)
  
  as_grouped_data(
    data %>% 
      select(model, everything()) %>%
      mutate(
        aOR = ifelse(
          is.na(aOR),
          yes = "",
          no = sprintf(
            "%.2g (%.2g - %.2g)", aOR, conf.low, conf.high
          )
        ),
        p.val = ifelse(
          is.na(p.val),
          yes = "",
          no = sprintf(
            "%.2g", p.val
          )
        ),
        model = factor(
          model,
          levels = c("no_prs", "prs_manual", "prs_auto", "MTHFR_2strat"),
          labels = c("no interaction", "interaction with PRS (manual)",
                     "interaction with PRS (automatic)",
                     "interaction with rs1801133")
        )
      ) %>%
      arrange(model, age) %>%
      select(-conf.low, -conf.high),
    groups = c("model", "age")
  ) %>%
    as_flextable() %>% 
    bold(bold = TRUE, part = "header") %>%
    align(i = ~ !is.na(age) | !is.na(model), align = "center") %>% 
    bold(i = ~ !is.na(age) | !is.na(model)) %>%
    set_header_labels(
      label = "term",
      aOR = "OR (95% CI)",
      p.val = "p-value",
      levels = "case-control group",
      `0` = levels[["0"]],
      `1` = levels[["1"]]
    ) %>%
    add_header_row(
      values = c("", "N (% of total)", ""),
      colwidths = c(2, 2, 2)
    ) %>%
    border_remove() %>%
    hline_bottom(
      j = 1:6,
      border = header_border
    ) %>%
    hline(
      border = header_border,
      j = 3:4,
      part = "header"
    ) %>%
    bg(
      i = ~!is.na(age),
      bg = "gray60"
    ) %>%
    bg(
      i = ~!is.na(model),
      bg = "gray80"
    ) %>%
    # THIS DIDN'T WORK:
    # hline(
    #   border = header_border,
    #   i = ~ !is.na(age)
    # ) %>%
    hline_top(
      j = 1:6,
      border = header_border
    ) %>%
    border_inner_h(
      border = small_border,
      part = "body"
    ) %>%
    fix_border_issues() %>%
    width(
      j = "aOR",
      width = 4,
      unit = "cm"
    ) %>%
    set_caption(
      caption = title
    )
}

#' Creating nice table with results: showing impact of interaction term;
#'not stratified by use of ASM and epilepsy!
#' 
#' @param data - current estimates
#' @param subtitle - title for this table
#' @param levels - levels of the dependent variable
#' 
#' @return flextable that looks nice
#' 
make_table2_more_compact <- function(
    data,
    title,
    levels = list(
      `0`= "no lang.impair.",
      `1` = "lang.impair."
    )
){
  small_border <- officer::fp_border(color = "gray", width = 1)
  header_border <- officer::fp_border(color = "gray30", width = 2)
  
  as_grouped_data(
    data %>% 
      select(-label) %>%
      select(adj, everything()) %>%
      mutate(
        estimate = ifelse(
          is.na(estimate),
          yes = "",
          no = sprintf(
            "%.2f (%.2f - %.2f)", estimate, conf.low, conf.high
          )
        )
      ) %>%
      arrange(age, adj) %>%
      select(-conf.low, -conf.high, -formula) %>%
      select(levels, `0`, `1`, adj, estimate, age),
    groups = "age"
  ) %>%
    as_flextable() %>% 
    bold(bold = TRUE, part = "header") %>%
    align(i = ~ !is.na(age), align = "center") %>% 
    bold(i = ~ !is.na(age)) %>%
    set_header_labels(
      adj = "model name",
      estimate = "OR (95% CI)",
      levels = "case-control group",
      `0` = levels[["0"]],
      `1` = levels[["1"]]
    ) %>%
    add_header_row(
      values = c("", "N (% of total)", ""),
      colwidths = c(1, 2, 2)
    ) %>%
    merge_v(
      j = "adj"
    ) %>%
    merge_v(
      j = c("levels", "0", "1")
    ) %>%
    border_remove() %>%
    hline_bottom(
      j = 1:5,
      border = header_border
    ) %>%
    hline(
      border = header_border,
      j = 2:3,
      part = "header"
    ) %>%
    bg(
      i = ~!is.na(age),
      bg = "gray80"
    ) %>%
    # THIS DIDN'T WORK:
    # hline(
    #   border = header_border,
    #   i = ~ !is.na(age)
    # ) %>%
    hline_top(
      j = 1:5,
      border = header_border
    ) %>%
    border_inner_h(
      border = small_border,
      part = "body"
    ) %>%
    fix_border_issues() %>%
    width(
      j = "estimate",
      width = 4,
      unit = "cm"
    ) %>%
    set_caption(
      caption = title
    )
}

# plot nice estimates + CI
or_plot <- function(cur_data){
  cur_data %>%
    select(estimate, conf.high, conf.low) %>%
    filter(!is.na(estimate)) %>%
    ggplot(aes(estimate, 1)) +
    geom_vline(xintercept = 1, color = "grey80") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
    scale_x_continuous(breaks = c(0, 1, 2.5, 5, 7.5, 10)) +
    coord_cartesian(expand = FALSE, xlim = c(0, 10)) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank()
    )
}


#' Creating nice table with results: showing impact of interaction term;
#'not stratified by use of ASM and epilepsy!
#' 
#' @param data - current estimates
#' @param subtitle - title for this table
#' @param levels - levels of the dependent variable
#' 
#' @return flextable that looks nice, with both data and plots
#' 
make_table2_plot_compact <- function(
    data,
    title,
    levels = list(
      `0`= "no lang.impair.",
      `1` = "lang.impair."
    )
){
  small_border <- officer::fp_border(color = "gray", width = 1)
  header_border <- officer::fp_border(color = "gray30", width = 2)
  
  tmp_data_table <- as.data.table(
    data %>% 
      select(-label) %>%
      select(adj, everything()) %>%
      arrange(age, adj)
  )
  plots_per_age_adj <- tmp_data_table[
    , list(gg_plot = list(or_plot(.SD))), by = c("age", "adj")
  ]
  
  # where to add the plots: depends on the amount of age categories and
  #  adj categories
  n_age <- length(unique(plots_per_age_adj$age))
  n_adj <- length(unique(plots_per_age_adj$adj))
  idx_age_rows <- seq.int(from = 1, by = n_adj*2 + 1, length.out = n_age)
  idx_add_plot <- map(
    idx_age_rows,
    function(cur_age_row){
      seq.int(from = cur_age_row + 2, by = 2, length.out = n_adj)
    }
  ) %>% do.call(what = c, args = .)
  
  as_grouped_data(
    data %>% 
      select(-label) %>%
      select(adj, everything()) %>%
      arrange(age, adj) %>%
      mutate(
        estimate = ifelse(
          is.na(estimate),
          yes = "",
          no = sprintf(
            "%.2f (%.2f - %.2f)", estimate, conf.low, conf.high
          )
        )
      ) %>%
      select(-conf.low, -conf.high, -formula) %>%
      select(levels, `0`, `1`, adj, estimate, age),
    groups = "age"
  ) %>%
    as_flextable() %>% 
    bold(bold = TRUE, part = "header") %>%
    align(i = ~ !is.na(age), align = "center") %>% 
    bold(i = ~ !is.na(age)) %>%
    set_header_labels(
      adj = "model name",
      estimate = "OR (95% CI)",
      levels = "case-control group",
      `0` = levels[["0"]],
      `1` = levels[["1"]]
    ) %>%
    add_header_row(
      values = c("", "N (% of total)", ""),
      colwidths = c(1, 2, 2)
    ) %>%
    merge_v(
      j = "adj"
    ) %>%
    merge_v(
      j = c("levels", "0", "1")
    ) %>%
    border_remove() %>%
    hline_bottom(
      j = 1:5,
      border = header_border
    ) %>%
    hline(
      border = header_border,
      j = 2:3,
      part = "header"
    ) %>%
    bg(
      i = ~!is.na(age),
      bg = "gray80"
    ) %>%
    hline_top(
      j = 1:5,
      border = header_border
    ) %>%
    border_inner_h(
      border = small_border,
      part = "body"
    ) %>%
    fix_border_issues() %>%
    width(
      j = "estimate",
      width = 10,
      unit = "cm"
    ) %>%
    set_caption(
      caption = cur_title
    ) %>%
    append_chunks(
      i = idx_add_plot,
      j = "estimate",
      gg_chunk(
        value = plots_per_age_adj$gg_plot,
        width = 2,
        height = 0.8
      )
    )
}
