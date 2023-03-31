# DESCRIPTION: Are the groups of women comparable between the different
#   age groups (i.e., the age of the child when the questionnaire was returned)
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-08-15
# DATE MODIFIED: 2022-08-16

# SETUP ----
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(sandwich)
library(lmtest)
library(here)

# READ DATA ----
input_data <- read_delim(
  here("DATA", "data_ready_analysis_table2.txt"),
  delim = "\t"
) %>%
  mutate(
    across(
      c(minor_al_load, a1, a2, CC, IK_EPI, EPILEPSY2_c, PARITET_MFR:ALC_PREG, 
        SCL1_depr22, AED_POLY, ends_with("_MONO"), FOL_PRE_or_TRIM1,
        VIT_FOL_PRETRIM1TRIM2, starts_with("lang_imp"), starts_with("autistic"),
        AED_PREG_c, NOAED_PREG_c, CC_STATUS, SES_sum, BMIbefore_cat,
        ALDER_MOR_cat, SVLEN_cat),
      ~ factor(.x)
    )
  )
all_ages <- levels(input_data$lang_imp_age)


# ANALYZE ----
input_data %>%
  group_by(lang_imp_age) %>%
  summarise(lang_imp_missing = sum(is.na(lang_imp)))

## general epilepsy status vs. age ----
cur_data <- input_data %>%
  mutate(CC = as.factor(ifelse(
    CC_STATUS == "no_epi",
    yes = CC_STATUS,
    no = "epi"
  ))) %>%
  select(CC, CC_STATUS, m_id, lang_imp_age, lang_imp) %>%
  filter(!is.na(lang_imp))
epi_status_age_OR <- glm(
  CC ~ lang_imp_age,
  family = "binomial",
  data = cur_data
)
robust_ci <- as_tibble(
  exp(coefci(
    epi_status_age_OR,
    vcov. = sandwich::vcovCL(epi_status_age_OR, cluster = cur_data$m_id)
  )),
  rownames = "term"
)
colnames(robust_ci) <- c("term", "conf.low", "conf.high")

tidy_res_epi_status <- broom::tidy(epi_status_age_OR, exp = TRUE) %>%
  left_join(robust_ci, by = "term")
tidy_res_epi_status

## women who took AED during pregnancy vs controls ----
cur_data <- input_data %>%
  filter(CC_STATUS != "epi_no_AED") %>%
  mutate(CC_STATUS = forcats::fct_drop(CC_STATUS)) %>%
  mutate(CC_STATUS = forcats::fct_relevel(CC_STATUS, "no_epi", "epi_AED")) %>%
  select(CC_STATUS, m_id, lang_imp_age, lang_imp) %>%
  filter(!is.na(lang_imp))
janitor::tabyl(cur_data, CC_STATUS, lang_imp_age)

epi_aed_age_OR <- glm(
  CC_STATUS ~ lang_imp_age,
  family = "binomial",
  data = cur_data
)
robust_ci <- as_tibble(
  exp(coefci(
    epi_aed_age_OR,
    vcov. = sandwich::vcovCL(epi_aed_age_OR, cluster = cur_data$m_id)
  )),
  rownames = "term"
)
colnames(robust_ci) <- c("term", "conf.low", "conf.high")

tidy_res_epi_aed <- broom::tidy(epi_aed_age_OR, exp = TRUE) %>%
  left_join(robust_ci, by = "term")
tidy_res_epi_aed

## women who did not take AED during pregnancy vs. controls ----
cur_data <- input_data %>%
  filter(CC_STATUS != "epi_AED") %>%
  mutate(CC_STATUS = forcats::fct_drop(CC_STATUS)) %>%
  mutate(CC_STATUS = forcats::fct_relevel(CC_STATUS, "no_epi", "epi_no_AED")) %>%
  select(CC_STATUS, m_id, lang_imp_age, lang_imp) %>%
  filter(!is.na(lang_imp))
janitor::tabyl(cur_data, CC_STATUS, lang_imp_age)

epi_no_aed_age_OR <- glm(
  CC_STATUS ~ lang_imp_age,
  family = "binomial",
  data = cur_data
)
robust_ci <- as_tibble(
  exp(coefci(
    epi_no_aed_age_OR,
    vcov. = sandwich::vcovCL(epi_no_aed_age_OR, cluster = cur_data$m_id)
  )),
  rownames = "term"
)
colnames(robust_ci) <- c("term", "conf.low", "conf.high")

tidy_res_epi_status <- broom::tidy(epi_no_aed_age_OR, exp = TRUE) %>%
  left_join(robust_ci, by = "term")
tidy_res_epi_status

## language impairment (disregarding the epilepsy status) vs age ----
janitor::tabyl(input_data, lang_imp_age, lang_imp) %>%
  mutate(case_ratio = `1`/`0`)

lang_imp_OR <- glm(
  lang_imp ~lang_imp_age,
  family = "binomial",
  data = input_data
)
robust_ci <- as_tibble(
  exp(coefci(
    lang_imp_OR,
    vcov. = sandwich::vcovCL(lang_imp_OR, cluster = input_data$m_id)
  )),
  rownames = "term"
)
colnames(robust_ci) <- c("term", "conf.low", "conf.high")

tidy_res_lang_imp <- broom::tidy(lang_imp_OR, exp = TRUE) %>%
  left_join(robust_ci, by = "term")
tidy_res_lang_imp

ggplot(tidy_res_lang_imp, aes(term, estimate)) +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high)
  )
