# DESCRIPTION: Analysis, table 2: PRS for low folate concentration and risk of
#    language impairment in the child; checking relations among variables;
#    calculating crude models and with each variable separately, to check
#    its impact on the crude estimate
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-04-21
# DATE MODIFIED: 2022-12-22

# SETUP ----
library(here)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(patchwork)
library(finalfit)

# READ DATA ----
input_data <- read_delim(
  here("DATA", "cleaned_data.txt"),
  delim = "\t"
  ) %>%
  select(iid:NOAED_PREG_c, PLURAL, DUPLIKAT, FOL_PRE_or_TRIM1,
         FOLAT:autistic_trait_8yr) %>%
  mutate(
    across(
      c(minor_al_load, a1, a2, CC, IK_EPI, EPILEPSY2_c, FOL_PRE_or_TRIM1,
        starts_with("lang_imp"), starts_with("autistic"),
        AED_PREG_c, NOAED_PREG_c, CC_STATUS),
      ~ factor(.x)
    )
  )
skimr::skim(input_data)

# LANGUAGE IMPAIRMENT DATA ----
# how does the distribution of language impariment scores look like?
# how does it relate to use of ASM?
input_data %>%
  janitor::tabyl(lang_imp_18mn, CC_STATUS)

input_data %>%
  janitor::tabyl(lang_imp_36mn, CC_STATUS)

input_data %>%
  janitor::tabyl(lang_imp_5yr, CC_STATUS)

input_data %>%
  janitor::tabyl(lang_imp_8yr, CC_STATUS)

# I need to first aggregate the measurements
## pivot data ----
input_data_long <- input_data %>%
  pivot_longer(
    cols = starts_with("lang_imp"),
    names_to = "lang_imp_age",
    values_to = "lang_imp"
  ) %>%
  mutate(
    lang_imp_age = factor(
      stringr::word(
        lang_imp_age,
        start = -1,
        end = -1,
        sep = stringr::fixed("_")
      )
    )
  )

input_data_long %>%
  janitor::tabyl(lang_imp, CC_STATUS, lang_imp_age)

input_data_long %>%
  filter(EPILEPSY2_c == "1") %>%
  janitor::tabyl(lang_imp, CC_STATUS, lang_imp_age)

ggplot(input_data_long, aes(prs_auto)) +
  geom_histogram(aes(fill = lang_imp), alpha = 0.7) +
  facet_wrap(facets = vars(lang_imp_age)) +
  theme_light()

ggplot(input_data_long, aes(prs_manual)) +
  geom_histogram(aes(fill = lang_imp), alpha = 0.7) +
  facet_wrap(facets = vars(lang_imp_age)) +
  theme_light()

## COAVRIATES ----
# we choose only those that give use info about folate

# CHECKING STRATIFICATION BY MTHFR ALLELE ----
input_data_long %>%
  filter(!is.na(lang_imp)) %>%
  janitor::tabyl(minor_al_load, EPILEPSY2_c, lang_imp_age)

# we should have 2 strata, instead of 3
input_data_long <- input_data_long %>%
  mutate(MTHFR_2strat = ifelse(
    minor_al_load %in% c("CT", "TT"),
    yes = "CT/TT",
    no = "CC"
  ))

input_data_long %>%
  filter(!is.na(lang_imp)) %>%
  janitor::tabyl(MTHFR_2strat, EPILEPSY2_c, lang_imp_age)


# entire dataset, with interactions ----
# I will take an entire dataset, but now case/control status will be only
# epilepsy/no epilepsy (EPILEPSY2_c); then, PRS will interact with the var
# that indicates taking AED (independent from CC_STATUS) - this one I need
# to create
input_data_long <- input_data_long %>%
  mutate(any_AED = ifelse(
    CC_STATUS == "epi_AED",
    yes = TRUE,
    no = FALSE
  ))

input_data_long <- input_data_long %>%
  mutate(m_id = as.numeric(stringr::str_sub(M_ID_1111, 2, -1)))

cur_data <- input_data_long
janitor::tabyl(cur_data, CC_STATUS, any_AED, lang_imp_age)
janitor::tabyl(cur_data, EPILEPSY2_c, CC_STATUS,  lang_imp_age)

## SAVE data ready for analysis ----
write_delim(
  input_data_long,
  file = here("DATA", "data_ready_analysis_table2.txt"),
  delim = "\t"
)
