# DESCRIPTION: Check the data from Elisabeth; filter the women that have PRS
# AUTHOR: Julia ROmanowska
# DATE CREATED: 2022-04-01
# DATE MODIFIED: 2022-12-28

# SETUP ----
library(here)
library(tidyr)
library(readr)
library(haven)
library(dplyr)
library(ggplot2)
library(gtsummary)
# library(janitor)


#' @title Create the new covariate
#' @description The new covariate will be based on FOL_PRECONC, FOL2, FOL3, FOL4,
#'   FOL6, FOL7, and will include the difference in weeks from week 18 (when the
#'   blood was taken) to last week of reported folate intake.
#'   The FOL_* variables show whether any folic acid supplement was taken during:
#'   FOL_PRECONC - gw -4 to 0
#'   FOL2 - gw 0-4
#'   FOL3 - gw 5-8
#'   FOL4 - gw 9-12
#'   FOL6 - gw 13-16
#'   FOL7 - gw 17-20
#' @param data_in input data (tibble or data.frame)
#' @return tibble
#' @export
create_cov <- function(data_in){
  chosen_fol_vars <- c("FOL_PRECONC", "FOL2", "FOL3", "FOL4", "FOL6", "FOL7")
  out <- 
    data_in %>%
    select(M_ID_1111, tidyselect::all_of(chosen_fol_vars)) %>%
    group_by(M_ID_1111) %>%
    rowwise() %>%
    # collapse all the variables into one string of 0 and 1 (and NA if any)
    mutate(all_fol_str = paste0(
      c_across(tidyselect::all_of(chosen_fol_vars)),
      collapse = ""
    )
    ) %>%
    # calculate when was the last reported folate intake
    mutate(last_fol_intake = ifelse(
      length(
        stringr::str_locate_all(all_fol_str, "1")[[1]][,1]
      ) == 0,
      yes = 0,
      no = max(
        stringr::str_locate_all(all_fol_str, "1")[[1]][,1]
      )
    )
    ) %>%
    # calculate the difference with week 18
    mutate(
      diff_folate_intake = switch(
        as.character(last_fol_intake),
        "0" = 18,
        "1" = 18,
        "2" = 18 - 4,
        "3" = 18 - 8,
        "4" = 18 - 12,
        "5" = 18 - 16,
        "6" = 0
      )
    ) %>%
    ungroup()
  return(
    out %>%
      select(M_ID_1111, diff_folate_intake)
  )
}

# READ DATA ----
all_vars <- read_sav(
  here("DATA", "WORKFILE_v12_sprakogautisme_200921_withPRSonly_FINAL.sav")
  )
all_vars

prs_data <- read_delim(
  here("DATA", "PRS_MTHFR_allele_status_women_epilepsy.txt"),
  delim = "\t"
)
prs_data

# MERGE AND CHECK ----
merged_data_important_vars <- prs_data %>%
  left_join(
    all_vars %>%
      select(
        PREG_ID_1111, M_ID_1111, IK_EPI, # IDs and indicator of epilepsy
        EPILEPSY2_c, AED_PREG_c, NOAED_PREG_c, # epilepsy diagnosis and use of AEDs
        ALDER_MOR, # mother's age (at birth?)
        SVLEN, # length of pregnancy (in weeks)
        PARITET_MFR, # parity
        KJONN_DIKOTOM, # sex
        PLURAL, # multiple vs. singleton pregnancy
        DUPLIKAT, # whether the child has a sibling in the dataset
        UNPLANNED_PREG, # unplanned
        ROYK_SVSKAP, # smoking during pregnancy
        ALC_PREG, # alcohol during 1st trimester
        BMIbefore, # BMI before pregnancy
        SCL1_depr22, # depression/anxiety during pregnancy
        SES_sum, # socio-economic status (0-3)
        AT_PREG_TOT, # epileptic seizures during pregnancy
        andel_med_GTK2, # tonic-clonic seizures during pregnancy
        AED_POLY, VPA_MONO, CBZ_MONO, LTG_MONO, LTA_MONO, TPX_MONO, OXC_MONO, #use of ASM
        AED_kons_sum_mor, AED_kons_sum, # maternal and maternal + umb.cord
        FOL_PRE_or_TRIM1, # folic acid suppl.periconceptionally
        VIT_FOL_PRETRIM1TRIM2, # folic acid suppl. use during gest.weeks -4 through 20
        FOLAT, # level of intake of folic acid through diet (mg/day)
        s_folat, # level of intake of folic acid through supplements (ug/day)
        FOL_hydro, # folate concentration in blood
        lang_imp_18mn = ASQ18_LANGUAGE_1.5SD, # language impairment at 1.5 years
        lang_imp_36mn = dårlig_språk_36, # language impairment at 3 years
        lang_imp_5yr = total_lang_delay_5år_i, # language impairment at 5 years
        lang_imp_8yr = SP_OTT_SEMAN_8ÅR_SPRÅKVANSKER_verC_i, # language impairment at 8 years
        autistic_trait_3yr = SCQ39M11_i, # autistic traits at age 3 years
        autistic_trait_8yr = SCQ39M11_8Y_i # autistic traits at age 8 years
        ),
    by = "M_ID_1111"
  ) %>%
  mutate(
    across(
      c(a1:CC, IK_EPI:NOAED_PREG_c, PARITET_MFR:ALC_PREG,
        SCL1_depr22:OXC_MONO, FOL_PRE_or_TRIM1:VIT_FOL_PRETRIM1TRIM2,
        lang_imp_18mn:autistic_trait_8yr),
      ~ as.factor(.x))
  )
skimr::skim(merged_data_important_vars)

# find duplicates
dim(merged_data_important_vars)
length(unique(merged_data_important_vars$PREG_ID_1111))

# there are some missing values in PREG_ID_1111 and some duplicates
merged_data_important_vars %>%
  filter(is.na(PREG_ID_1111)) %>% glimpse()
# these were in the PRS file, but not in the 'main file'
merged_data_important_vars <- merged_data_important_vars %>%
  filter(!is.na(PREG_ID_1111))

dim(merged_data_important_vars)
length(unique(merged_data_important_vars$PREG_ID_1111))

# Elisabeth said that the variable "CTR" has surely all those with epilepsy
table(
  (all_vars %>%
     filter(PREG_ID_1111 %in% merged_data_important_vars$PREG_ID_1111) %>%
     filter(is.na(EPILEPSY2_c))
   )$CTR)
# thus, these with missing i EPILEPSY2_c are actually controls?
# 2022-04-25: ELisabeth said that no - only those with EPILEPSY2_c non-missing!

# Elisabeth said that NOAED_PREG_c and AED_PREG_c are part of one variable!
merged_data_important_vars %>% 
  filter(PLURAL == 1) %>%
  janitor::tabyl(NOAED_PREG_c, AED_PREG_c)

merged_data_important_vars %>% 
  filter(PLURAL == 1) %>%
  janitor::tabyl(EPILEPSY2_c, AED_PREG_c)

# filter those we don't want
merged_data_important_vars %>%
  count(IK_EPI)
merged_data_important_vars %>%
  count(PLURAL)
# 2022-06-08: Elisabeth says that we don't need to exclude those that don't have PLURAL var!
merged_data_important_vars %>%
  count(EPILEPSY2_c)

dim(merged_data_important_vars)
merged_data_important_vars <- merged_data_important_vars %>%
  filter(IK_EPI == 0) %>% # only verified epilepsy cases
  # filter(PLURAL == 1) %>% # only singular pregnancies
  filter(!is.na(EPILEPSY2_c)) %>%
  mutate(
    CC_STATUS = factor( # merge the two mutual variables
      ifelse(
        is.na(as.character(NOAED_PREG_c)),
        yes = "epi_AED", #women with epilepsy and taking AED
        no = ifelse(
          NOAED_PREG_c == 1,
          yes = "epi_no_AED", #women with epilepsy and not taking AED
          no = "no_epi" #women without epilepsy
        )
      )
    ),
    .before = EPILEPSY2_c
  )
dim(merged_data_important_vars)

skimr::skim(merged_data_important_vars)

## checking all covariates ----
# pregnancy length (SVLEN)
merged_data_important_vars %>%
  ggplot(aes(CC_STATUS, SVLEN)) +
  geom_boxplot(varwidth = TRUE)

#AED_PREG + prs_auto_dec + ALC_PREG + 
  #     ROYK_SVSKAP + SCL1_depr22 + ALDER_MOR + PARITET_MFR + SVLEN + BMIbefore + 
  #     SES_sum + FOLAT + s_folat + VIT_FOL_PRETRIM1TRIM2 + FOL_PRE_or_TRIM1 + 
  #     KJONN_DIKOTOM

# PRS - automatically generated
merged_data_important_vars %>%
  ggplot(aes(CC_STATUS, prs_auto)) +
  geom_boxplot(varwidth = TRUE)

# PRS - manually calculated
merged_data_important_vars %>%
  ggplot(aes(CC_STATUS, prs_manual)) +
  geom_boxplot(varwidth = TRUE)

# smoking (ROYK_SVSKAP)
merged_data_important_vars %>%
  janitor::tabyl(CC_STATUS, ROYK_SVSKAP)

# mother's age
merged_data_important_vars %>%
  ggplot(aes(CC_STATUS, ALDER_MOR)) +
  geom_boxplot(varwidth = TRUE)

# parity (PARITET_MFR)
merged_data_important_vars %>%
  janitor::tabyl(CC_STATUS, PARITET_MFR)

# BMIbefore
merged_data_important_vars %>%
  ggplot(aes(CC_STATUS, BMIbefore)) +
  geom_boxplot(varwidth = TRUE)

# FOLAT (folate intake from diet)
merged_data_important_vars %>%
  ggplot(aes(CC_STATUS, FOLAT)) +
  geom_boxplot(varwidth = TRUE)

# s_folat (folate intake from supplements)
merged_data_important_vars %>%
  ggplot(aes(CC_STATUS, s_folat)) +
  geom_boxplot(varwidth = TRUE)

# socioeconomical status:
merged_data_important_vars %>%
  janitor::tabyl(CC_STATUS, SES_sum)

# folate intake and vitamins periconceptionally
merged_data_important_vars %>%
  janitor::tabyl(CC_STATUS, VIT_FOL_PRETRIM1TRIM2)

merged_data_important_vars %>%
  janitor::tabyl(CC_STATUS, FOL_PRE_or_TRIM1)

# baby's sex
merged_data_important_vars %>%
  janitor::tabyl(CC_STATUS, KJONN_DIKOTOM)

# how many sibling pairs?
merged_data_important_vars %>%
  janitor::tabyl(DUPLIKAT)

merged_data_important_vars %>%
  janitor::tabyl(CC_STATUS, DUPLIKAT)

siblings <- merged_data_important_vars %>%
  filter(DUPLIKAT == 1)
famillies <- merged_data_important_vars %>%
  filter(M_ID_1111 %in% siblings$M_ID_1111) %>%
  count(M_ID_1111) %>%
  count(n) %>%
  arrange(nn)
famillies

# CREATE THE NEW COVARIATE ----
# this is the same as for the control group, showing how long time from the
#   last folate supplement intake to the week 18, where the blood was drawn
chosen_fol_vars <- c("FOL_PRECONC", "FOL2", "FOL3", "FOL4", "FOL6", "FOL7")
fol_suppl_details <- all_vars %>%
  select(M_ID_1111, tidyselect::all_of(chosen_fol_vars)) %>%
  filter(M_ID_1111 %in% merged_data_important_vars$M_ID_1111)
# some of the mothers did not answer all the questionnaires, so there
#   are missing variables:
fol_suppl_details %>%
  filter(if_any(FOL_PRECONC:FOL7, ~is.na(.x)))
# I will assume that these mothers had '0' instead of 'NA'
replace_NA_vec <- as.list(rep(0, length(chosen_fol_vars)))
names(replace_NA_vec) <- chosen_fol_vars

fol_suppl_details <- fol_suppl_details %>%
  tidyr::replace_na(replace_NA_vec)

fol_suppl_diff <- create_cov(fol_suppl_details)
fol_suppl_diff

# SAVE CLEANED DATA ----
write_delim(
  merged_data_important_vars,
  here("DATA", "cleaned_data.txt"),
  delim = "\t"
)

write_delim(
  fol_suppl_diff,
  here("DATA", "folate_suppl_diff_covariate.txt"),
  delim = "\t"
)

# DISTRIBUTIONS ----
ggplot(merged_data_important_vars, aes(prs_auto)) +
  geom_histogram() +
  facet_wrap(facets = vars(EPILEPSY2_c), scales = "free_y") +
  labs(
    title = "Histogram of PRS values for women with epilepsy (1) and without (0)",
    subtitle = "(best-fit PRS, 2 SNPs)"
    )

ggplot(merged_data_important_vars, aes(prs_manual)) +
  geom_histogram() +
  facet_wrap(facets = vars(EPILEPSY2_c), scales = "free_y") +
  labs(
    title = "Histogram of PRS values for women with epilepsy (1) and without (0)",
    subtitle = "(2nd best-fit PRS, 76 SNPs)"
    )

# CHECKING FOLATE CONC. stratified by MTHFR SNP ----
plot_folate_status_strat_genotype_suppl <- merged_data_important_vars %>%
  select(iid, minor_al_load, CC, M_ID_1111, EPILEPSY2_c, FOL_hydro,
         FOL_PRE_or_TRIM1:FOL_hydro) %>%
  distinct() %>%
  mutate(s_folat_cat = cut(
    s_folat/1000,
    breaks = c(0, 0.4, 1, max(s_folat/1000, na.rm = TRUE) + 1),
    labels = c("low", "medium", "high"),
    include.lowest = TRUE)
  ) %>%
  ggplot() +
  aes(minor_al_load, FOL_hydro) +
  geom_boxplot(aes(color = s_folat_cat), varwidth = TRUE) +
  scale_color_discrete("supplement\n intake dose") +
  labs(
    title = "Folate conc. in blood measured at GW 18 and maternal rs1801133 genotype,",
    subtitle = "stratified by dose of folate supplement",
    caption = "width of the boxplot is proportional to sample size;\n low dose = <0.4 mg/day, medium = 0.4-1 mg/day, high = >1 mg/day"
  ) +
  xlab("rs1801133 genotype") +
  ylab("maternal plasma folate conc.(nmol/L)")
plot_folate_status_strat_genotype_suppl

ggsave(
  plot = plot_folate_status_strat_genotype_suppl,
  filename = here(
    "FIGURES", "folate_status_stratified_genotype_and_suppl_dose.png"
  )
)

# strange high values of FOL_hydro for those with TT-genotype and low suppl.dose
merged_data_important_vars %>%
  filter(s_folat < 400 & minor_al_load == "TT" & !is.na(FOL_hydro)) %>%
  skimr::skim()
# that's only one value!

merged_data_important_vars %>%
  filter(s_folat > 1000 & minor_al_load == "TT" & !is.na(FOL_hydro)) %>%
  skimr::skim()
# and only 3 persons who had TT and high doses of supplement

plot_folate_status_strat_genotype_suppl +
  facet_wrap(vars(FOL_PRE_or_TRIM1), labeller = label_both) +
  labs(subtitle = "stratified by dose of folate supplement and whether they took any folate preconceptionally")

ggsave(
  filename = here(
    "FIGURES", "folate_status_stratified_genotype_and_suppl_dose_and_intake_preconc.png"
  )
)

# DESCRIPTIVE STATS ----
# some women were taking ASM which was not in one of the groups that are
#   in the *_MONO variables, so they were not included in the tables!
data_for_descr_tbl <- merged_data_important_vars %>%
  select(CC_STATUS, AED_PREG_c, AED_POLY, ends_with("_MONO"), s_folat,
         FOLAT, FOL_hydro, FOL_PRE_or_TRIM1) %>%
  pivot_longer(cols = ends_with("_MONO"),
               names_to = "AED_MONO_name",
               values_to = "AED_MONO_val")

# were there any who got many ASM_MONO?
data_for_descr_tbl %>%
  janitor::tabyl(AED_MONO_name, AED_MONO_val) %>%
  janitor::adorn_totals(where = "col")
# no :)

# to add the "other" category, we need to do this at an earlier stage
data_for_descr_tbl <- merged_data_important_vars %>%
  select(CC_STATUS, AED_PREG_c, AED_POLY, ends_with("_MONO"), s_folat,
         FOLAT, FOL_hydro, FOL_PRE_or_TRIM1, minor_al_load) %>%
  # check who is already in any "_MONO" category
  rowwise() %>%
  mutate(any_mono_ASM = any(c_across(ends_with("MONO")) == 1)) %>%
  ungroup() %>%
  mutate( 
    # the AED_PREG_c variable does not include values for women who had
    #   epilepsy but did not use ASM
    AED_PREG_c = if_else(
      is.na(AED_PREG_c),
      true = 0,
      false = as.numeric(levels(AED_PREG_c)[AED_PREG_c])
    )
  ) %>%
  # those who used a monotherapy but were not specified in any other '_MONO'
  mutate(
    other_MONO = as.factor(if_else(
      AED_PREG_c == 1 & any_mono_ASM == 0 & AED_POLY == 0,
      true = 1,
      false = 0))
  ) %>% select(-AED_PREG_c, -any_mono_ASM)

data_for_descr_tbl %>%
  pivot_longer(cols = ends_with("_MONO"),
               names_to = "AED_MONO_name",
               values_to = "AED_MONO_val") %>%
  janitor::tabyl(AED_MONO_name, AED_MONO_val) %>%
  janitor::adorn_totals(where = "col")


descriptive_tbl <- tbl_summary(
  data_for_descr_tbl %>%
    select(-minor_al_load) %>%
    mutate(
      CC_STATUS = forcats::fct_recode(
        CC_STATUS,
        "Children of mothers without epilepsy" = 'no_epi',
        "AED-exposed children of mothers with epilepsy" = 'epi_AED',
        "AED-unexposed children of mothers with epilepsy" = 'epi_no_AED'
      )) %>%
    mutate(across(
      c(AED_POLY:OXC_MONO, other_MONO, FOL_PRE_or_TRIM1), 
      ~ forcats::fct_recode(
        .x, "No" = '0', "Yes" = '1'
      )
    )),
  by = CC_STATUS,
  label = list(
      AED_POLY ~ "ASM polytherapy",
      s_folat ~ "folate intake from supplements (mg)",
      FOLAT ~ "folate intake from diet",
      FOL_hydro ~ "folate concentration in blood",
      FOL_PRE_or_TRIM1 ~ "use of folate before or in 1st trim."
    ),
    missing_text = "(missing)"
  )

descriptive_tbl_gt <- descriptive_tbl %>%
  as_gt(rowname_col = "label") %>%
  gt::tab_row_group(
    label = "Perinatal factors",
    rows = variable %in% c("s_folat", "FOLAT", "FOL_PRE_or_TRIM1")
  ) %>%
  gt::tab_row_group(
    label = "Use of AED",
    rows = ((variable == "AED_POLY") | stringr::str_ends(variable, "_MONO"))
  )

gt::gtsave(
  data = descriptive_tbl_gt,
  filename = "Table01_descriptive_stats.html",
  path = here("RESULTS")
)

writeLines(
  text = gt::as_rtf(data = descriptive_tbl_gt),
  con = here("RESULTS", "table_1_descriptive.rtf")
)

# DESCRIPTIVE TABLE - stratified by MTHFR-snp ----
mthfr_descriptive_tbl <- tbl_summary(
  data_for_descr_tbl %>%
    mutate(MTHFR_2strat = as.factor(if_else(
      minor_al_load == "CC",
      "CC",
      "CT/TT"
    ))) %>%
    mutate(
      CC_STATUS = forcats::fct_recode(
        CC_STATUS,
        "Children of mothers without epilepsy" = 'no_epi',
        "AED-exposed children of mothers with epilepsy" = 'epi_AED',
        "AED-unexposed children of mothers with epilepsy" = 'epi_no_AED'
      )) %>%
    mutate(across(
      c(AED_POLY:OXC_MONO, other_MONO, FOL_PRE_or_TRIM1), 
      ~ forcats::fct_recode(
        .x, "No" = '0', "Yes" = '1'
      )
    )) %>%
    select(-minor_al_load),
  by = MTHFR_2strat,
  label = list(
    CC_STATUS ~ "case-control status",
    AED_POLY ~ "ASM polytherapy",
    s_folat ~ "folate intake from supplements (mg)",
    FOLAT ~ "folate intake from diet",
    FOL_hydro ~ "maternal folate concentration in blood",
    FOL_PRE_or_TRIM1 ~ "use of folate before or in 1st trim."),
  missing_text = "(missing)"
)

mthfr_descriptive_tbl_gt <- mthfr_descriptive_tbl %>%
  as_gt(rowname_col = "label") %>%
  gt::tab_row_group(
    label = "Perinatal factors",
    rows = variable %in% c("s_folat", "FOLAT", "FOL_PRE_or_TRIM1", "FOL_hydro")
  ) %>%
  gt::tab_row_group(
    label = "Use of AED",
    rows = ((variable == "AED_POLY") | stringr::str_ends(variable, "_MONO"))
  ) %>%
  gt::tab_header(
    title = "Characteristics stratified by SNP rs1801133"
  )

gt::gtsave(
  data = mthfr_descriptive_tbl_gt,
  filename = "Table01_MTHFR_SNP_descriptive_stats.html",
  path = here("RESULTS")
)

writeLines(
  text = gt::as_rtf(data = mthfr_descriptive_tbl_gt),
  con = here("RESULTS", "table_1_MTHFR_SNP_descriptive.rtf")
)


# Checking the siblings ----

# how many times various number of siblings?
merged_data_important_vars %>%
  group_by(M_ID_1111) %>%
  count(DUPLIKAT) %>%
  ungroup() %>%
  count(DUPLIKAT, n)

# the DUPLIKAT == 1 & n == 1 is strange - maybe there was a sibling, but
#   the data was not good enough for the study here? Or these can be half-
#   siblings... but how to account for that? I don't know which ones are
#   related
merged_data_important_vars %>%
  group_by(M_ID_1111) %>%
  count(DUPLIKAT) %>%
  ungroup()  %>%
  filter(DUPLIKAT == 1 & n == 1)
