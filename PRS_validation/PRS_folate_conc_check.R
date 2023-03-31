# DESCRIPTION: Create figure 2 - PRS vs folate conc.
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-12-01
# DATE MODIFIED: 2022-12-28

# SETUP ----
library(targets)
library(readr)
library(dplyr)
library(ggplot2)
library(here)
library(patchwork)

# READ DATA ----
out_prefix <- "RESULTS/MoBa_folat_prsice_with-cov_selected"
mthfr_status_file <- here("RESULTS", "mthfr_status")

# best-best PRS: only 2 SNPs
tar_load(best_prs_fit)

# 2nd-best PRS: 76 SNPs
tar_load(best_prs_fit_manual)
# best_prs_fit_manual

# covariates:
tar_load(cntrl_ids_folate_cov)
tar_load(prs_covariate)

# here, there are also 'FOLAT' and 's_folat':
cntrl_ids_extra_folate_cov <- read_table(here("DATA", "folate_cov_merged.txt"))
cntrl_ids_folate_cov <- cntrl_ids_folate_cov %>%
  left_join(
    cntrl_ids_extra_folate_cov %>%
      select(M_ID_1111, SENTRIXID, FOLAT, s_folat)
  ) %>% left_join(
    prs_covariate %>%
      select(IID, PC1, PC2, PC3),
    by = c("SENTRIXID" = "IID")
  )

# check the output from PLINK recoding
mthfr_status <- read_delim(
  paste0(mthfr_status_file, ".ped"),
  delim = " ",
  col_names = c("fid", "iid", "mother", "father", "sex", "cc", "a1", "a2")
) %>%
  select(iid, a1, a2)
# mthfr_status


# PREPARE DATA ----
## SNP
mthfr_status <- mthfr_status %>%
  mutate(
    # how many doses of minor allele in this SNP?
    minor_al_load = case_when(
      a1 + a2 == 2 ~ 2,
      a1 + a2 == 3 ~ 1,
      a1 + a2 == 4 ~ 0
    )) %>%
  # strata for analyses
  mutate(stratum = if_else(
    minor_al_load == 0,
    "CC",
    "CT/TT"
  ))
mthfr_status %>%
  ggplot(aes(minor_al_load)) +
  geom_bar()

mthfr_status %>%
  count(a1)
mthfr_status %>%
  count(a2)
mthfr_status %>%
  count(a2) %>%
  mutate(prcnt = n/nrow(mthfr_status)*100)

mthfr_status %>%
  count(minor_al_load) %>%
  mutate(prcnt = n/nrow(mthfr_status)*100)

mthfr_status %>%
  janitor::tabyl(stratum)

mthfr_status_folate <- mthfr_status %>%
  left_join(
    cntrl_ids_folate_cov,
    by = c("iid" = "SENTRIXID")
  ) %>%
  filter(!is.na(SPFOLATE_F))
# mthfr_status_folate

ggplot(mthfr_status_folate, aes(SPFOLATE_F)) +
  geom_density(aes(fill = stratum), alpha = 0.5) +
  theme_light()

ggplot(mthfr_status_folate, aes(SPFOLATE_F)) +
  geom_density(aes(fill = as.factor(minor_al_load)), alpha = 0.5) +
  theme_light()

mthfr_status_folate %>%
  janitor::tabyl(stratum)

## PRS
best_prs_fit <- best_prs_fit %>%
  select(-FID) %>%
  left_join(
    cntrl_ids_folate_cov,
    by = c("IID" = "SENTRIXID")
  )

best_prs_fit_all <- best_prs_fit
best_prs_fit <- best_prs_fit %>%
  filter(!is.na(SPFOLATE_F))


best_prs_fit_manual <- best_prs_fit_manual %>%
  select(-FID) %>%
  left_join(
    cntrl_ids_folate_cov,
    by = c("IID" = "SENTRIXID")
  )

best_prs_fit_manual_all <- best_prs_fit_manual
best_prs_fit_manual <- best_prs_fit_manual %>%
  filter(!is.na(SPFOLATE_F))

# CORRELATION AND PLOT ----

## PRS - automatically chosen:
cor_test_prs_auto <- correlation::cor_test(
  data = best_prs_fit,
  x = "PRS",
  y = "orig_SPFOLATE_F",
  method = "spearman"
)
cor_test_prs_auto

# partial - i.e., adjusted for other (linear) covariates
cor_test_prs_auto_adj <- correlation::cor_test(
  data = best_prs_fit %>%
    select(PRS, orig_SPFOLATE_F, diff_folate_intake, FOLAT, s_folat),
  x = "PRS",
  y = "orig_SPFOLATE_F",
  method = "spearman",
  partial = TRUE
)
cor_test_prs_auto_adj

# partial - i.e., adjusted for other (linear) covariates - nowwith genetic PC
cor_test_prs_auto_adj_pc <- correlation::cor_test(
  data = best_prs_fit %>%
    select(PRS, orig_SPFOLATE_F, diff_folate_intake, PC1:PC2),
  x = "PRS",
  y = "orig_SPFOLATE_F",
  method = "spearman",
  partial = TRUE
)
cor_test_prs_auto_adj_pc

prs_auto_cor_plot <- plot(
  cor_test_prs_auto,
  show_statistic = FALSE, show_ci = FALSE, stars = TRUE) +
  # annotate("text", x = -0.001, y = 50, label = paste0("N = ", cor_test_prs$n_Obs)) +
  xlab("polygenic risk score (PRS)") +
  ylab("maternal plasma folate conc. (nmol/L)") +
  theme_minimal()


## PRS - manually chosen: ----
cor_test_prs <- correlation::cor_test(
  data = best_prs_fit_manual,
  x = "PRS",
  y = "orig_SPFOLATE_F",
  method = "spearman"
)
cor_test_prs

# partial - i.e., adjusted for other (linear) covariates
cor_test_prs_adj <- correlation::cor_test(
  data = best_prs_fit_manual %>%
    select(PRS, orig_SPFOLATE_F, diff_folate_intake, FOLAT, s_folat),
  x = "PRS",
  y = "orig_SPFOLATE_F",
  method = "spearman",
  partial = TRUE
)
cor_test_prs_adj

# partial - i.e., adjusted for other (linear) covariates - nowwith genetic PC
cor_test_prs_adj_pc <- correlation::cor_test(
  data = best_prs_fit_manual %>%
    select(PRS, orig_SPFOLATE_F, diff_folate_intake, PC1:PC2),
  x = "PRS",
  y = "orig_SPFOLATE_F",
  method = "spearman",
  partial = TRUE
)
cor_test_prs_adj_pc


prs_cor_plot <- plot(
  cor_test_prs,
  show_statistic = FALSE, show_ci = FALSE, stars = TRUE) +
  # annotate("text", x = -0.001, y = 50, label = paste0("N = ", cor_test_prs$n_Obs)) +
  xlab("polygenic risk score (PRS)") +
  ylab("maternal plasma folate conc. (nmol/L)") +
  theme_minimal()

# density vs MTHFR-status
mthfr_strat_plot <- ggplot(mthfr_status_folate, aes(x = stratum, y = orig_SPFOLATE_F)) + 
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_minimal() +
  xlab("rs1801133 stratum") +
  ylab("maternal plasma folate conc. (nmol/L)")

(prs_cor_plot / mthfr_strat_plot) +
  plot_annotation(tag_levels = 'A', tag_suffix = ")")

ggsave(here("FIGURES", "prs_manual_mthfr_vs_folate_conc_cntrls.png"))

(prs_auto_cor_plot / mthfr_strat_plot) +
  plot_annotation(tag_levels = 'A', tag_suffix = ")")

ggsave(here("FIGURES", "prs_auto_mthfr_vs_folate_conc_cntrls.png"))

#median values for the two categories:
summary(mthfr_status_folate %>%
          filter(stratum == "CC") %>%
          pull(orig_SPFOLATE_F))

summary(mthfr_status_folate %>%
          filter(stratum == "CT/TT") %>%
          pull(orig_SPFOLATE_F))

# t-test needs the data to be normally-distributed, so I need to take the 
#    transformed values here (not 'orig')
t.test(SPFOLATE_F ~ stratum,
       data = mthfr_status_folate) %>%
  broom::tidy()

aov(SPFOLATE_F ~ stratum,
    data = mthfr_status_folate) %>%
  broom::tidy()

# ANCOVA - checking whether including 'diff_folate_intake' has impact on the correlation!
aov(SPFOLATE_F ~ stratum + as.factor(diff_folate_intake) + PC1 + PC2 + PC3,
    data = mthfr_status_folate) %>%
  broom::tidy()


# check also for the original categories, minor allele load
aov(SPFOLATE_F ~ as.factor(minor_al_load), data = mthfr_status_folate) %>%
  broom::tidy()

# check dose of folate intake in the two strata
mthfr_strat_plot_folate_diet <- ggplot(
  mthfr_status_folate,
  aes(x = stratum, y = FOLAT)
  ) + 
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  # coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_minimal() +
  xlab("rs1801133 stratum") +
  ylab("diet folate intake (ug/day)")
mthfr_strat_plot_folate_diet

t.test(FOLAT ~stratum, data = mthfr_status_folate)

# check also for the original categories, minor allele load
mthfr_3strat_plot_folate_diet <- ggplot(
  mthfr_status_folate,
  aes(x = as.factor(minor_al_load), y = FOLAT)
) + 
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  # coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_minimal() +
  xlab("rs1801133 T-allele load") +
  ylab("diet suppl.dose (ug/day)")
mthfr_3strat_plot_folate_diet

aov(FOLAT ~ as.factor(minor_al_load), data = mthfr_status_folate) %>%
  broom::tidy()


mthfr_strat_plot_folate_suppl <- ggplot(
  mthfr_status_folate,
  aes(x = stratum, y = s_folat)
) + 
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  # coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_minimal() +
  xlab("rs1801133 stratum") +
  ylab("folate suppl.dose (mg/day)")
mthfr_strat_plot_folate_suppl

t.test(s_folat ~stratum, data = mthfr_status_folate)

# check also for the original categories, minor allele load
mthfr_3strat_plot_folate_suppl <- ggplot(
  mthfr_status_folate,
  aes(x = as.factor(minor_al_load), y = s_folat)
) + 
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  # coord_cartesian(xlim = c(1.2, NA), clip = "off") +
  theme_minimal() +
  xlab("rs1801133 T-allele load") +
  ylab("folate suppl.dose (mg/day)")
mthfr_3strat_plot_folate_suppl

aov(s_folat ~ as.factor(minor_al_load), data = mthfr_status_folate) %>%
  broom::tidy()

# SAVE DATA FOR LATER ----
write_delim(
  best_prs_fit_manual,
  here("RESULTS", "PRS_manual_covs_folate_conc.txt"),
  delim = "\t"
)

write_delim(
  best_prs_fit,
  here("RESULTS", "PRS_auto_covs_folate_conc.txt"),
  delim = "\t"
)
