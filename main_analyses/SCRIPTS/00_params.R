# Constants and parameters

# types of analyses
all_res_types <- c(
  "res_prelim_aed", "res_prelim_no_aed", #stratification by use of AED
  "res_prelim_interact", # no stratification, interaction term
  #additional strat. by MTHFR-SNP genotype
  "res_prelim_mthfr_no_aed", "res_prelim_mthfr_aed",
  #additional strat. use of folate suppl.
  "res_prelim_fol_int_no_aed", "res_prelim_fol_int_aed"
)

# all the genotype-scores:
# 'prs_auto' - automatically chosen PRS (best fit, 2 SNPs),
# 'prs_manual' - manually chosen 2nd best PRS (76 SNPs),
# 'MTHFR_2strat' - categorizing by the genotype on MTHFR-SNP
all_prs_types <- c("prs_auto", "prs_manual", "MTHFR_2strat")

# covariates
explanatory_vars <- c("FOL_PRE_or_TRIM1", "FOLAT", "s_folat")
explanatory_names <- c("FOL_PRE_or_TRIM1", "FOLAT", "s_folat")
names(explanatory_vars) <- explanatory_names
