# DESCRIPTION: recode MTHFR allele with PLINK
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-12-01
# DATE MODIFIED:

# SETUP ----
## LIBRARIES ----
library(readr)
library(here)
library(dplyr)

## PARAMETERS ----
target_data_prefix <- here("DATA", "cntrl_all_chr_QC")
mthfr_snp_file <- here("DATA", "MTHFR_SNP.txt")
mthfr_status_file <- here("RESULTS", "mthfr_status")


# RECODE MTHFR allele ----
# check whether rs1801133 SNP is CC, CT, or TT
write_lines(x = "rs1801133", file = mthfr_snp_file)

# NOTE: '--recode 12' option means that the minor allele (T) will be coded
#       as 1, while major allele (C), 2
plink_command <- sprintf(
  "./plink --bfile %s --extract %s --recode 12 --allow-no-sex --out %s",
  target_data_prefix,
  mthfr_snp_file,
  mthfr_status_file
)
error_code <- system(
  plink_command,
  intern = FALSE,
  ignore.stdout = TRUE
)
if(error_code != 0){
  stop("There was an error running PLINK! Check the log!")
}

