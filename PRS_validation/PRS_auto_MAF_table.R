# DESCRIPTION: Check MAF for all SNPs in the PRS
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-12-12
# DATE MODIFIED: 

# SETUP ----
library(readr)
library(dplyr)
library(here)
library(flextable)

# READ DATA ----
prs_auto_snp_file <- here("RESULTS", "MoBa_folat_prsice_with-cov_selected.snp")
snps_prs_auto <- read_table(prs_auto_snp_file)
snps_prs_auto

moba_data_qced <- here("DATA", "cntrl_all_chr_QC")

# write the snps in a list - format needed for PLINK run
snp_list_4plink <- here("RESULTS", "snps_list_prs_auto_4PLINK.txt")
write_delim(
  snps_prs_auto %>%
    select(SNP),
  file = snp_list_4plink,
  delim = "\t",
  col_names = FALSE
)

# RUN PLINK ----
## 1. EXTRACTING SNPs ----
extracted_snps_base_file <- here("DATA", "cntrl_prs_auto_snps")
plink_command_extract <- paste0(
  "./plink --bfile ", moba_data_qced, " --allow-no-sex --extract ", snp_list_4plink,
  " --make-bed --out ", extracted_snps_base_file
)
plink_command_extract
system(plink_command_extract, intern = TRUE)

list.files(path = here("DATA"), pattern = "cntrl_prs_auto_snps")

## 2. COUNTING FREQUENCY AND HWE TEST ----
count_freq_file <- here("RESULTS", "cntrl_prs_auto_snps_stats")
plink_command_stats <- paste0(
  "./plink --bfile ", extracted_snps_base_file,
  " --allow-no-sex --nonfounders --freq --hardy --out ", count_freq_file
)
plink_command_stats
system(plink_command_stats, intern = TRUE)

list.files(path = here("RESULTS"), pattern = "cntrl_prs_auto")

# CHECK RESULTS ----
prs_auto_hwe <- read_table(
  paste0(count_freq_file, ".hwe"),
  col_names = c("CHR", "SNP", "TEST", "A1", "A2", "O_HET", "E_HET", "P", "nothing"),
  skip = 1
) %>%
  select(-nothing)
prs_auto_hwe

prs_auto_freq <- read_table(
  paste0(count_freq_file, ".frq")
)
prs_auto_freq

# NICE TABLE ----
output_table <- prs_auto_hwe %>%
  select(-TEST, -E_HET) %>%
  left_join(prs_auto_freq)

small_border <- officer::fp_border(color = "gray", width = 1)

nice_table <- output_table %>%
  select(-CHR) %>%
  mutate(
    P = round(P, digits = 2),
    MAF = round(MAF, digits = 2)
  ) %>%
  select(SNP, A1, A2, everything()) %>%
  flextable() %>%
  set_header_labels(
    A1 = "minor allele",
    A2 = "reference allele",
    P = "HWE test p-value",
    N_OBS = "n.observations",
    O_HET = "heterozygote freq."
  ) %>%
  merge_v(
    j = c("SNP", "A1", "A2"),
    combine = TRUE
  ) %>%
  border_inner_h(
    border = small_border,
    part = "body"
  )
nice_table

## SAVE ----
save_as_docx(nice_table, path = here("RESULTS", "MAF_table_PRS_auto_validation.docx"))

write_delim(
  output_table,
  path = here("RESULTS", "MAF_HWE_PRS_auto_validation.txt"),
  delim = "\t"
)
