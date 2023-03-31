# DESCRIPTION: checking base data, i.e., GWAS summary data
# AUTHOR: Julia Romanowska
# DATE CREATED: 2021-11-02
# DATE MODIFIED: 2021-11-29


# Read data -----------
read_gwas <- function(file_in){
  read_delim(
    file = file_in,
    delim = " "
  )
}

# manifest file
read_manifest <- function(file_in){
  info_all_snps <- read_csv(
    here("DATA", "HumanOmni1-Quad_v1-0_H.csv"),
    skip = 7
    
  )
  # this file is not really one csv, so there are some problematic rows
  #   - these can be removed
  info_all_snps <- info_all_snps %>%
    slice_head(n = nrow(info_all_snps) - nrow(problems(info_all_snps)))
#  info_all_snps
  
  # just checking whether the 'Name' column includes SNP names with 'rs'
  # info_all_snps %>%
  #   transmute(rs_name = (stringr::str_sub(Name, 1, 2) == "rs")) %>%
  #   count(rs_name)
  return(info_all_snps)
}

# Tidy data -----------
clean_gwas <- function(gwas_data, snps_data){
  # # need to check duplicates and rename the columns
  # identical(gwas_summary$CHR...10, gwas_summary$CHR...2)
  # identical(gwas_summary$A1...4, gwas_summary$A1...11)
  
  gwas_summary <- gwas_data %>%
    select(-CHR...10, -A1...4) %>%
    rename(CHR = CHR...2, A1 = A1...11)

  # QC -------
  # filter out SNPs with MAF <= 0.01
  gwas_summary_okMAF <- gwas_summary %>%
    filter(MAF > 0.01)

  # check for duplicated SNPs
  dupl_snps <- gwas_summary_okMAF %>%
    count(SNP) %>%
    filter(n > 1)
  if(length(dupl_snps) > 0){
    gwas_summary_okMAF <- gwas_summary_okMAF %>%
      group_by(SNP) %>%
      slice(1) %>%
      ungroup()
  }

  # remove ambiguous SNPs, i.e., those with complementary alleles
  gwas_summary_okMAF <- gwas_summary_okMAF %>% 
    filter(!((A1 == 'C' & A2 == 'G') |
             (A1 == 'G' & A2 == 'C') |
             (A1 == 'A' & A2 == 'T') |
             (A1 == 'T' & A2 == 'A')))

  # # check p-values
  # ggplot(gwas_summary_okMAF) +
  #   geom_qq(aes(sample = P))
  
  # recode SNP names
  gwas_summary_okMAF <- gwas_summary_okMAF %>%
    mutate(SNP = tolower(SNP))
  
  # amend genomic positions
  #    the positions from the summary stats. were in GRCh36,
  #    while in the manifest file, the positions were in GRCh37;
  #    this is also the version in MoBa v1.0 merged
  gwas_summary_okMAF <- gwas_summary_okMAF %>%
    left_join(
      snps_data %>%
        select(Name, alleles = SNP, GenomeBuild, Chr, MapInfo) %>%
        filter(Chr %in% as.character(1:22)) %>%
        mutate(Chr = as.numeric(Chr)),
      by = c("SNP" = "Name", "CHR" = "Chr")
    )

  return(gwas_summary_okMAF)
}

# Save data ----
write_and_return_gwas_cleaned <- function(data_in){
  write_delim(
    data_in,
    file = here("DATA", "base_data_cleaned.txt"),
    delim = "\t"
  )
  return(here("DATA", "base_data_cleaned.txt"))
}
