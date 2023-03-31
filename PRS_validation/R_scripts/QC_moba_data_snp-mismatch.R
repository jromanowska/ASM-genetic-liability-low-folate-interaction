# DESCRIPTION: Check SNP mismatch between the GWAS summary stat.data and
#    MoBa Genetics: there might be some SNPs where the strands are flipped.
# AUTHOR: Julia Romanowska
# DATE CREATED: 2021-11-09
# DATE MODIFIED: 2021-11-29


create_merged_data <- function(moba_file, gwas_summ_file, snps_moba_qc_file){
  # READ DATA -----
  # Read in bim file  - all SNPs in MoBa Genetics data
  MoBa_bim <- fread(
      moba_file, #here("DATA", "22_combined.bim")
      col.names = c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
    ) %>%
    .[,c("B.A1","B.A2") := list(toupper(B.A1), toupper(B.A2))]

  # Read in summary stats
  gwas_summ <- fread(
      gwas_summ_file #here("DATA", "base_data_cleaned.txt")
    ) %>%
    .[,c("A1","A2") := list(toupper(A1), toupper(A2))]

  # Read in QCed SNPs
  snps_moba_qc <- fread(
    snps_moba_qc_file, #here("DATA", "cntrl_all_chr_QC.snplist"),
    header = FALSE,
    col.names = "SNP"
  )

  # WHICH SNPS REQUIRE STRAND FLIPPING? -----
  # Merge summary statistic with target
  merged_data <- merge(
      MoBa_bim,
      gwas_summ,
      by = c("SNP", "CHR")
    ) %>%
    # And filter out QCed SNPs
    .[SNP %in% snps_moba_qc$SNP]

  # CHECK THE GENOMIC POSITIONS
  merged_data_noNAs <- merged_data[!is.na(MapInfo)]
#  waldo::compare(merged_data_noNAs$MapInfo, merged_data_noNAs$BP.x)
  
  # very few positions that are different - I'll write it down and check 
  # afterwards if these SNps are in the PRS
  diff_merged_gen_positions <- waldo::compare(merged_data_noNAs$MapInfo, merged_data_noNAs$BP.x)
  writeLines(
    diff_merged_gen_positions,
    con = here("DATA", "diff_merged_gen_positions.txt")
  )
  
  return(merged_data)
}

update_snppos_gwas_summ <- function(merged_data){
  # thus, I can use the positions from MoBa, in GRCh37
  new_gwas_summ <- merged_data[,
    .(CHR, BP = BP.x, SNP, A1, A2, TEST, NMISS, BETA, STAT, P, MAF)
  ]
  
  out_file <- here("DATA", "base_data_cleaned.txt")
  readr::write_delim(
    new_gwas_summ,
    file = out_file,
    delim = "\t"
  )
  
  return(out_file)
}

# Function for calculating the complementary allele
complement <- function(x){
  switch (x,
          "A" = "T",
          "C" = "G",
          "T" = "A",
          "G" = "C",
          return(NA)
  )
} 

check_complementary_recode <- function(merged_data, moba_bim_file){
  # Read in bim file  - all SNPs in MoBa Genetics data
  moba_bim <- fread(
    moba_bim_file, #here("DATA", "22_combined.bim")
    col.names = c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
  ) %>%
    .[,c("B.A1","B.A2") := list(toupper(B.A1), toupper(B.A2))]
  
  # Get SNPs that have the same alleles across base and target
  # A1, A2 is the base data
  # B.A1, B.A2 is the MoBa data
  merged_snps_match <- merged_data[A1 == B.A1 & A2 == B.A2, SNP]
  # Identify SNPs that are complementary between base and target
  merged_snps_complementary <- merged_data[
    sapply(B.A1, complement) == A1 &
      sapply(B.A2, complement) == A2,
    SNP
  ]
  # Now update the bim file
  moba_bim[
    SNP %in% merged_snps_complementary,
    c("B.A1", "B.A2") :=
        list(sapply(B.A1, complement),
             sapply(B.A2, complement))
  ]
  
  
  # identify SNPs that need recoding -----
  merged_snps_recode <- merged_data[B.A1==A2 & B.A2==A1, SNP]
  # Update the bim file
  moba_bim[
    SNP %in% merged_snps_recode,
    c("B.A1", "B.A2") :=
        list(B.A2, B.A1)
  ]
  
  # identify SNPs that need recoding & complement -----
  merged_snps_recorde_compl <- merged_data[
    sapply(B.A1, complement) == A2 &
    sapply(B.A2, complement) == A1,
    SNP
  ]
  # Now update the bim file
  moba_bim[
    SNP %in% merged_snps_recorde_compl,
    c("B.A1", "B.A2") :=
        list(sapply(B.A2, complement),
             sapply(B.A1, complement))
  ]

  # save also info on which SNPs do not match at all
  #  this is usually due to difference in genome build
  snps_mismatch <- moba_bim[
    !(SNP %in% merged_snps_match |
        SNP %in% merged_snps_recorde_compl |
        SNP %in% merged_snps_complementary |
        SNP %in% merged_snps_recode),
    SNP
  ]
  length(snps_mismatch)
  write.table(
    snps_mismatch,
    here("DATA", "22_combined_mismatch.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  return(moba_bim)
}

# Write the updated bim file to disk ----
save_and_return_moba_bim <- function(moba_bim){
  out_file <- here("DATA", "22_combined_recoded.bim")
  # write the updated .bim file
  fwrite(
    x = moba_bim,
    file = out_file,
    sep = "\t",
    col.names = FALSE
  )
  
  return(out_file)
}

