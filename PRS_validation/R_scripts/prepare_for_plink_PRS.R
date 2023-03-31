#' @description Extract the SNPs from the GWAS summary stats
#' 
#' @param base_data filename of the file with GWAS summary statistics
#' @param SNP_list characte vector with the SNP names
#' 
#' @return part of the GWAS summary stat
extract_defined_SNPs <- function(
  base_data,
  SNP_list
){
  data_in <- read_table(base_data)
  non_matched_SNPs <- SNP_list[which(
    !(SNP_list %in% data_in$SNP)
  )]
  if(length(non_matched_SNPs) != 0){
    warning("These SNPs were not found: ", paste(non_matched_SNPs, collapse = ","))
  }
  
  data_out <- data_in %>%
    filter(SNP %in% SNP_list)
  if(nrow(data_out) == 0){
    stop("No SNPs found!", call. = FALSE)
  }
  
  out_file <- here("DATA", "extracted_gwas_summ_stats.txt")
  write_delim(data_out, file = out_file, delim = "\t")
  return(out_file)
}

#' @description Run PLINK to calculate PRS that includes only the chosen SNPs,
#'   based on literature search.
#' 
#' @param base_data filename of the file with extracted SNPs from the GWAS
#'   summary statistics
#' @param target_data the core of the file name for .bed, .bim, and .fam files
#'   with genotypes for the target population
#' @param pheno_data filename of the file containing the phenotype data
#'   (columns: FID, IID, pheno)
#' @param out_prefix character for naming all the output files
#' @param nq number of quantiles to categorize the best fit PRS into
#'   (deafult: 10)
#' 
#' @return data with the best PRS fit: tibble with the following columns:
#'   \itemize{
#'     \item FID - family ID
#'     \item IID - individual ID
#'     \item CNT, CNT2 - output from PLINK
#'     \item SCORE - calculated PRS
#'     \item SPFOLATE_F - measured folate concentration
#'     \item PRS_quant - the PRS variable categorized into nq quantiles
#'   }
run_plink_manual <- function(
  base_data,
  target_data,
  pheno_data,
  out_prefix,
  nq = 10
){
  run_command <- sprintf(
    "./plink --bfile %s --score %s 3 4 8 header --out %s",
    target_data,
    base_data,
    out_prefix
  )
  error_code <- system(
    run_command,
    intern = FALSE,
    ignore.stdout = TRUE
  )
  if(error_code != 0){
    stop("There was an error running PLINK! Check the log!")
  }
  
  # get the best fit and create quantiles
  pheno <- read_table(pheno_data)
  best <- read_table(paste0(out_prefix, ".profile")) %>%
    select(-PHENO) %>%
    left_join(
      pheno,
      by = c("FID" = "fam_id", "IID" = "indiv_id")
    ) %>%
    mutate(
      PRS_quant = cut(
        SCORE,
        breaks = nq,
        include.lowest = TRUE,
        labels = 1:nq
      )
    )
  
  return(best)
}
