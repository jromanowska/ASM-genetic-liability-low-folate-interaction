# DESCRIPTION: Functions that prepare for running PRSice2, after all QC
# AUTHOR: Julia Romanowska
# DATE CREATED: 2021-11-30


#' @param gen_pca tibble with genetic PCs (first two columns are assumed to
#'   be family ID and personal ID, named 'FID' and 'IID', respectively)
#' @param npc how many PCs to take into account? (default: 6)
#' @param moba_data covariates from MoBa (with individual IDs named 'SENTRIXID')
#' 
#' @return tibble with both covariates together
combine_covs <- function(gen_pca, npc = 6, moba_data){
  max_npc <- ncol(gen_pca) - 2
  if(npc > max_npc){
    stop("The requested number of PCs to copy to the final covariate set is larger than the number of available PCs!")
  }
  
  # the first two columns of 'gen_pca' are assumed to be indiv.id. and fam.id
  cov_all <- full_join(
    gen_pca[, 1:(npc + 2)],
    moba_data %>%
      # select(-M_ID_1111) %>%
      select(SENTRIXID:diff_folate_intake),
    by = c("FID" = "SENTRIXID")
  )
  return(cov_all)
}

save_and_return_covs <- function(cov_all){
  out_file <- here("DATA", "covariates_4prs.txt")
  write_delim(
    cov_all,
    file = out_file,
    delim = "\t"
  )
  return(out_file)
}

#' @description Run the PRSice2 analysis and get the best fit, then categorize
#'   it into quantiles.
#' 
#' @param base_data filename of the file with GWAS summary statistics
#' @param target_data the core of the file name for .bed, .bim, and .fam files
#'   with genotypes for the target population
#' @param pheno_data filename of the file containing the phenotype data
#'   (columns: FID, IID, pheno)
#' @param all_cov_data filename of the file containing all covariates together
#'   (columns: FID, IID, cov1, cov2, ...)
#' @param out_prefix character for naming all the output files
#' @param nq number of quantiles to categorize the best fit PRS into
#'   (deafult: 10)
#' 
#' @return data with the best PRS fit: tibble with the following columns:
#'   \itemize{
#'     \item FID - family ID
#'     \item IID - individual ID
#'     \item In_Regression - whether the sample was included in the regression
#'     \item PRS - calculated best fit
#'     \item PRS_quant - the PRS variable categorized into nq quantiles
#'   }
run_prsice2 <- function(
  base_data,
  target_data,
  pheno_data,
  all_cov_data,
  out_prefix,
  nq = 10
){
  run_command <- sprintf(
    "Rscript PRSice.R --prsice PRSice_linux --base %s --target %s --binary-target F --pheno %s --base-maf MAF:0.01 --stat BETA --beta --cov %s --cov-factor 'diff_folate_intake' --out %s",
    base_data,
    target_data,
    pheno_data,
    all_cov_data,
    out_prefix
  )
  error_code <- system(
    run_command,
    intern = FALSE,
    ignore.stdout = TRUE
  )
  if(error_code != 0){
    stop("There was an error running PRSice2! Check the log!")
  }
  
  # get the best fit and create quantiles
  best <- read_delim(
    paste0(out_prefix, ".best"),
    delim = " ")
  best <- best %>%
    mutate(PRS_quant = cut(
        PRS,
        breaks = nq,
        include.lowest = TRUE,
        labels = 1:nq
      )
    )
  
  return(best)
}

#' @description Run the PRSice2 analysis and output all the scores, and SNPs.
#'   Does not perform regression to choose the optimal fit!
#' 
#' @param base_data filename of the file with GWAS summary statistics
#' @param target_data the core of the file name for .bed, .bim, and .fam files
#'   with genotypes for the target population
#' @param pheno_data filename of the file containing the phenotype data
#'   (columns: FID, IID, pheno)
#' @param all_cov_data filename of the file containing all covariates together
#'   (columns: FID, IID, cov1, cov2, ...)
#' @param prs_levels numeric vector with thresholds at which to calculate PRS
#' @param out_prefix character for naming all the output files
#' @param nq number of quantiles to categorize the best fit PRS into
#'   (deafult: 10)
#' 
#' @return data with the best PRS fit: tibble with the following columns:
#'   \itemize{
#'     \item FID - family ID
#'     \item IID - individual ID
#'     \item In_Regression - whether the sample was included in the regression
#'     \item PRS - calculated best fit
#'     \item PRS_quant - the PRS variable categorized into nq quantiles
#'   }
run_prsice2_manual <- function(
  base_data,
  target_data,
  pheno_data,
  all_cov_data,
  prs_levels,
  out_prefix,
  nq = 10
){
  run_command <- sprintf(
    "Rscript PRSice.R --prsice PRSice_linux --base %s --target %s --binary-target F --pheno %s --base-maf MAF:0.01 --stat BETA --beta --cov %s --cov-factor 'diff_folate_intake' --bar-levels %s --fastscore --no-full --print-snp --out %s",
    base_data,
    target_data,
    pheno_data,
    all_cov_data,
    paste0(prs_levels, collapse = ","),
    out_prefix
  )
  error_code <- system(
    run_command,
    intern = FALSE,
    ignore.stdout = TRUE
  )
  if(error_code != 0){
    stop("There was an error running PRSice2! Check the log!")
  }
  
  # get the best fit and create quantiles
  best <- read_delim(
    paste0(out_prefix, ".best"),
    delim = " ")
  best <- best %>%
    mutate(PRS_quant = cut(
      PRS,
      breaks = nq,
      include.lowest = TRUE,
      labels = 1:nq
    )
    )
  
  return(best)
}
