# calling final PLINK scripts

run_plink_final_qc <- function(
  script_in,
  filtered_het_samples,
  moba_bim_updated_file
){
  if(!file.exists(filtered_het_samples)){
    stop("File ", filtered_het_samples, " doesn't exist!", call. = FALSE)
  }
  # create final .bim, .bed, .fam files
  if(!file.exists(moba_bim_updated_file)){
    stop("File ", moba_bim_updated_file, " doesn't exist!", call. = FALSE)
  }
  
  system(
    paste0("mv ", here("DATA", "22_combined.bim"), " ",
           here("DATA", "22_combined.bim.backup"))
  )
  system(
    paste0("ln -s ", moba_bim_updated_file, " ", here("DATA", "22_combined.bim"))
  )
  
  system(paste0("./", script_in, " >& ", script_in, ".out"))
  return(here("DATA", "cntrl_all_chr_QC.bed"))
}

run_plink_pca <- function(script_in, qced_gen_file){
  if(!file.exists(qced_gen_file)){
    stop("File ", qced_gen_file, " doesn't exist!")
  }
  system(paste0("./", script_in, " >& ", script_in, ".out"))
  return(here("DATA", "cntrl_all_chr_QC_PCA.eigenvec"))
}
