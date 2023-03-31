# DESCRIPTION: checking folate measurements
# AUTHOR: Julia Romanowska
# DATE CREATED: 2021-11-02
# DATE MODIFIED: 2022-02-07


# Read data -----------
# Since the folate is not normally-distributed, we need to re-normalize it
read_and_transform_moba_cov <- function(file_in){
  data_out <- read_delim(
    file = file_in,
    delim = "\t"
  )
  # folate measurements are called SPFOLATE_F
  #  - create a new variable to hold the original values while the
  #    transformed values will be overwritten in the old variable
  # source("transformation_functions.R")
  data_out <- data_out %>%
    dplyr::mutate(
      orig_SPFOLATE_F = SPFOLATE_F
    ) %>%
    dplyr::mutate(
      SPFOLATE_F = inttransform(orig_SPFOLATE_F, data = .)
    )
  return(data_out)
}

# Extract data ----
# I will need another format for PLINK and PRSice2
write_and_return_moba_IDs <- function(data_in){
  # Write sample IDs to use when filtering in plink
  write_delim(
    data_in %>%
      dplyr::select(fam_id = SENTRIXID, indiv_id = SENTRIXID),
    file = here("DATA", "mothers_IDs_4plink.txt"),
    delim = "\t",
    col_names = FALSE
  )
  return(here("DATA", "mothers_IDs_4plink.txt"))
}

# 
write_and_return_moba_IDs_spfolate <- function(data_in){
  # Write sample IDs and folate measuremetns
  write_delim(
    data_in %>%
      dplyr::select(fam_id = SENTRIXID, indiv_id = SENTRIXID, SPFOLATE_F),
    file = here("DATA", "mothers_IDs_folate_meas_blood.txt"),
    delim = "\t"
  )
  return(here("DATA", "mothers_IDs_folate_meas_blood.txt"))
}

# Run the PLINK scripts ----
run_plink_combine <- function(script_in, ids_file){
  if(!file.exists(ids_file)){
    stop("File ", ids_file, " doesn't exist!")
  }
  system(paste0("./", script_in, " >& ", script_in, ".out"))
  return(here("DATA", "22_combined.bim"))
}

run_plink_snplist <- function(script_in, combined_file){
  if(!file.exists(combined_file)){
    stop("File ", combined_file, " doesn't exist!")
  }
  system(paste0("./", script_in, " >& ", script_in, ".out"))
  return(here("DATA", "cntrl_all_chr_QC.snplist"))
}

run_plink_prune <- function(script_in, snplist_file){
  if(!file.exists(snplist_file)){
    stop("File ", snplist_file, " doesn't exist!")
  }
  system(paste0("./", script_in, " >& ", script_in, ".out"))
  return(here("DATA", "cntrl_all_chr_QC.prune.in"))
}

run_plink_het <- function(script_in, prune_file){
  if(!file.exists(prune_file)){
    stop("File ", prune_file, " doesn't exist!")
  }
  system(paste0("./", script_in, " >& ", script_in, ".out"))
  return(here("DATA", "cntrl_all_chr_QC.het"))
}
