library(targets)

all_sources <- list.files("R_scripts", pattern = ".R$")
# read functions in all the scripts
purrr::walk(
  all_sources,
  ~ source(file.path("R_scripts", .x))
)

tar_option_set(packages = c("here", "patchwork", "dplyr", "ggplot2",
                            "readr", "tidyr", "data.table", "magrittr", "waldo"))

list(
  # FILES -----
  # base data
  tar_target(
    gwas_summ_stats_file, # GWAS summary statistics - original
    here("DATA", "Sfol_covar.assoc.add.MAF..txt"),
    format = "file"
  ),
  tar_target(
    info_all_snps_file, # illumina manifest file
    here("DATA", "HumanOmni1-Quad_v1-0_H.csv"),
    format = "file"
  ),
  # MoBa data
  tar_target(
    cntrl_ids_folate_cov_file,
    here("DATA", "folate_cov_merged.txt"),
    format = "file"
  ),
  # bash scripts
  tar_target(
    QC_plink_script,
    "QC_moba_data_plink.sh",
    format = "file"
  ),
  tar_target(
    QC_plink_snplist_script,
    "QC_moba_data_plink_snplist.sh",
    format = "file"
  ),
  tar_target(
    QC_plink_prune_script,
    "QC_moba_data_plink_snp-prune.sh",
    format = "file"
  ),
  tar_target(
    QC_plink_het_script,
    "QC_moba_data_plink_het.sh",
    format = "file"
  ),
  tar_target(
    QC_plink_final_script,
    "QC_moba_data_final.sh",
    format = "file"
  ),
  tar_target(
    QC_plink_pca_script,
    "QC_moba_data_calc_PCA.sh",
    format = "file"
  ),
  tar_target(
    SNP_list_4test_file, # list of SNPs that we want to test in PRS
    here("DATA", "SNP_list_4test.txt"),
    format = "file"
  ),
  # FUNCTIONS ----
  # base data
  tar_target(
    gwas_summ_stats,
    read_gwas(gwas_summ_stats_file)
  ),
  tar_target(
    info_all_snps,
    read_manifest(info_all_snps_file)
  ),
  tar_target(
    gwas_summ_cleaned,
    clean_gwas(gwas_summ_stats, info_all_snps)
  ),
  tar_target(
    gwas_summ_cleaned_file,
    write_and_return_gwas_cleaned(gwas_summ_cleaned),
    format = "file"
  ),
  tar_target(
    SNP_list_4test,
    scan(SNP_list_4test_file, what = "character")
  ),
  # MoBa data
  tar_target(
    cntrl_ids_folate_cov,
    read_and_transform_moba_cov(cntrl_ids_folate_cov_file)
  ),
  tar_target(
    mothers_ids_4_plink_file,
    write_and_return_moba_IDs(cntrl_ids_folate_cov),
    format = "file"
  ),
  tar_target(
    mothers_ids_folate_meas_blood_file, # "DATA/mothers_IDs_folate_meas_blood.txt"
    write_and_return_moba_IDs_spfolate(cntrl_ids_folate_cov),
    format = "file"
  ),
  # run PLINK script to combine all chromosomes - output: "22_combined.bim"
  tar_target(
    output_plink_combine,
    run_plink_combine(QC_plink_script, mothers_ids_4_plink_file),
    format = "file"
  ),
  # create .snplist file with all the SNPs that passed QC
  tar_target(
    output_plink_snplist,
    run_plink_snplist(QC_plink_snplist_script, output_plink_combine),
    format = "file"
  ),
  # Create list .prune.in with the pruned SNPs
  tar_target(
    output_plink_prune,
    run_plink_prune(QC_plink_prune_script, output_plink_snplist),
    format = "file"
  ),
  # Calculate heterozygosity
  tar_target(
    output_plink_het,
    run_plink_het(QC_plink_het_script, output_plink_prune),
    format = "file"
  ),
  # Check heterozygosity
  tar_target(
    het_calculated,
    read_table(output_plink_het)
  ),
  tar_target(
    output_filter_het, #cntrl_all_chr_QC_sample_OK.txt
    filter_het(het_calculated),
    format = "file"
  ),
  # Check SNP mismatches
  tar_target(
    merged_snps,
    create_merged_data(
      moba_file = output_plink_combine,
      gwas_summ_file = gwas_summ_cleaned_file,
      snps_moba_qc_file = output_plink_snplist
    )
  ),
  tar_target(
    # update the positions of SNPs in gwas_summary
    gwas_summ_cleaned_upd,
    update_snppos_gwas_summ(merged_snps),
    format = "file"
  ),
  tar_target(
    moba_bim_updated,
    check_complementary_recode(merged_snps, output_plink_combine)
  ),
  tar_target(
    moba_bim_upd_file,
    save_and_return_moba_bim(moba_bim_updated),
    format = "file"
  ),
  # FINAL QC STEPS
  tar_target(
    moba_data_qced_file, #here("DATA", "cntrl_all_chr_QC.bed")
    run_plink_final_qc(
      QC_plink_final_script,
      output_filter_het,
      moba_bim_upd_file
    ),
    format = "file"
  ),
  # Calculate PCA on genotype data
  tar_target(
    gen_pca_file, # here("DATA", "cntrl_all_chr_QC_PCA.eigenvec")
    run_plink_pca(QC_plink_pca_script, moba_data_qced_file),
    format = "file"
  ),
  tar_target(
    gen_pca,
    read_table(gen_pca_file)
  ),
  # combine genotype PCA with the covariate based on MoBa questions
  tar_target(
    prs_covariate,
    combine_covs(gen_pca, npc = 3, moba_data = cntrl_ids_folate_cov)
  ),
  tar_target(
    prs_cov_file,
    save_and_return_covs(prs_covariate),
    format = "file"
  ),
  # run PRSice2
  tar_target(
    best_prs_fit,
    run_prsice2(
      base_data = gwas_summ_cleaned_upd, #"DATA/base_data_cleaned.txt"
      target_data = stringr::word(
        string = moba_data_qced_file,
        start = 1,
        end = -2,
        sep = stringr::fixed(".")
      ), # "DATA/cntrl_all_chr_QC"
      pheno_data = mothers_ids_folate_meas_blood_file, # "DATA/mothers_IDs_folate_meas_blood.txt"
      all_cov_data = prs_cov_file, #"DATA/covariates_4prs.txt"
      out_prefix = "RESULTS/MoBa_folat_prsice_with-cov",
      nq = 5
    )
  ),
  tar_target(
    best_prs_fit_manual,
    run_prsice2_manual(
      base_data = gwas_summ_cleaned_upd, #"DATA/base_data_cleaned.txt"
      target_data = stringr::word(
        string = moba_data_qced_file,
        start = 1,
        end = -2,
        sep = stringr::fixed(".")
      ), # "DATA/cntrl_all_chr_QC"
      pheno_data = mothers_ids_folate_meas_blood_file, # "DATA/mothers_IDs_folate_meas_blood.txt"
      all_cov_data = prs_cov_file, #"DATA/covariates_4prs.txt"
      prs_levels = 0.0001500, # this is manually extracted after checking the automatically calculated PRS in 'Run_PRSice2.Rmd'
      out_prefix = "RESULTS/MoBa_folat_prsice_with-cov_selected"
    )
  ),
  # calculating PRS using PLINK and a chosen set of SNPs
  tar_target(
    extracted_gwas_summ_stats_file,
    extract_defined_SNPs(
      base_data = gwas_summ_cleaned_upd, #"DATA/base_data_cleaned.txt"
      SNP_list = SNP_list_4test
    ),
    format = "file"
  ),
  tar_target(
    extracted_gwas_summ_stats, # only the SNPs that we want to test in PRS
    read_table(extracted_gwas_summ_stats_file)
  ),
  tar_target(
    plink_prs_manual,
    run_plink_manual(
      base_data = extracted_gwas_summ_stats_file, #DATA/extracted_gwas_summ_stats.txt
      target_data = stringr::word(
        string = moba_data_qced_file,
        start = 1,
        end = -2,
        sep = stringr::fixed(".")
      ), # "DATA/cntrl_all_chr_QC"
      pheno_data = mothers_ids_folate_meas_blood_file, # "DATA/mothers_IDs_folate_meas_blood.txt"
      # all_cov_data = prs_cov_file, #"DATA/covariates_4prs.txt"
      out_prefix = "RESULTS/MoBa_folat_plink_selected"
    )
  )
)
