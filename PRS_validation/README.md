This directory contains files that were used to calculate and check PRS for
folate levels in blood among women _without_ epilepsy.

The process is based on the tutorial:
 https://choishingwan.github.io/PRS-Tutorial/

## FOLDERS

- `R_scripts` - all the R scripts
- `_targets` - automatically generated objects and states (`targets` package)
- `DATA` - all the raw data and created data and results _(not available here)_
- `FIGURES` - all the figures created on the go
- `RESULTS` - results from PRSice2 run _(not available here)_

## ANALYSIS PROCESS

> **!!!NOTE ON targets PACKAGE USE!!!**    
> The entire process is in `_targets.R` script.
> The updates, checks, etc. can be run through the commands in `_targets.Rmd`.
> Read documentation of the targets package before modifying.

### R_scripts/QC_base_data.R

Quality Control of the GWAS summary data

**INPUT:**    
  `Sfol_covar.assoc.add.MAF..txt` - gwas summary stats

**OUTPUT:**    
  `base_data_cleaned.txt`

### R_scripts/QC_moba_data.R

Intro to Quality control of the MoBa folate measurements - extracting the genetic data IDs

**INPUT:**    
  `mothers_folate_suppl_genetic.txt` - IDs of mothers that have genetic data and
  folate measurements in blood

**OUTPUT:**    
`mothers_IDs_4plink.txt` - IDs from genetic data

### QC_moba_data_plink.sh

Quality control of the genotype data, using plink; part one: combining all
chromosomes, filtering only the samples in the input file.

**INPUT:**    
  genotypes, separately per each chromosome    
  `DATA/mothers_IDs_4plink.txt`

**OUTPUT:**
  `22_combined.*` - merged genotypes    
  `cntrl_all_chr_QC.*` - QCed data

### QC_moba_data_plink_snplist.sh

Quality control of the genotype data, using plink; part two: generating
the list of SNPs

### QC_moba_data_plink_snp-prune.sh

Quality control of the genotype data, using plink; part three: pruning

### QC_moba_data_plink_het.sh

Quality control of the genotype data, using plink; part four: calculating
heterozygosity
    
### R_scripts/QC_moba_data_hetero.R

Checking heterozygosity calculated using plink in the previous step

**INPUT:**    
  `cntrl_all_chr_QC.het` - information about hetegozygosity of each sample that survived QC

**OUTPUT:**    
  `cntrl_all_chr_QC_sample_OK.txt` - IDs of samples that have ok heterozygosity

### QC_moba_data_snp-mismatch.R

Checking which SNPs need to be strand-flipped and which are completely
not matching between the summary stat. data and MoBa Gen.

**INPUT:**    
  `22_combined.bim` - information on markers in MoBa chr1-22    
  `cntrl_all_chr_QC.snplist` - list of markers that survived QC    
  `base_data_cleaned.txt` - cleaned GWAS summary stat

**OUTPUT:**    
  `22_combined_recoded.bim` - information on which allele should be A1 in each of the SNPs    
  `22_combined_mismatch.txt` - list of the SNPs that did not have match in GWAS summary stat

**NOTE:** in the tutorial, one checks correctness of samples by checking
the given sex vs. the genotype sex, but I don't have sex chromosomes here
and this check was done in QC for all MoBa

### QC_moba_data_final.sh

Checking relatedness and generating the complete data after QC
output of commands: QC_moba_data_final.out

**INPUT:**    
  `22_combined_recoded(.bim, .bam, .fam)` - genotypes from MoBa    
  SNPs filtered by: `cntrl_all_chr_QC.prune.in`    
  samples filtered by: `cntrl_all_chr_QC_sample_OK.txt`

**OUTPUT:**    
  `cntrl_all_chr_QC.rel.id` - people that are not excluded by relatedness check

### QC_moba_data_calc_PCA.sh

Calculating PCA of genotypes. This will be used as covariates in PRS.

**OUTPUT:**    
  `cntrl_all_chr_QC_PCA.*` - eigenvectors and eigenvalues for first 10 PCA


### OTHER FILES:

- `PRSice.R` - original file from PRSice2, used to run the analysis

- `Run_PRSice2.Rmd` - creating the covariates and running analysis.

- `MTHFR_SNP_recode_controls.R` - recode the rs1801133 SNP in controls to check 
   minor allele (T) loading

- `PRS_folate_conc_check.R` - showing the relation btwn calculated PRS and folate
   concentration in blood in controls; also folate conc. depending on the
   rs1801133 SNP genotype

- `PRS_auto_MAF_table.R:`

    - checking the frequency counts and Hardy-Weinberg Eq.tests for each SNP in
      the PRS chosen by PRSice2
    - this script uses PLINK to produce the files:    
    
        - in DATA folder:    
    `cntrl_prs_auto_snps.*` - files produced by PLINK when extracting SNPs from
    the combined, 'cntrl_all_chr_QC' files
    
        - in the RESULTS folder:    
    `.frqx` (genotype count report) - Produced by `--freqx`. Valid input for `--read-freq`.    
    .hwe (Hardy-Weinberg equilibrium exact test statistic report) - Produced by `--hardy`.

