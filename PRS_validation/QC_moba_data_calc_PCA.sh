#!/bin/bash

cd DATA

# First, we need to perform prunning
../plink --bfile cntrl_all_chr_QC --indep-pairwise 200 50 0.25 --out cntrl_all_chr_QC_forPCA

# calculate PCA based on genotypes
../plink --bfile cntrl_all_chr_QC --extract cntrl_all_chr_QC_forPCA.prune.in --pca 10 'header' --out cntrl_all_chr_QC_PCA

cd ..
