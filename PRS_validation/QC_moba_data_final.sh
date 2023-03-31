#!/bin/bash
# this file will collect all the logs:
all_logs_file='PLINK_all_logs.log'

cd DATA

# calculate relatedness
../plink --bfile 22_combined --extract cntrl_all_chr_QC.prune.in --keep cntrl_all_chr_QC_sample_OK.txt --rel-cutoff 0.125 --out cntrl_all_chr_QC

cat cntrl_all_chr_QC.log >> ${all_logs_file}

../plink --bfile 22_combined --make-bed --keep cntrl_all_chr_QC.rel.id --out cntrl_all_chr_QC --extract cntrl_all_chr_QC.snplist --exclude 22_combined_mismatch.txt

cat cntrl_all_chr_QC.log >> ${all_logs_file}

cd ..
