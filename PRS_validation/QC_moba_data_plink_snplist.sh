#!/bin/bash

# DESCRIPTION: Extracting .snplist
#    Part 2 of simple QC of MoBa Genotype data; based on
#    https://choishingwan.github.io/PRS-Tutorial/target/
# AUTHOR: Julia Romanowska

# main directory with the genetic files:
moba_plink_dir=''

# I will need to do that in the no-backup directory because these are large files
working_dir=$(pwd)
cd 'MOBA_plink_cntrl_group'

# this is the base of the names for the files that combine all chromosomes
out_fname_base='cntrl_all_chr_QC'

# this file will collect all the logs:
all_logs_file='PLINK_all_logs.log'

# now - the samples are extracted and the SNPs are QCed, but the samples are not QCed
./plink --bfile 22_combined --mind 0.01 --write-snplist --make-just-fam --out $out_fname_base

if [ $? -ne 0 ]
then
  cd $working_dir
  exit 1
fi

# synchronise newest data (one should replace '~/DATA' with a correct path)
rsync -auvz 22_combined.* ~/DATA/
rsync -auvz ${out_fname_base}.* ~/DATA/

# get back to R-project directory
cd $working_dir

cd DATA
cat ${out_fname_base}.log >> ${all_logs_file}
cd ..
