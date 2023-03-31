#!/bin/bash

# DESCRIPTION: simple QC of MoBa Genotype data; based on
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

if [ ! -f 22_combined.bed ] # if there is no combined .bed file
then
  # EXTRACT SAMPLES - taking only the control group: women without epilepsy
  #   who did not take folate supplements
  for i in {1..22}
  do
    # create links to the original files
    ln -s ${moba_plink_dir}/${i}.* .
    ./plink --bfile $i --keep mothers_IDs_4plink.txt --maf 0.01 --hwe 1e-6 --geno 0.01 --make-bed --out ${i}_SNP_QC
    
    if [ $? -ne 0 ]
    then
      cd $working_dir
      exit 1
    fi
    
    # merge one by one from chr.2
    if [ $i -gt 1 ]
    then
      prev_files=$(echo "$i - 1" | bc -l)
      if [ ! -f ${prev_files}_combined ] # if there is no such file, this is first merging
      then
        # create symlinks
        ln -s ${prev_files}_SNP_QC.bed ${prev_files}_combined.bed
        ln -s ${prev_files}_SNP_QC.bim ${prev_files}_combined.bim
        ln -s ${prev_files}_SNP_QC.fam ${prev_files}_combined.fam
      fi
      #perform merging
      ./plink --bfile ${i}_SNP_QC --bmerge ${prev_files}_combined --make-bed --out ${i}_combined
      if [ $? -ne 0 ]
      then
        cd $working_dir
        exit 1
      fi
      # delete unnecessary copies
      rm ${prev_files}_combined.b* ${i}_SNP_QC.b*
    fi
  done
fi

# synchronise newest data (one should replace '~/DATA' with a correct path)
rsync -auvz 22_combined.* ~/DATA/
rsync -auvz ${out_fname_base}.* ~/DATA/

# get back to R-project directory
cd $working_dir

cd DATA
touch ${all_logs_file}
cat ${out_fname_base}.log >> ${all_logs_file}
cd ..
