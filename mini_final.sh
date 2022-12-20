#!/bin/bash
################################################################################
# Data preparation and analysis for predicting ancestry
# 
# Created by David de la Cerda
# script name: mini_final.sh
# 
# 
# 
# input: VCF and ID files provided by Dr. Tiao
# output:  raw and filtered stats, binary files from Plink used for PCA visualization, predicted ancestry text file 
# required software: python 3.7 samtools 1.10, R 4.0.4, bcftools-1.16, VCFtools 0.1.16, and fraposa
# fraposa is hosted by Github user: daviddaiweizhang and its Python dependencies are listed
# fraposa was cloned at: https://github.com/daviddaiweizhang/fraposa.git
################################################################################
#checking for unfiltered variants
bcftools view -H acs_mini_project.vcf.bgz | wc -l #2203614 variants unfiltered
# index subset vcf in case it's needed
bcftools index acs_mini_project.vcf.bgz
#make directory for raw stats
mkdir vcf_stats
RAW_VCF=/Users/David/brd_mini/acs_mini_filled.vcf.bgz
OUT=/Users/David/brd_mini/vcf_stats_unfiltered/mini_project
#calculate allele freq, only include bi-allelic sites
vcftools --gzvcf $RAW_VCF --freq2 --out $OUT --max-alleles 2
#quality scores
vcftools --gzvcf $RAW_VCF --site-quality --out $OUT
#calculate proportion of missing data per individ
vcftools --gzvcf $RAW_VCF --missing-indv --out $OUT
#missing data per site
vcftools --gzvcf $RAW_VCF --missing-site --out $OUT
#calculate heterozygosity 
vcftools --gzvcf $RAW_VCF --het --out $OUT
mkdir /Users/David/brd_mini/vcf_stats_filtered
#run stats on filtered data
FILTER_VCF=/Users/David/brd_mini/acs_mini_filtered.vcf.bgz
FILTER_OUT=/Users/David/brd_mini/vcf_stats_filtered/mini_project
#calculate allele freq, only include bi-allelic sites
vcftools --gzvcf $FILTER_VCF --freq2 --out $FILTER_OUT --max-alleles 2
vcftools --gzvcf $FILTER_VCF --site-quality --out $FILTER_OUT
#calculate proportion of missing data per individ
vcftools --gzvcf $FILTER_VCF --missing-indv --out $FILTER_OUT
#missing data per site
vcftools --gzvcf $FILTER_VCF --missing-site --out $FILTER_OUT
#calculate heterozygosity 
vcftools --gzvcf $FILTER_VCF --het --out $FILTER_OUT
#after visualizing in mini_acs_PCA.R and from guidance from Dr. Tiao instructions
#using BCF tools for filtering
VCF_IN=/deac/bio/peaseGrp/ddelacer/brd_mini/acs_mini_project.vcf.bgz
VCF_OUT=/deac/bio/peaseGrp/ddelacer/brd_mini/acs_mini_project_filtered.vcf.bgz
#fill
bcftools +setGT $VCF_IN -Oz -o acs_mini_filled.vcf.bgz -- -t . -n 0
#add allele freq
bcftools +fill-tags acs_mini_filled.vcf.bgz -Ov -o acs_mini_filled_alfq.vcf.bgz -- -t AF
# set filters
#applying filters, removing F missing gt .1
bcftools view -i 'F_MISSING < 0.1' acs_mini_filled_alfq.vcf.bgz -o acs_mini_missing.vcf.bgz
#MAF threshold based on what was asked to filter on
bcftools view -i 'MAF > 0.01' acs_mini_missing.vcf.bgz -Ov -o acs_mini_maf.vcf.bgz
#monomorphic sites removed
bcftools view -c 1 acs_mini_maf.vcf.bgz -Ov -o acs_mini_mono.vcf.bgz
#multiallelic sites removed
bcftools norm -d all acs_mini_mono.vcf.bgz -Ov -o acs_mini_filtered.vcf.bgz
#Subset the filtered data by known and unknown individuals
#these data are needed for prediction tool fraposa
#study_ids.txt and reference_ids.txt generated using mini_acs_PCA.R
bcftools view -S study_ids.txt acs_mini_filtered.vcf.bgz > acs_mini_filtered_study.vcf.bgz
bcftools view -S reference_ids.txt acs_mini_filtered.vcf.bgz > acs_mini_filtered_reference.vcf.bgz
###Perform PCA and pruning in Plink 
#PCA with unknown ancestry
#performing pruning in window of 50 with .1 R-squared
#allowing for extra chromosomes in case any unknown but informative sites from unknown contigs included
VCF=/Users/David/brd_mini/acs_mini_filtered.vcf.bgz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out miniproject
#PCA step
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--extract miniproject.prune.in \
--make-bed --pca --out acsminiproject
###Using prediction tool fraposa to estimate ancestry
#plink steps
VCF_study=/Users/David/brd_mini/predict_ancestry_pruned/acs_mini_filtered_pruned_study.vcf.bgz
VCF_ref=/Users/David/brd_mini/predict_ancestry_pruned/acs_mini_filtered_pruned_reference.vcf.bgz
plink --vcf $VCF_study --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--make-bed --out study
plink --vcf $VCF_ref --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--make-bed --out reference
python3 fraposa/fraposa_runner.py --stu_filepref study reference
#predstupopu file generated using R, it is individual, family, and ancestry of known individuals, it generates a file called stupref.pcs which assigns ancestry and gives percent assignment to that ancestry
python3 fraposa/predstupopu.py reference study

