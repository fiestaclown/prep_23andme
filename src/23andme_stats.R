#23andme statistics

#load librarys
library(tidyverse)
library(data.table)
library(glue)
source("src/PRS_plink.R")
GWAS <- "resources/GWAS/"

#create a folder with information that will back to me for evaluation
if(!dir.exists("data/output")){
  dir.create("data/output")
}
###################################### frequencies and geno counts ##################################################################
#lets start with allele frequencies for both the imputed and non-imputed variants
#pre-imp
#ensure that no duplicates exist and convert to bim
call = glue("plink2 --vcf data/valid/phase/full_dataset.vcf.gz --max-alleles 2  --rm-dup 'force-first' --set-all-var-ids @:# --make-bed --out data/pre_imp_dup")
system(call)


call = glue("plink2 --bfile data/pre_imp_dup --freq  --out data/output/pre_imp_dup")
system(call)

call = glue("plink2 --bfile data/pre_imp_dup --geno-counts --rm-dup 'force-first' --out data/output/pre_imp_dup")
system(call)

#post-imp qual filter
call = glue("plink2 --bfile data/full_dataset_imp_filter --freq  --out data/output/full_dataset_imp_filter")
system(call)

call = glue("plink2 --bfile data/full_dataset_imp_filter --geno-counts  --out data/output/full_dataset_imp_filter")
system(call)
###################################### PRS #####################################################################################

## Lee et al models 
#just with plink
#loop across all the variants in Lee et al, score each, and then read all scores to a single file for ease of use
#load Lee's combined file
lee2018_combined10k = read_tsv(glue("{GWAS}COMBINED.to10K.txt"),
                               #have to explicitly give the types, otherwise it gusses these are dates
                               col_types = cols(chrpos37 = col_character(),
                                                chrpos38 = col_character()))
#determine names of cols of betas
(lee2018_models = lee2018_combined10k %>% names() %>% str_subset("Beta") %>% str_replace("Beta_", ""))

#fill in 0's for NA's
if (F) {
  beta_cols = names(lee2018_combined10k) %>% str_subset("Beta_")
  for (v in beta_cols) {
    lee2018_combined10k[[v]] = lee2018_combined10k[[v]] %>% mapvalues(from = NA, to = 0)
  }
  
  #ensure no NA
  assert_that(!anyNA(lee2018_combined10k[beta_cols]))
  
  #write back
  lee2018_combined10k %>% write_tsv(glue("{gwas}COMBINED_10K_fixed.tsv"))
}


#pre_imp
#EA
quick_prs(prs_column = 19,
          allele_column = 4,
          pos_column = 31,
          dataset = "data/pre_imp_dup",
          score_file = glue("{GWAS}COMBINED_10K_fixed.tsv"),
          frequency_file ="data/output/pre_imp_dup.gcount",
          prs_score_file = 'data/output/pre_imp_dup_EA'
)
#height
quick_prs(prs_column = 15,
          allele_column = 6,
          pos_column = 16,
          dataset = "data/pre_imp_dup",
          score_file = glue("{GWAS}CTPR_beta_coefficients_hg38.tsv"),
          frequency_file ="data/output/pre_imp_dup.gcount",
          prs_score_file = 'data/output/pre_imp_dup_height'
)

#post_imp
#EA
quick_prs(prs_column = 19,
          allele_column = 4,
          pos_column = 32,
          dataset = "data/full_dataset_imp_filter",
          score_file = glue("{GWAS}COMBINED_10K_fixed.tsv"),
          frequency_file ="data/output/full_dataset_imp_filter.gcount",
          prs_score_file = 'data/output/full_dataset_imp_EA'
)
#height
quick_prs(prs_column = 15,
allele_column = 6,
pos_column = 2,
dataset = "data/full_dataset_imp_filter",
score_file = glue("{GWAS}CTPR_beta_coefficients_hg38.tsv"),
frequency_file ="data/output/full_dataset_imp_filter.gcount",
prs_score_file = 'data/output/full_dataset_imp_height'
)

