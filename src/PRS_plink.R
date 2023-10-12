#calculate prs scores 

quick_prs <- function(prs_column,
                      pos_column,
                      allele_column, 
                      score_file = "resources/GWAS/CTPR_beta_coefficients_hg38.tsv",  
                      prs_score_file,
                      frequency_file,
                      dataset
                      ) { 
  


# Calculate PRS using plink2
call <- glue::glue("plink2 --bfile {dataset} --score {score_file}  {pos_column} {allele_column} {prs_column} 'header' --read-freq {frequency_file} --out {prs_score_file}")
print(call)
system(call)


}