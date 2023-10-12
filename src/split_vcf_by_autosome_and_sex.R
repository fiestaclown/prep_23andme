#split vcf by chr
library(tidyverse)
library(glue)

split_vcf_by_autosome_and_sex <- function(input, output = '', prefix = F, unzip = F){
  #if output is empty expect output should be written to the same folder as the input
  if(output == ''){output = input}
 #create vector with chromosome contigs
  chr <- 1:22 %>% as.character %>% c(.,'X','Y')
  if(prefix == T){
    chr <- 1:22 %>% as.character %>% c(.,'X','Y') %>% str_c('chr',.)  
    
    if(unzip == T){
      for(i in chr){
        #extract chromosome i from vcf
        call = glue("bcftools view {input} -r {i} -o {str_replace(output, '.vcf.gz',glue('_{i}.vcf'))}")
        system(call)
      }
      }else{
        for(i in chr){
    #extract chromosome i from vcf
    call = glue("bcftools view {input} -r {i} -o {str_replace(output, '.vcf.gz',glue('_{i}.vcf.gz'))}")
    system(call)
    
    #index file
    call = glue("bcftools index -f {str_replace(output, '.vcf.gz',glue('_{i}.vcf.gz'))}")
    system(call)
    }
        }
      
  
  }else{
    if(unzip == T){
      for(i in chr){
        #extract chromosome i from vcf
        call = glue("bcftools view {input} -r {i} -o {str_replace(output, '.vcf.gz',glue('_chr{i}.vcf'))}")
        system(call)
      }
      }else{
    
  for(i in chr){
    #extract chromosome i from vcf
    call = glue("bcftools view {input} -r {i} -o {str_replace(output, '.vcf.gz',glue('_chr{i}.vcf.gz'))}")
    system(call)
    
    #index file
    call = glue("bcftools index -f {str_replace(output, '.vcf.gz',glue('_chr{i}.vcf.gz'))}")
    system(call)
    }
        }
      }
  }
