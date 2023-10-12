#23andme data preperation
## Goal is to convert 23andme files to VCF without broken data or incorrect ref alleles

convert_twothreemee <- function(file, name, threads){
#create temp file for cleaned files and delete afterwards
if(!dir.exists(glue("data/temp"))){
  dir.create(glue("data/temp"))
}


# filter bad data ---------------------------------------------------------
#lot of non-acgt data
#ok genotypes
ok_genotypes = expand_grid(
  allele1 = c("A", "C", "G", "T"),
  allele2 = c("A", "C", "G", "T"),
) %>% mutate(
  genotype = glue("{allele1}{allele2}")
)



file_out <- file %>% str_replace("data/",'data/temp/')

  
  x = read_tsv(glue("{file}.txt"), comment = "#", col_names = F)
  
  #keep autosome only
  x2 = x %>% filter(X2 %in% (1:22))
  
  #throw away bad data
  x3 = x2 %>% filter(X4 %in% ok_genotypes$genotype)
  
  #remove dups
  x_chrpos = glue("{x3$X2}:{x3$X3}")
  x4 = x3[!duplicated(x_chrpos), ]
  
  x4
  #save again
  write_tsv(x4, file = glue("{file_out}.fixed.txt") , col_names = F)



# convert with bcftools ---------------------------------------------------
#https://samtools.github.io/bcftools/howtos/convert.html
#notice that they don't add chr prefix



  system(str_glue("bcftools convert --tsv2vcf {file_out}.fixed.txt -f resources/hg19.23andme.fa -s {name} -o {file_out}.vcf.gz"))
  


#bgzip to minimize filesize

  # system(str_glue("bgzip {file_out}.vcf.gz"))
   system(str_glue("bcftools index {file_out}.vcf.gz"))

# 
# meand23_vcf.gz_files = list.files('data',pattern = "genome_",full.names = T) %>% str_subset("fixed.vcf.gz$")



split_vcf_by_autosome_and_sex(input = glue("{file_out}.vcf.gz"))

tictoc::tic()
# Execute the shell script with the input VCF file as a parameter
system(glue("src/update_vcf.sh  {file_out} {threads}"))
tictoc::toc()

list.files("data/temp", full.names = T) %>% str_extract(".*_updateALT.vcf.gz$") %>% na.omit %>% gtools::mixedsort() %>% write_lines("data/temp/concat")

system(glue("bcftools concat -f data/temp/concat -o {file_out}_concat.vcf"))

# Remove variants where REF == ALT
call = glue("awk -F'\t' '{{ if ($1 ~ /^#/) print; else if ($4 != $5) print }}' {file_out}_concat.vcf > {file_out}_REFALT.vcf")
system(call)

#Finally remove indels
call = glue("plink2 --vcf {file_out}_REFALT.vcf  --snps-only just-acgt --recode vcf bgz --out {file_out}_final") 
system(call)
#index vcf


#check if vcf is valid
system(str_glue("vcf_validator -i {file_out}_final.vcf.gz"))

#move to valid or not valid folder 
sum_file <- list.files("data/temp", full.names = T) %>% str_extract('.*.errors_summary.*') %>% str_extract(glue(".*{file_out}.*")) %>% na.omit
vcf_sum <-  read_lines(sum_file)
if(str_detect(vcf_sum[1],'According to the VCF specification, the input file is valid')){
  system(glue("mv {file_out}_final.vcf.gz {str_replace(file_out,'data/temp/','data/valid/')}.vcf.gz"))
}else{
  system(glue("mv {file_out}_final.vcf.gz {str_replace(file_out,'data/temp/','data/not_valid/')}.vcf.gz"))
  system(glue("mv sum_file {str_replace(file_out,'data/temp/','data/not_valid/')}.vcf.gz"))
}

#nuke temp folder
unlink(x = "data/temp/*",recursive = T)


}



