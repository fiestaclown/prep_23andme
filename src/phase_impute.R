#combine,phase,impute - prep for scoring
#create new folder in valid for phasing results and imputation results
prep_data <- function(threads){
if(!dir.exists("data/valid/phase")){
  dir.create("data/valid/phase")
}

if(!dir.exists("data/valid/impute")){
  dir.create("data/valid/impute")
}

#set path and file prefix
#for phasing
file_out <- "data/valid/phase/full_dataset"

#for imputation
file_out_imp <- "data/valid/impute/full_dataset"

#combine all files in valid
vcf_files <- list.files("data/valid", full.names = T) %>% str_extract('.*vcf.gz$') %>% na.omit 
vcf_files %>%  write_lines("merge")

#first index files
for(i in vcf_files){
  call = glue("bcftools index -f {i}")
  system(call)
}


#merge all files and save to /phase
call = glue("bcftools merge -l merge -o {file_out}.vcf.gz")
system(call)

#index file
call = glue("bcftools index {file_out}.vcf.gz")
system(call)

#set env path for plugins
Sys.setenv("BCFTOOLS_PLUGINS" = "software/bcftools_plugins/")

#we need to add the AN,AC field for phasing 
call = glue("bcftools +fill-tags {file_out}.vcf.gz -o {file_out}_AC.vcf.gz  -- -t AN,AC ")
system(call)

#index file
call = glue("bcftools index {file_out}_AC.vcf.gz")
system(call)

#lift 
source("src/picard_liftover.R")
picard_liftover(input = glue("{file_out}_AC.vcf.gz"), output = glue("{file_out}_AC_hg38.vcf.gz"), target_reference_file = "resources/GRCh38_full_analysis_set_plus_decoy_hla.fa", chain_file = "resources/b37ToHg38.over.chain.gz")

#split by chromosome
source("src/split_vcf_by_autosome_and_sex.R")
split_vcf_by_autosome_and_sex(input = glue("{file_out}_AC_hg38.vcf.gz"), prefix = T)


#phase
source("src/SHAPEIT5.R")
for(i in 22:1){
phase_SHAPEIT5(phasing_panel = glue('resources/ref_panel/phase/1kGP_high_coverage_Illumina.chr{i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz'), chr = glue("chr{i}"), input_file = glue("{file_out}_AC_hg38_chr{i}.vcf.gz"), mode = 'phase_panel', maps = "resources/genetic_maps.b38/",output_file = glue("{file_out_imp}_AC_chr{i}.phased.vcf"), build = 38, threads = threads)
}


#impute
source("src/impute.R")
for(i in  22:1){
imputePhasedData(output_dir = file_out_imp,
                 target_dir = file_out_imp,
                 ref_dir = "resources/ref_panel/impute/",
                 map_dir = "resources/genetic_maps.b38/",
                 imp5Chunker = "software/imp5Chunker_static_file",
                 imp5Converter = "software/imp5Converter_static_file",
                 impute5 = "software/impute5_static_file",
                 chr = glue("{i}"))

}

#combine all files and convert to plink
list.files(file_out_imp %>% str_extract('.*/'), full.names = T) %>% str_extract(".*chr\\d+_imputed.vcf.gz$") %>% na.omit %>% gtools::mixedsort() %>% write_lines('data/valid/impute/merge')

#concat
call = glue("bcftools concat -a -f data/valid/impute/merge -o data/full_dataset_imputed.vcf.gz")
system(call)

#we might need to sort the data
call = glue("plink2 --vcf data/full_dataset_imputed.vcf.gz --sort-vars --set-all-var-ids @:# --rm-dup force-first --make-pgen --out data/full_dataset_imputed_sort")
system(call)


#convert to plink and filter for acceptable imputation level
call = glue('plink2 --pfile data/full_dataset_imputed_sort  --max-alleles 2  --extract-if-info "INFO >= 0.4" --make-bed --out data/full_dataset_imp_filter')
system(call)
}



