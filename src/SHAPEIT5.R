phase_SHAPEIT5 <- function(phasing_panel, chr, input_file, output_file, mode = 'phase_panel', pedigree_file = '', maps = '', SHAPEIT5 = '', build = '38', threads){
 if(mode == 'phase_panel'){
   #phase common variants with phase panel
  call = glue("software/phase_common_static --input {input_file} --region {chr} {maps}chr{chr}.b{build}.gmap.gz --reference {phasing_panel} --output {str_replace(output_file,'vcf','bcf')}  --thread {threads}")
  system(call)
}else if(mode == 'phase_panel_ped'){
#phase common variants with phase panel and pedigree info
system(glue("software/phase_common_static --input {input_file} --region chr{chr} --pedigree {pedigree_file} --map {maps}chr{chr}.b{build}.gmap.gz --reference {phasing_panel}  --output {str_replace(output_file,'vcf','bcf')} --thread {threads}"))
} else if(mode == 'samples'){
#phase common variants without phasing panel but at least 50 samples
system(glue("software/phase_common_static --input {input_file} --region chr{chr} --map {maps}chr{chr}.b{build}.gmap.gz --output {str_replace(output_file,'vcf','bcf')} --thread {threads}"))
}else if(mode == 'samples_ped'){
#phase common variants without phasing panel but at least 50 samples with pedigree information
system(glue("software/phase_common_static --input {input_file} --pedigree {pedigree_file} --region chr{chr} --map {maps}chr{chr}.b{build}.gmap.gz --output {str_replace(output_file,'vcf','bcf')} --thread {threads}"))
}
#convert to vcf.gz
system(glue("software/plink2 --bcf {str_replace(output_file,'vcf','bcf')} --recode vcf bgz --output-chr chrM --out {str_remove(output_file,'.vcf')} "))

#create index
system(glue("software/bcftools index {output_file}.gz"))
}

