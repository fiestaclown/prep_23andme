#goal is to rename the contigs of dbsnp GCF_000001405.40

#set download path and wd for dbsnp
#download latest dbsnp for hg38 
system(glue("wget ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.25.gz -P resources/dbsnp_hg37_contig/"))
system(glue("wget ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.25.gz.gz.tbi -P resources/dbsnp_hg37_contig/"))


#first download Assembly Unit table from  https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25
html <- rvest::read_html("https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25")

#extract table with contig labels from page and convert to tibble and select only relevant columns
contig_table <- html %>%
  rvest::html_element(css = '#asm_Primary_Assembly > div:nth-child(2) > div:nth-child(1) > table:nth-child(1)') %>%
  rvest::html_table() %>% select(1,4)


#select dbsnp contig name and add mutate new contig name
contig_table_new_contig <- contig_table %>% rename('Molecule_name' = `Molecule name` , old_contig = `RefSeq sequence`) %>% dplyr::mutate(new_contig = str_c(str_remove(Molecule_name,'Chromosome '))) %>% select(-Molecule_name)

#get rid of last line which contains no useful information (n/a)
contig_table_new_contig_final <- contig_table_new_contig %>% filter(row_number() <= n()-1)

#write lines with contig changes
fwrite(contig_table_new_contig_final, 'resources/dbsnp_hg37_contig/dbsnp_hg37_contig.txt', col.names = F, sep =' ')


#adjust contigs in dbnsp
system("bcftools annotate resources/dbsnp_hg37_contig/GCF_000001405.25.gz --rename-chrs resources/dbsnp_hg37_contig/dbsnp_hg37_contig.txt -o resources/dbsnp_hg37_contig/GCF_000001405.25.gz_hg37_contigs.vcf.gz")

#remove old file
file.remove("resources/dbsnp_hg37_contig/GCF_000001405.25.gz")

#split by chromosome
for(i in 1:22){
  #split by chromosome, convert to pvar
   call = glue("plink2 --vcf resources/dbsnp_hg37_contig/GCF_000001405.25.gz_hg37_contigs.vcf.gz --chr {i} --recode vcf bgz --out resources/dbsnp_hg37_contig/GCF_000001405.25_hg37_contigs_chr{i}")
   system(call)
   call = glue("bcftools index resources/dbsnp_hg37_contig/GCF_000001405.25_hg37_contigs_chr{i}.vcf.gz")
}

#remove not splitted version
file.remove("resources/dbsnp_hg37_contig/GCF_000001405.25.gz_hg37_contigs.vcf.gz")