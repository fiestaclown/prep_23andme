#Main script
# Get the arguments passed to the script
args <- commandArgs(trailingOnly = TRUE)
red_cores <- as.integer(args[1])
install_logic <- as.integer(args[2])

#download phasing panel and reference fasta files
#get number of avaible threads
threads_all <- as.integer(system("nproc", intern = T))
threads_to_use <- threads_all - red_cores

#0 make scripts usable and check for package installations
system("chmod +x src/dl_large_resources.sh")
system("chmod +x src/update_vcf.sh")

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(tidyverse, glue, gtools, rms, tictoc , rvest, data.table)

#1 optional if needed
#download necessary resources already prepared for the scripts
#those include reference files needed for data preparation. Unpacked this will be around 100GB
if(install_logic == 1){

call = (glue("sh src/dl_large_resources.sh {threads_to_use}"))
system(call)

#download and prep newest dbsnp version 
source("src/annotate_contig_dbsnp_37.R")
}






#2
#start with 23andme text file conversion to valid vcf
#files
meand23_files = list.files( 'data', pattern = "genome_",full.names = T) %>% str_subset("Full_\\d{14,}.txt$") %>% str_remove('.txt')


#unique sample names
names <- str_c("andme_0",1:length(meand23_files))

#save sample_names and file_names combo
tibble(file_names = meand23_files, names = names) %>% write_csv('files_sample_names.dict')
source("src/split_vcf_by_autosome_and_sex.R")
source("src/convert_23andme.R")
for(i in seq_along(meand23_files)){
  tictoc::tic()
  convert_twothreemee(file = meand23_files[i], name = names[i], threads = threads_to_use)
  tictoc::toc()
}

#3 
#prep for scoring - lift, phase, impute
source("src/phase_impute.R")
prep_data(threads = threads_to_use)

#4
#calculate scores and frequencies
source("src/23andme_stats.R")