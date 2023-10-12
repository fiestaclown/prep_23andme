#!/bin/bash

#--------Paths to files and directories---------
input_file=data/valid/combined.vcf.gz	#edit  bcf/vcf format
output_dir=data/valid/combined.phased.vcf.gz	#edit
shapeit4=software/shapeit4	#edit
map_dir=resources/genetic_maps.b37/ 	#edit


#Input file need to be indexed before phasing 
tabix -p vcf ${input_file}

#Phasing is done for each chromosome seperately so edit arguments accordingly.    
${shapeit4} \
	--input ${input_file} \
	--map ${map_dir}chr#.b37.gmap.gz \
	--region chr# \
	--output ${output_dir}chr#_phased_file.vcf.gz \
	--log ${output_dir}chr#_phased_file.log \

#Indexing phased files for imputation	
tabix -p vcf ${output_dir}chr#_phased_file.vcf.gz 