imputePhasedData <- function(output_dir, target_dir, ref_dir, map_dir, imp5Chunker, imp5Converter, impute5, chr) {


  # Chunking step
  chunk_cmd <- glue("{imp5Chunker} --h {ref_dir}1kGP_high_coverage_Illumina.chr{i}.filtered.SNV_INDEL_SV_phased_panel_xcf.bcf --g {target_dir}_AC_chr{chr}.phased.vcf.gz --r chr{chr} --o {output_dir}_chr{chr}_chunk_coordinates.txt --l {output_dir}chr{chr}_chunker.log")
  system(chunk_cmd)

  # Extracting imputation regions from coordinates file
  coords_file <- glue("{output_dir}_chr{chr}_chunk_coordinates.txt")
  if (!file.exists(coords_file)) stop(glue("Coordinates file not found"))

  regions <- read.table(coords_file, header = FALSE)[,4]
  write.table(regions, glue("{output_dir}chr{chr}_imputation_regions.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Imputation step
  for(buff in regions){
  impute_cmds <- glue("{impute5} --h {ref_dir}1kGP_high_coverage_Illumina.chr{i}.filtered.SNV_INDEL_SV_phased_panel_xcf.bcf --m {map_dir}chr{chr}.b38.gmap.gz --g {target_dir}_AC_chr{chr}.phased.vcf.gz --r {buff} --buffer-region {buff} --o {output_dir}_chr_{chr}_imputed.vcf.gz --l {output_dir}_chr_{i}_imputed.log")
system(impute_cmds)
  }


  # Create list of imputed files to be used in ligation step
  imputed_files <- list.files(path = output_dir %>% str_extract('.*/'), full.names = TRUE) %>% str_extract(pattern = glue(".*chr_{chr}.*imputed.vcf.gz$")) %>% na.omit %>% gtools::mixedsort()
  write_lines(imputed_files, glue("{output_dir}_chr{chr}_imputed_chunk_files.txt"))

  # Ligation step
  call <- glue("software/bcftools concat  --ligate -f {output_dir}_chr{chr}_imputed_chunk_files.txt  -o {output_dir}chr{chr}_imputed.vcf.gz")
  system(call)
  call = glue("bcftools index {output_dir}chr{chr}_imputed.vcf.gz") 
  system(call)

}

