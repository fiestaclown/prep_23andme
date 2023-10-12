#picard liftover
library(glue)
picard_liftover <- function(input, output, target_reference_file, chain_file, rejected_variants = 'rejected_variants.vcf', recover_ref_swap = 'true', options ='' ){
  #ensure that input is not equal to output
  if(input == output){
    stop('Input filename equals output filename. Original file would get destroyed! Please use an unique output filename.')
  }  
  
  #try liftover with picard
  call = glue('python3 software/gatk-4.4.0.0/gatk LiftoverVcf I={input} O={output} CHAIN={chain_file} REJECT={rejected_variants} R={target_reference_file} RECOVER_SWAPPED_REF_ALT={recover_ref_swap} {options}')
  system(call)  
  
}

