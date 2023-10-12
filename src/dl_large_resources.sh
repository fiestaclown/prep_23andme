#!/bin/bash
threads_to_use=$1

wget -m ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/ -P resources/

mv  -v resources/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/* resources/ref_panel/phase

total_threads=$(nproc) # Get total number of CPU cores


echo "Total threads: $total_threads"
echo "Threads to use: $threads_to_use"

for ((i=1; i<=22; i+=threads_to_use)); do
    end=$((i+threads_to_use-1))
    [ $end -gt 22 ] && end=22 # Ensure the end value doesn't exceed 22

    echo "Processing chromosomes $i to $end"

    seq $i $end | parallel -j $threads_to_use "software/xcftools_static view -i resources/ref_panel/phase/1kGP_high_coverage_Illumina.chr{}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -o resources/ref_panel/impute/1kGP_high_coverage_Illumina.chr{}.filtered.SNV_INDEL_SV_phased_panel_xcf.bcf -O sh -r chr{} -T$threads_to_use"
done


wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P resources/

software/samtools faidx resources/GRCh38_full_analysis_set_plus_decoy_hla.fa

python3 software/gatk-4.4.0.0/gatk CreateSequenceDictionary R=resources/GRCh38_full_analysis_set_plus_decoy_hla.fa O=resources/GRCh38_full_analysis_set_plus_decoy_hla.dict

pip install gdown

gdown https://drive.google.com/uc?id=13jGKvtooQz9ASgtOrKjBgPfp4QvG6FYM -O resources/hg19.23andme.fa


software/samtools faidx resources/hg19.23andme.fa
 
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip

unzip -d software/ software/gatk-4.4.0.0.zip

rm software/gatk-4.4.0.0.zip