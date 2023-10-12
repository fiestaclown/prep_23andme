#!/bin/bash

input_vcf=$1
threads_to_use=$2    
# A function to execute bcftools annotate on multiple files
annotate_vcf() {
    chr=$1
    input_vcf=$2
    bcftools annotate -i 'ALT="."' -k --collapse all \
        -a resources/dbsnp_hg37_contig/GCF_000001405.25_hg37_contigs_chr${chr}.vcf.gz \
        -c ALT ${input_vcf}_chr${chr}.vcf.gz \
        -o ${input_vcf}_chr${chr}_updateALT.vcf.gz
}

# Export the function so it is available to GNU parallel
export -f annotate_vcf

input_vcf_base=$1

# Run the function in parallel for each chromosome
total_threads=$(nproc)  # Get total number of CPU cores

echo "Total threads: $total_threads"
echo "Threads to use: $threads_to_use"

for ((i=1; i<=22; i+=threads_to_use)); do
    end=$((i+threads_to_use-1))
    [ $end -gt 22 ] && end=22  # Ensure the end value doesn't exceed 22

    echo "Processing chromosomes $i to $end"

    seq $i $end | parallel -j $threads_to_use "annotate_vcf {} ${input_vcf_base}"
done