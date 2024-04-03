#!/bin/bash

input_lnc_dir="/data/lyc_data/lyc_data/AS_Ribo_RAW/GDRT22110123/GDRT22110123_std_2-lnc/cleandata/"
input_ribo_dir="/data/lyc_data/lyc_data/AS_Ribo_RAW/GDRT22110123/GDRT22110123_std_3-riboseq/cleandata/"
output_dir="/data/lyc_data/lyc_data/AS_Ribo_RAW/GDRT22110123/circpro_out/"

for file in $input_lnc_dir/*_1.clean.fq.gz; do
    filename=$(basename "$file" _1.clean.fq.gz)
    output_prefix="$output_dir${filename%_*}_plus"
    mkdir -p "$output_prefix"

    lnc_1="${input_lnc_dir}${filename}_1.clean.fq.gz"
    lnc_2="${input_lnc_dir}${filename}_2.clean.fq.gz" 
    ribo_1="${input_ribo_dir}${filename}_1.clean.fq.gz"
    ribo_2="${input_ribo_dir}${filename}_2.clean.fq.gz"

    echo "Running on: $lnc_1 and $lnc_2,\n Aligning on: $ribo_1 and $ribo_2,\n Output to: $output_prefix"
    perl circProplus.pl -c "$lnc_1,$lnc_2" -o "$output_prefix" -g /data/ref_data/mmu/mm10/Annotation/mm10_ref.gtf \
        -ref /data/ref_data/mmu/mm10/Index/circpro/mm10.fa -m ve -r "$ribo_1,$ribo_2" -rRNA /data/ref_data/mmu/mm10/fa/filtered.fa -T 32 -D
done
