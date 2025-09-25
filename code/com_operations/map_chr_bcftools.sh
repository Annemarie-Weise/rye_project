#!/bin/bash
# usage: ./filter_vcf.sh /path/to/folder

FOLDER="$1"

bcftools annotate \
         --rename-chrs "$FOLDER/data/ID_data/map_chr.txt" \
         -o "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted_chr.vcf" \
         "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted.vcf"
#tabix -p vcf output.vcf.gz

