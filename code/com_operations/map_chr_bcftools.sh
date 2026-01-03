#!/bin/bash
# usage: ./map_chr_bcftools.sh /path/to/folder

FOLDER="$1"

bcftools annotate \
         --rename-chrs "$FOLDER/data/ID_data/map_chr.txt" \
         -o "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.splitted_chr.vcf" \
         "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.splitted.vcf"

bcftools annotate \
         --rename-chrs "$FOLDER/data/ID_data/map_chr.txt" \
         -o "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.splitted_chr.vcf" \
         "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.splitted.vcf"
