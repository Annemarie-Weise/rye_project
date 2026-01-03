#!/bin/bash
# usage: ./split_multiallelic_sites_in_vcf_bcftools.sh /path/to/folder

FOLDER="$1"

bcftools norm -m \
    -any "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.recode.vcf" \
    -o "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.splitted.vcf"

bcftools norm -m \
    -any "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.recode.vcf" \
    -o "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.splitted.vcf"

