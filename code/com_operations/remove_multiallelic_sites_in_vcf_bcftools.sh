#!/bin/bash
# usage: ./remove_multiallelic_sites_in_vcf_bcftools.sh /path/to/folder

FOLDER="$1"

bcftools view -m2 -M2 \
    "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.recode.vcf" \
    -o "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.biallelic.vcf"

bcftools view -m2 -M2 \
    "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.recode.vcf" \
    -o "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.biallelic.vcf"
