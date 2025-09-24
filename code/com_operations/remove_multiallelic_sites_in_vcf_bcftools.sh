#!/bin/bash
# usage: ./filter_vcf.sh /path/to/folder

FOLDER="$1"

bcftools \
    view -m2 -M2 \
    "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.recode.vcf" \
    -o "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.biallelic.vcf"
