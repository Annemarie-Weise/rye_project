#!/bin/bash
# usage: ./filter_vcf.sh /path/to/folder

FOLDER="$1"

vcftools --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.biallelic.vcf" \
         --out "$FOLDER/data/vcftools/freq_biallelic_vcf" \
         --freq



