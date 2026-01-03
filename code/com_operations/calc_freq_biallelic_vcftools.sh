#!/bin/bash
# usage: ./calc_freq_biallelic_vcftools.sh /path/to/folder

FOLDER="$1"

vcftools --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.biallelic.vcf" \
         --out "$FOLDER/data/pca_PLINK/freq_biallelic.90_vcf" \
         --freq



