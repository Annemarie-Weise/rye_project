#!/bin/bash
# usage: ./calc_freq_vcftools.sh /path/to/folder

FOLDER="$1"

vcftools --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.recode.vcf" \
         --out "$FOLDER/data/vcftools/freq_vcf.90" \
         --freq



