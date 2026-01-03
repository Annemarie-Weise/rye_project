#!/bin/bash
# usage: ./filter_vcf_non-sylvestre_vcftools.sh /path/to/folder

FOLDER="$1"

vcftools \
  --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.recode.vcf" \
  --remove "$FOLDER/data/ID_data/Secale_sylvestre_IDs.txt" \
  --recode --recode-INFO-all \
  --mac 1 \
  --out "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.non-sylvestre"


