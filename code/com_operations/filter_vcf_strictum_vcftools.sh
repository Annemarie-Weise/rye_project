#!/bin/bash
# usage: ./filter_vcf.sh /path/to/folder

FOLDER="$1"

vcftools \
  --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.recode.vcf" \
  --keep "$FOLDER/data/ID_data/Secale_strictum_IDs.txt" \
  --mac 1 \
  --recode --recode-INFO-all \
  --out "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.strictum"

