#!/bin/bash
# usage: ./filter_vcf.sh /path/to/folder

FOLDER="$1"

plink --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted.vcf" \
      --pca \
      --double-id \
      --allow-extra-chr \
      --out "$FOLDER/data/pca_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted"



