#!/bin/bash
# usage: ./calc_eigen_biallelic_plink2.sh /path/to/folder

FOLDER="$1"

plink2 --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.biallelic.vcf" \
    --pca 20 meanimpute --allow-extra-chr --double-id \
    --out "$FOLDER/data/pca_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.90.biallelic"

plink2 --vcf "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.biallelic.vcf" \
    --pca 20 meanimpute --allow-extra-chr --double-id \
    --out "$FOLDER/data/pca_PLINK/maf.05_minDP20_maxDP100_minQ40_missing.90.biallelic"
