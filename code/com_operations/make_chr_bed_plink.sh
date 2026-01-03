#!/bin/bash
# usage: ./make_chr_bed_plink.sh /path/to/folder

FOLDER="$1"

plink --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.splitted_chr.vcf" \
    --make-bed \
    --out "$FOLDER/data/bed_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.90.splitted_chr" \
    --double-id \
    --allow-extra-chr

plink --vcf "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.splitted_chr.vcf" \
    --make-bed \
    --out "$FOLDER/data/bed_PLINK/maf.05_minDP20_maxDP100_minQ40_missing.90.splitted_chr" \
    --double-id \
    --allow-extra-chr
