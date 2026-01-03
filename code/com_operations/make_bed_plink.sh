#!/bin/bash
# usage: ./make_bed_plink.sh /path/to/folder

FOLDER="$1"

plink --vcf "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.splitted.vcf" \
    --make-bed \
    --out "$FOLDER/data/bed_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.90.splitted" \
    --double-id \
    --allow-extra-chr
