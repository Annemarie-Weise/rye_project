#!/bin/bash
# usage: ./find_annotation.sh /path/to/folder

FOLDER="$1"
SNP_LIST="$FOLDER/data/ID_data/only_sylvestre_biallelic_SNPs_maf.01.txt"
VCF="$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.recode.vcf"

REGIONS="$FOLDER/data/vcftools/snps.txt"

# chr1R_311940740  -> chr1R:311940740-311940740
awk -F'_' '{print $1 ":" $2 "-" $2}' "$SNP_LIST" > "$REGIONS"

bcftools view -R "$REGIONS" "$VCF" \
  | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ANN\n' \
  > "$FOLDER/results/annotation_snps.tsv"
