#!/bin/bash
# usage: ./create_hapmap_tassel.sh /path/to/folder

FOLDER="$1"

run_pipeline.pl -Xmx5g -fork1 -vcf \
                "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.90.biallelic.vcf" \
                -sortPositions \
                -ImputationPlugin -ByMean true -endPlugin \
                -export "$FOLDER/data/hmp_TASSEL/mygenotypes_maf0.01" \
                -exportType Hapmap

run_pipeline.pl -Xmx5g -fork1 -vcf \
                "$FOLDER/data/vcftools/maf.05_minDP20_maxDP100_minQ40_missing.90.biallelic.vcf" \
                -sortPositions \
                -ImputationPlugin -ByMean true -endPlugin \
                -export "$FOLDER/data/hmp_TASSEL/mygenotypes_maf0.05" \
                -exportType Hapmap



