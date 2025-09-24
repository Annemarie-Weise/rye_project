#!/bin/bash
# usage: ./filter_vcf.sh /path/to/folder

FOLDER="$1"

run_pipeline.pl -Xmx5g -fork1 -vcf \
                "$FOLDER/data/vcftools/maf.01_minDP20_maxDP100_minQ40_missing.10.biallelic.vcf" \
                -sortPositions \
                -ImputationPlugin -ByMean true -endPlugin \
                -export "$FOLDER/data/hmp_TASSEL/mygenotypes" \
                -exportType Hapmap



