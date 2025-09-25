#!/bin/bash
# usage: ./filter_vcf.sh /path/to/folder

FOLDER="$1"

OUT_DIR="$FOLDER/results/admixture/output"
LOG_DIR="$FOLDER/results/admixture/logs"
mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

for K in {2..4}; do
    LOG_FILE="$LOG_DIR/log${K}.out"
    OUTPUT_FILE="$OUT_DIR/admixture_output.${K}"
    admixture --cv "$FOLDER/data/bed_PLINK/maf.01_minDP20_maxDP100_minQ40_missing.10.splitted_chr.bed" $K | tee "$LOG_FILE" > "$OUTPUT_FILE"
done
