#!/bin/bash
# usage: ./pop_structure_admixture.sh /path/to/folder

FOLDER="$1"

MAFS=("0.05") #"0.01"

for MAF in "${MAFS[@]}"; do
    MAF_TAG="${MAF#0}"

    BED_FILE="$FOLDER/data/bed_PLINK/maf${MAF_TAG}_minDP20_maxDP100_minQ40_missing.90.splitted_chr.bed"

    OUT_DIR="$FOLDER/results/admixture/maf_${MAF}/output"
    LOG_DIR="$FOLDER/results/admixture/maf_${MAF}/logs"
    mkdir -p "$OUT_DIR" "$LOG_DIR"
    cd "$OUT_DIR"

    for K in {2..6}; do
        LOG_FILE="$LOG_DIR/log_K${K}.out"
        admixture --cv "$BED_FILE" $K | tee "$LOG_FILE"
    done
done
