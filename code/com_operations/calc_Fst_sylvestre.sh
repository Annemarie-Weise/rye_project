#!/bin/bash
# usage: ./calc_Fst_sylvestre.sh /path/to/folder

FOLDER="$1"

MAFS=("01" "05")
MBS=(25 1)

POP_SYLVESTRE="$FOLDER/data/ID_data/Secale_sylvestre_IDs.txt"
POP_OTHER="$FOLDER/data/ID_data/Secale_other_plus_strictum_IDs.txt"

OUTDIR="$FOLDER/results/vcftools"

for maf in "${MAFS[@]}"; do
  VCF="$FOLDER/data/vcftools/maf.${maf}_minDP20_maxDP100_minQ40_missing.90.recode.vcf"

  for mb in "${MBS[@]}"; do
    WIN=$((mb * 1000000))
    TAG="${mb}Mb"

    #sylvestre vs other
    vcftools \
      --vcf "$VCF" \
      --weir-fst-pop "$POP_SYLVESTRE" \
      --weir-fst-pop "$POP_OTHER" \
      --fst-window-size "$WIN" \
      --fst-window-step "$WIN" \
      --out "$OUTDIR/Fst_sylvestre_other_maf0.${maf}_${TAG}"
  done
done
