#!/bin/bash
# usage: ./calc_Fst_3_pop.sh /path/to/folder

FOLDER="$1"

MAFS=("01" "05")
MBS=(25 1)

POP_STRICTUM="$FOLDER/data/ID_data/Secale_strictum_IDs.txt"
POP_SYLVESTRE="$FOLDER/data/ID_data/Secale_sylvestre_IDs.txt"
POP_OTHER="$FOLDER/data/ID_data/Secale_other_IDs.txt"

OUTDIR="$FOLDER/results/vcftools"

for maf in "${MAFS[@]}"; do
  VCF="$FOLDER/data/vcftools/maf.${maf}_minDP20_maxDP100_minQ40_missing.90.recode.vcf"

  for mb in "${MBS[@]}"; do
    WIN=$((mb * 1000000))
    TAG="${mb}Mb"

    echo "Running FST scans: maf.${maf}, window=${TAG} (${WIN} bp)"

    # 1) strictum vs sylvestre
    vcftools \
      --vcf "$VCF" \
      --weir-fst-pop "$POP_STRICTUM" \
      --weir-fst-pop "$POP_SYLVESTRE" \
      --fst-window-size "$WIN" \
      --fst-window-step "$WIN" \
      --out "$OUTDIR/Fst_strictum_sylvestre_maf0.${maf}_${TAG}"

    # 2) strictum vs other
    vcftools \
      --vcf "$VCF" \
      --weir-fst-pop "$POP_STRICTUM" \
      --weir-fst-pop "$POP_OTHER" \
      --fst-window-size "$WIN" \
      --fst-window-step "$WIN" \
      --out "$OUTDIR/Fst_strictum_other_maf0.${maf}_${TAG}"

    # 3) sylvestre vs other
    vcftools \
      --vcf "$VCF" \
      --weir-fst-pop "$POP_SYLVESTRE" \
      --weir-fst-pop "$POP_OTHER" \
      --fst-window-size "$WIN" \
      --fst-window-step "$WIN" \
      --out "$OUTDIR/Fst_sylvestre_other_minus_strictum_maf0.${maf}_${TAG}"
  done
done

