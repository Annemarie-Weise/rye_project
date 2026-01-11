#!/bin/bash

#SBATCH --job-name=vcf_filtering1                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=1                               # number of cpus for this task (default=1)
#SBATCH --mem=20G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job/job_%j.out                             # File to which standard out will be written
#SBATCH --error=job/job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=weisea@ipk-gatersleben.de  # Email to which notifications will be sent

# load your module
module load vcftools


for K in {1..5}; do
        MAF=$(echo "0.01 * $K" | bc)
        for N in {90..95}; do
                MISS=$(echo "0.01 * $N" | bc)
                vcftools --gzvcf ~/rye_evg2023t0_divpanel2022_schreiber2022_ER22_q30Q30_DPAD_Q40.ANN.vcf.gz \
                        --maf $MAF --min-meanDP 20 --max-meanDP 100 --minQ 40 --max-missing $MISS \
                        --out ~/maf${MAF}_minDP20_maxDP100_minQ40_missing${MISS} \
                        --recode --recode-INFO-all
        done
done

