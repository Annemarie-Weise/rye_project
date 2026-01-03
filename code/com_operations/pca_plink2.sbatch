#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --job-name=pca_job                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=8                               # number of cpus for this task (default=1)
#SBATCH --mem=20G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job/job_%j.out                             # File to which standard out will be written
#SBATCH --error=job/job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=FAIL,END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=weisea@ipk-gatersleben.de # Email to which notifications will be sent

# init module environment
. /etc/profile.d/modules.sh

# load your module
module load plink2

VCF_DIR="/filer-5/user/weisea/project/vcf_new_miss"
OUT_DIR="/filer-5/user/weisea/project/pca_new"

for vcf_file in "$VCF_DIR"/*.vcf; do
    base_name=$(basename "$vcf_file" .vcf)
    plink2 --vcf "$vcf_file" --threads 8 --pca 20 meanimpute --allow-extra-chr --double-id --out "${OUT_DIR}/${base_name}"
done

