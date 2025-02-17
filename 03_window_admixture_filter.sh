#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition nocona
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-26

source activate bcftools

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} scaffolds.txt | tail -n1 )

# define main working directory
workdir=/lustre/scratch/jmanthey/05_certhia_hybrids


# filter for windowed admixture
# up to 20% missing
# mac = 2
# no outgroup
vcftools --vcf ${workdir}/04_vcf/${input_array}.vcf --keep ingroup.txt --max-missing 0.8 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/09_admixture_windows/${input_array}

# bgzip and index the vcf files
bgzip -c ${workdir}/09_admixture_windows/${input_array}.recode.vcf > ${workdir}/09_admixture_windows/${input_array}.recode.vcf.gz

tabix -p vcf ${workdir}/09_admixture_windows/${input_array}.recode.vcf.gz

