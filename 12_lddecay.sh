#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=lddecay
#SBATCH --nodes=1 --ntasks=2
#SBATCH --partition nocona
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-26

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/05_certhia_hybrids

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# filter for ld in admixed individuals
vcftools --vcf ${workdir}/04_vcf/${input_array}.vcf --keep mixed.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10_lddecay/${input_array}_mixed

# gzip
gzip ${workdir}/10_lddecay/${input_array}_mixed.recode.vcf

# run PopLDdecay
~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10_lddecay/${input_array}_mixed.recode.vcf.gz \
-MaxDist 200 -OutStat ${input_array}__admixed

# unzip
gunzip ${input_array}__admixed.stat.gz

# filter for ld in west individuals
vcftools --vcf ${workdir}/04_vcf/${input_array}.vcf --keep west.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10_lddecay/${input_array}_west

# gzip
gzip ${workdir}/10_lddecay/${input_array}_west.recode.vcf

# run PopLDdecay
~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10_lddecay/${input_array}_west.recode.vcf.gz \
-MaxDist 200 -OutStat ${input_array}__west

# unzip
gunzip ${input_array}__west.stat.gz

# filter for ld in east individuals
vcftools --vcf ${workdir}/04_vcf/${input_array}.vcf --keep east.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10_lddecay/${input_array}_east

# gzip
gzip ${workdir}/10_lddecay/${input_array}_east.recode.vcf

# run PopLDdecay
~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10_lddecay/${input_array}_east.recode.vcf.gz \
-MaxDist 200 -OutStat ${input_array}__east

# unzip
gunzip ${input_array}__east.stat.gz







