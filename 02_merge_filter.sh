#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=merge_filter
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-26

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/05_certhia_hybrids

# define variables
region_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/scaffolds.txt | tail -n1 )

# run bcftools to merge the vcf files
bcftools merge -m id --regions ${region_array} ${workdir}/03_vcf/*vcf.gz > \
${workdir}/04_vcf/${region_array}.vcf

# filter for structure (PCA, ADMIXTURE)
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep ingroup.txt \
--max-missing 0.9 --mac 3 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/05_structure/${region_array}

# filter for network (splitstree)
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf \
--max-missing 0.9 --mac 3 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/11_splitstree/${region_array}

# filter for diagnostic site identification
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep parentals.txt \
--max-missing 0.8 --mac 8 --max-alleles 2 --recode \
--recode-INFO-all --out ${workdir}/06_diagnostic/${region_array}

# filter for GADMA of parentals
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep parentals.txt \
--max-missing 1.0 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/12_gadma/${region_array}

# filter for windowed ADMIXTURE
# up to 20% missing, mac = 4, no outgroup
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep ingroup.txt --max-missing 0.8 \
--mac 4 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/09_admixture_windows/${region_array}

# bgzip and index the windowed admixture files
bgzip -c ${workdir}/09_admixture_windows/${region_array}.recode.vcf > \
${workdir}/09_admixture_windows/${region_array}.recode.vcf.gz

tabix -p vcf ${workdir}/09_admixture_windows/${region_array}.recode.vcf.gz

# filter for stats
vcftools --vcf ${workdir}/04_vcf/${region_array}.vcf --keep ingroup.txt \
--max-missing 0.8 --max-alleles 2 --max-maf 0.49 --recode \
--recode-INFO-all --out ${workdir}/07_window_stats/${region_array}

# bgzip and index
bgzip ${workdir}/07_window_stats/${region_array}.recode.vcf

tabix ${workdir}/07_window_stats/${region_array}.recode.vcf.gz



