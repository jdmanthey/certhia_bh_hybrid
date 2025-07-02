interactive -p nocona -c 4 -m 8G

source activate bcftools

# define main working directory
workdir=/lustre/scratch/jmanthey/05_certhia_hybrids

############################################
############################################
# prep
############################################
############################################

cd $workdir/04_vcf

# bgzip chromosome of interest
bgzip -c ${workdir}/04_vcf/Ca_0001_Tg_2.vcf > ${workdir}/04_vcf/Ca_0001_Tg_2.vcf.gz

bgzip -c ${workdir}/04_vcf/Ca_0004_Tg_Z.vcf > ${workdir}/04_vcf/Ca_0004_Tg_Z.vcf.gz

#tabix both files
tabix ${workdir}/04_vcf/Ca_0001_Tg_2.vcf.gz

tabix ${workdir}/04_vcf/Ca_0004_Tg_Z.vcf.gz

# pull out two regions of interest
# add header
gunzip -cd ${workdir}/04_vcf/Ca_0001_Tg_2.vcf.gz | head -n100000 | grep "^#" > chr2_subset.vcf
gunzip -cd ${workdir}/04_vcf/Ca_0004_Tg_Z.vcf.gz | head -n100000 | grep "^#" > chrZ_subset.vcf
# add region
tabix ${workdir}/04_vcf/Ca_0001_Tg_2.vcf.gz Ca_0001_Tg_2:120000001-150000000 >> chr2_subset.vcf
tabix ${workdir}/04_vcf/Ca_0004_Tg_Z.vcf.gz Ca_0004_Tg_Z:12000001-38000000 >> chrZ_subset.vcf

############################################
############################################
# filter and run lddecay for each set of individuals
############################################
############################################

cd $workdir/10b_lddecay

# filter for ld in admixed individuals
vcftools --vcf ${workdir}/04_vcf/chr2_subset.vcf --keep mixed.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10b_lddecay/chr2

vcftools --vcf ${workdir}/04_vcf/chrZ_subset.vcf --keep mixed.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10b_lddecay/chrZ

# gzip 
gzip ${workdir}/10b_lddecay/chr2.recode.vcf

gzip ${workdir}/10b_lddecay/chrZ.recode.vcf

# run PopLDdecay
~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10b_lddecay/chr2.recode.vcf.gz \
-MaxDist 200 -OutStat chr2__admixed

~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10b_lddecay/chrZ.recode.vcf.gz \
-MaxDist 200 -OutStat chrZ__admixed

# unzip
gunzip chr2__admixed.stat.gz
gunzip chrZ__admixed.stat.gz

# remove the vcfs
rm ${workdir}/10b_lddecay/chr2.recode.vcf.gz
rm ${workdir}/10b_lddecay/chrZ.recode.vcf.gz



# filter for ld in west individuals
vcftools --vcf ${workdir}/04_vcf/chr2_subset.vcf --keep west.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10b_lddecay/chr2

vcftools --vcf ${workdir}/04_vcf/chrZ_subset.vcf --keep west.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10b_lddecay/chrZ

# gzip 
gzip ${workdir}/10b_lddecay/chr2.recode.vcf

gzip ${workdir}/10b_lddecay/chrZ.recode.vcf

# run PopLDdecay
~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10b_lddecay/chr2.recode.vcf.gz \
-MaxDist 200 -OutStat chr2__west

~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10b_lddecay/chrZ.recode.vcf.gz \
-MaxDist 200 -OutStat chrZ__west

# unzip
gunzip chr2__west.stat.gz
gunzip chrZ__west.stat.gz

# remove the vcfs
rm ${workdir}/10b_lddecay/chr2.recode.vcf.gz
rm ${workdir}/10b_lddecay/chrZ.recode.vcf.gz



# filter for ld in east individuals
vcftools --vcf ${workdir}/04_vcf/chr2_subset.vcf --keep east.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10b_lddecay/chr2

vcftools --vcf ${workdir}/04_vcf/chrZ_subset.vcf --keep east.txt --max-missing 1.0 \
--mac 2 --max-maf 0.49 --recode --recode-INFO-all --out ${workdir}/10b_lddecay/chrZ

# gzip 
gzip ${workdir}/10b_lddecay/chr2.recode.vcf

gzip ${workdir}/10b_lddecay/chrZ.recode.vcf

# run PopLDdecay
~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10b_lddecay/chr2.recode.vcf.gz \
-MaxDist 200 -OutStat chr2__east

~/PopLDdecay/bin/PopLDdecay -InVCF ${workdir}/10b_lddecay/chrZ.recode.vcf.gz \
-MaxDist 200 -OutStat chrZ__east

# unzip
gunzip chr2__east.stat.gz
gunzip chrZ__east.stat.gz

# remove the vcfs
rm ${workdir}/10b_lddecay/chr2.recode.vcf.gz
rm ${workdir}/10b_lddecay/chrZ.recode.vcf.gz

