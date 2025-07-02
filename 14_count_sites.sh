#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=count_sites
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G



for i in $( ls *vcf.gz ); do
	echo $i >> filtering_stats.txt
	
	echo "number sites" >> filtering_stats.txt

	gzip -cd $i | grep -v "^#" | cut -f5 |  wc -l >> filtering_stats.txt
	
	echo "number snps" >> filtering_stats.txt
	
	gzip -cd $i | grep -v "^#" | cut -f5 | grep -v "\\." | wc -l >> filtering_stats.txt
	
done

