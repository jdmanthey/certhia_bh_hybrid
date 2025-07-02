# define main working directory
workdir=/lustre/scratch/jmanthey/05_certhia_hybrids

# make output directories
cd ${workdir}

mkdir 00_fastq
mkdir 01_cleaned
mkdir 01_bam_files
mkdir 02_vcf
mkdir 03_vcf
mkdir 04_vcf
mkdir 05_structure
mkdir 06_diagnostic
mkdir 07_window_stats
mkdir 07_window_stats/windows
mkdir 08_bgchm
mkdir 09_admixture_windows
mkdir 09_admixture_windows/windows
mkdir 10_lddecay
mkdir 10b_lddecay
mkdir 11_splitstree
mkdir 20_filter
mkdir 21_filter_vcf
mkdir 22_certhia_stats
mkdir 23_bgc_script
mkdir 24_ld_script
