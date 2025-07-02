options(scipen=999)
library(bgchm)

####################################
####################################
# functions used in below sections
####################################
####################################

# read in a small vcf (don't use for large vcf files)
read_vcf <- function(input_file) {
  header <- readLines(input_file)
  header <- header[grep('^#C', header)]
  header <- strsplit(header, "\t")[[1]]
  vcf <- read.table(input_file, header=F)
  colnames(vcf) <- header
  return(vcf)
}

####################################
####################################
# find diagnostic snps
####################################
####################################

# read in vcf 
vcf <- read_vcf("diagnostic.vcf")

# define groups
east_inds <- c("Certhia_americana_NCarolina_JK05-134", "Certhia_americana_NCarolina_JK05-135", "Certhia_americana_NCarolina_JK05-138", "Certhia_americana_WVirginia_UWBM107061", "Certhia_americana_WVirginia_UWBM112054", "Certhia_americana_WVirginia_UWBM112055")
west_inds <- c("Certhia_americana_Colorado_GMS1712", "Certhia_americana_Colorado_GMS1718", "Certhia_americana_Colorado_JK05-305", "Certhia_americana_Utah_UWBM113162", "Certhia_americana_Utah_UWBM113167", "Certhia_americana_Utah_UWBM113168")
others_inds <- colnames(vcf[10:ncol(vcf)])[colnames(vcf[10:ncol(vcf)]) %in% c(east_inds, west_inds) == F]

# remove all non-genotype data
for(a in 10:ncol(vcf)) {
	vcf[,a] <- substr(vcf[,a], 1, 3)
}
# convert files for format needed in bgchm
for(a in 10:ncol(vcf)) {
	vcf[vcf[,a] == "0/0", a] <- "0"
	vcf[vcf[,a] == "0/1", a] <- "1"
	vcf[vcf[,a] == "1/1", a] <- "2"
	vcf[vcf[,a] == "./.", a] <- "NA"	
	vcf[,a] <- as.numeric(vcf[,a])
}

# subset vcf into separate objects
vcf_info <- vcf[,1:5]
others_bgc <- vcf[,colnames(vcf) %in% others_inds]
east_bgc <- vcf[,colnames(vcf) %in% east_inds]
west_bgc <- vcf[,colnames(vcf) %in% west_inds]

# convert to matrix
others_bgc <- as.matrix(t(others_bgc))
east_bgc <- as.matrix(t(east_bgc))
west_bgc <- as.matrix(t(west_bgc))

# sample random 2000 loci for estimating cline SDs and hybrid indices
rsamp <- sort(sample(1:ncol(others_bgc), 2000))
others_bgc_samp <- others_bgc[,rsamp]
east_bgc_samp <- east_bgc[,rsamp]
west_bgc_samp <- west_bgc[,rsamp]

# measure parental allele frequencies
parental_freqs <- est_p(G0=east_bgc_samp, G1=west_bgc_samp, model="genotype", ploidy="diploid", HMC=FALSE)

# estimate hybrid indices using default HMC settings
# and using point estimates of allele frequencies
hybrid_indices <- est_hi(Gx=others_bgc_samp, p0=parental_freqs$p0[,1], p1=parental_freqs$p1[,1], model="genotype",ploidy="diploid")

# check sampling of HMC (check the n_eff (> n_chains * 100), and Rhat (less than 1.1, close to 1))
hybrid_indices

# write to output
write.table(file="hybrid_indices.txt", hybrid_indices$hi,row.names=FALSE,quote=FALSE)

# estimate genomic clines for subset of loci with input from hybrid indices and parental allele freqs.
# using 10000 iterations -> used to estimate cline SDs
genomic_clines <- est_genocl(Gx=others_bgc_samp, p0=parental_freqs$p0[,1], p1=parental_freqs$p1[,1], H=hybrid_indices$hi[,1], model="genotype", ploidy="diploid", hier=TRUE, n_iters=10000, n_thin=10)

# get the SDs (used for step 2)
genomic_clines$SDc
genomic_clines$SDv

# write bgc files for next step
write.table(others_bgc, file="bgc_others.txt", sep="\t", quote=F)
write.table(east_bgc, file="bgc_east.txt", sep="\t", quote=F)
write.table(west_bgc, file="bgc_west.txt", sep="\t", quote=F)
write.table(vcf_info, file="vcf_info.txt", sep="\t", quote=F, row.names=F)






