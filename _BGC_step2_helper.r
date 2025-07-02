options(scipen=999)
library(bgchm)

# read in SNP files
bgc_east <- read.table("bgc_east.txt", sep="\t", header=T, row.names=1)
bgc_west <- read.table("bgc_west.txt", sep="\t", header=T, row.names=1)
bgc_others <- read.table("bgc_others.txt", sep="\t", header=T, row.names=1)

# read in pre-computed hybrid indices
hybrid_indices <- read.table("hybrid_indices.txt", header=T)

# point estimates of cline SDs inferred from the subset of loci
sdc <- 0.2929867
sdv <- 0.261480

# read in start and end number SNPs 
my_args <- commandArgs(trailingOnly=T)
start_number <- my_args[1]
end_number <- my_args[2]
job_number <- my_args[3]

# subset SNP files
bgc_east <- bgc_east[,start_number:end_number]
bgc_west <- bgc_west[,start_number:end_number]
bgc_others <- bgc_others[,start_number:end_number]

## number of loci and hybrids for this analysis
number_loci <- ncol(bgc_others)
number_mixed <- nrow(hybrid_indices)

# create empty matrices for output
gradients <- matrix(NA, nrow=number_loci, ncol=3)
centers <- matrix(NA, nrow=number_loci, ncol=3)

# loop over all loci, fitting genomic clines for each SNP at a time
# using the pre-computed estimated of sdc and sdv
for(a in 1:number_loci) {
	a_rep <- est_genocl(Gx=bgc_others[,a], G0=bgc_east[,a], G1=bgc_west[,a], H = hybrid_indices[,1],
						ploidy="diploid", hier=F, SDc=sdc, SDv=sdv, n_iters=10000)
	gradients[a, ] <- a_rep$gradient
	centers[a, ] <- a_rep$center
}

# write tables of the gradients and centers
write.table(gradients, file=paste0("gradients_", job_number, ".txt"), sep="\t", row.names=F, quote=F)
write.table(centers, file=paste0("centers_", job_number, ".txt"), sep="\t", row.names=F, quote=F)

