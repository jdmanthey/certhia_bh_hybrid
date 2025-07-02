options(scipen=999)

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

x_files <- list.files(pattern="*recode.vcf")
outnames <- paste0(sapply(strsplit(x_files, ".recode"), "[[", 1), ".diagnostic.vcf")

# define groups
east_inds <- c("Certhia_americana_NCarolina_JK05-134", "Certhia_americana_NCarolina_JK05-135", "Certhia_americana_NCarolina_JK05-138", "Certhia_americana_WVirginia_UWBM107061", "Certhia_americana_WVirginia_UWBM112054", "Certhia_americana_WVirginia_UWBM112055")
west_inds <- c("Certhia_americana_Colorado_GMS1712", "Certhia_americana_Colorado_GMS1718", "Certhia_americana_Colorado_JK05-305", "Certhia_americana_Utah_UWBM113162", "Certhia_americana_Utah_UWBM113167", "Certhia_americana_Utah_UWBM113168")

# loop for each vcf
for(a in 1:length(x_files)) {
	print(a)
	# read in vcf 
	vcf <- read_vcf(x_files[a])
	
	# tracker of what sites to keep
	keep_vcf <- seq(from=1, to=nrow(vcf), by=1)
	
	# subset 
	east <- vcf[,colnames(vcf) %in% east_inds]
	west <- vcf[,colnames(vcf) %in% west_inds]
	
	# keep only genotype info
	for(b in 1:ncol(east)) {
		east[,b] <- substr(east[,b], 1, 3)
	}
	for(b in 1:ncol(west)) {
		west[,b] <- substr(west[,b], 1, 3)
	}
	
	# find sites monomorphic in west
	keep <- list()
	for(b in 1:nrow(west)) {
		b_rep <- as.character(west[b,])
		b_rep <- b_rep[b_rep != "./."]
		# needs a minimum length
		if(length(b_rep) >= 4) {
			b_rep <- unique(b_rep)
			if(length(b_rep) == 1 & b_rep[1] != "0/1") { keep[[b]] <- b } 
		}
	}
	keep <- unlist(keep)
	# subset to only those sites identified
	keep_vcf <- keep_vcf[keep]
	east <- east[keep,]
	west <- west[keep,]
	
	# find sites monomorphic in east
	keep <- list()
	for(b in 1:nrow(east)) {
		b_rep <- as.character(east[b,])
		b_rep <- b_rep[b_rep != "./."]
		# needs a minimum length
		if(length(b_rep) >= 4) {
			b_rep <- unique(b_rep)
			if(length(b_rep) == 1 & b_rep[1] != "0/1") { keep[[b]] <- b } 
		}
	}
	keep <- unlist(keep)
	# subset to only those sites identified
	keep_vcf <- keep_vcf[keep]
	east <- east[keep,]
	west <- west[keep,]

	
	# find putatively fixed differences
	keep <- list()
	for(b in 1:nrow(east)) {
		b_rep1 <- unique(as.character(east[b,]))
		b_rep2 <- unique(as.character(west[b,]))
		b_rep1 <- b_rep1[b_rep1 != "./."]
		b_rep2 <- b_rep2[b_rep2 != "./."]
		if(b_rep1 != b_rep2) { keep[[b]] <- b } 
	}
	keep <- unlist(keep)
	# subset to only those sites identified
	keep_vcf <- keep_vcf[keep]
	east <- east[keep,]
	west <- west[keep,]

	# subset full vcf
	vcf <- vcf[keep_vcf,]
	
	# remove sites with missing data (if desired or comment out this loop)
	keep <- list()
	for(b in 1:nrow(vcf)) {
		b_rep <- substr(vcf[b,10:ncol(vcf)], 1, 3)
		if(length(b_rep) == length(b_rep[b_rep != "./."])) {
			keep[b] <- b
		}
	}
	keep <- unlist(keep)
	vcf <- vcf[keep,]
	
	# write table
	if(nrow(vcf) > 0) {
		write.table(vcf, file=outnames[a], row.names=F, quote=F, sep="\t")
	}
}


x_files <- list.files(pattern="*diagnostic.vcf")


# how many with and without missing data?
total <- 0
total2 <- 0
for(a in 1:length(x_files)) {
	print(a)
	a_rep <- read_vcf(x_files[a])
	a_rep2 <- a_rep[,10:ncol(a_rep)]
	n_inds <- ncol(a_rep2)
	for(b in 1:ncol(a_rep2)) {
		a_rep2[,b] <- substr(a_rep2[,b], 1, 3)
	}
	keep <- list()
	for(b in 1:nrow(a_rep2)) {
		b_rep <- a_rep2[b,]
		if(length(b_rep[b_rep != "./."]) == n_inds) {
			keep[[b]] <- b
		}
	}
	keep <- unlist(keep)
	a_rep2 <- a_rep2[keep,]
	total <- total + nrow(a_rep)
	total2 <- total2 + nrow(a_rep2)
}
# with some missing
total
# with no missing
total2


# write all to single file
for(a in 1:length(x_files)) {
	a_rep <- read_vcf(x_files[a])
	if(a == 1) {
		write.table(a_rep, file="diagnostic.vcf", row.names=F, quote=F, sep="\t")
	} else {
		write.table(a_rep, file="diagnostic.vcf", row.names=F, col.names=F, quote=F, sep="\t", append=T)
	}
}


