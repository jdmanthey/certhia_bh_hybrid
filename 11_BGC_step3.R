options(scipen=999)
library(bgchm)

locus_info <- read.table("vcf_info.txt", header=F)

# output from step 2 (centers and gradients files) in a directory called step2_output
centers <- list()
gradients <- list()
for(a in 1:100) {
	centers[[a]] <- read.table(paste0("step2_output/centers_", a, ".txt"), header=T)
	gradients[[a]] <- read.table(paste0("step2_output/gradients_", a, ".txt"), header=T)
}
centers <- do.call(rbind, centers)
gradients <- do.call(rbind, gradients)

# combine vcf info and centers and gradients
x <- data.frame(chrom=as.character(locus_info[,1]), site=as.numeric(locus_info[,2]), center_estimate=as.numeric(centers[,1]), center_low=as.numeric(centers[,2]), center_high=as.numeric(centers[,3]), gradient_estimate=as.numeric(gradients[,1]), gradient_low=as.numeric(gradients[,2]), gradient_high=as.numeric(gradients[,3]))

# number of loci per chromosome
table(x$chrom)

# bias of ancestry from east or west
table(x$chrom[x$center_low > 0.5])
table(x$chrom[x$center_high < 0.5])

# selection against introgression
table(x$chrom[x$gradient_low > 1])

# selection for introgression
table(x$chrom[x$gradient_high < 1])



# check clustering of outliers

sites <- seq(from=1,to=nrow(x))
x2 <- cbind(x, sites)


plot(0, 0, col="white", xlab="", ylab="", xlim=c(min(x2$sites), max(x2$sites)), ylim=c(min(x2$gradient_low), max(x2$gradient_high)))
segments(x2$sites, x2$gradient_low, x2$sites, x2$gradient_high, lwd=0.1, col="gray", lty=1)
abline(1,0)
points(x2$sites, x2$gradient_estimate, pch=19, cex=0.1)

plot(0, 0, col="white", xlab="", ylab="", xlim=c(min(x2$sites), max(x2$sites)), ylim=c(min(x2$center_low), max(x2$center_high)))
segments(x2$sites, x2$center_low, x2$sites, x2$center_high, lwd=0.1, col="gray", lty=1)
abline(0.5,0)
points(x2$sites, x2$center_estimate, pch=19, cex=0.1)



n_windows <- floor(nrow(x2) / 10)
sites <- seq(from=1, to=n_windows)
gl_mean <- list()
for(a in 1:n_windows) {
	a_rep <- x2[(a*10 - 9):(a*10),]
	gl_mean[[a]] <- mean(a_rep$gradient_low)
}
gl_mean <- unlist(gl_mean)
plot(sites, gl_mean, pch=19, cex=0.2)
abline(1,0)
plot(x2$sites, x2$gradient_low, pch=19, cex=0.2)
abline(1,0)







# gradients and centers per chromosome
chromosomes <- unique(x$chrom)
gradient_estimate <- c()
gradient_low <- c()
gradient_high <- c()
for(a in 1:length(chromosomes)) {
	a_rep <- x[x$chrom == chromosomes[a],]
	gradient_estimate <- c(gradient_estimate, mean(a_rep$gradient_estimate))
	gradient_low <- c(gradient_low, mean(a_rep$gradient_low))
	gradient_high <- c(gradient_high, mean(a_rep$gradient_high))
}
plot_windows <- seq(from=1, to=length(chromosomes))

plot(0, 0, col="white", xlab="", ylab="", xlim=c(min(plot_windows), max(plot_windows)), ylim=c(min(gradient_low), max(gradient_high)))
segments(plot_windows, gradient_low, plot_windows, gradient_high, lwd=1, col="gray", lty=1)
abline(1,0)
points(plot_windows, gradient_estimate, pch=19, cex=1)



# gradients and centers per 2 Mbp window

# read in index
fai <- read.table("06_certhia_reordered.fasta.fai")
fai <- fai[fai[,2] >= 10000000, ] # keep only chromosomes of at least 10 Mbp 

window <- list()
chromosome <- list()
start <- list()
end <- list()
gradient_low <- list()
gradient_high <- list()
gradient_estimate <- list()
center_estimate <- list()
center_low <- list()
center_high <- list()
n_snps <- list()
window_size <- 2000000
min_snps_per_window <- 5
counter <- 1
for(a in 1:nrow(fai)) {
	# subset to chromosome
	a_rep <- x[x$chrom == fai[a,1],]
	
	# if there are at least 10 SNPs on this chromosome, continue
	if(nrow(a_rep) > 10) {
		# determine how many windows
		a_windows <- floor(fai[a,2] / window_size)
		
		# loop for each window
		for(b in 1:a_windows) {
			b_min <- (window_size * b - (window_size - 1))
			b_max <- (window_size * b)
			b_rep <- a_rep[a_rep$site >= b_min & a_rep$site <= b_max,]
			# calculate means if at least min_snps_per_window
			if(nrow(b_rep) >= min_snps_per_window ) {
				window[[counter]] <- counter
				chromosome[[counter]] <- fai[a,1]
				start[[counter]] <- b_min
				end[[counter]] <- b_max
				gradient_estimate[[counter]] <- mean(b_rep$gradient_estimate)
				gradient_low[[counter]] <- mean(b_rep$gradient_low)
				gradient_high[[counter]] <- mean(b_rep$gradient_high)
				center_estimate[[counter]] <- mean(b_rep$center_estimate)
				center_low[[counter]] <- mean(b_rep$center_low)
				center_high[[counter]] <- mean(b_rep$center_high)
				n_snps[[counter]] <- nrow(b_rep)				
			}
			counter <- counter + 1 # count goes up whether or not there were SNPs
		}
	}	
}
window <- unlist(window)
chromosome <- unlist(chromosome)
start <- unlist(start)
end <- unlist(end)
gradient_low <- unlist(gradient_low)
gradient_high <- unlist(gradient_high)
gradient_estimate <- unlist(gradient_estimate)
center_estimate <- unlist(center_estimate)
center_low <- unlist(center_low)
center_high <- unlist(center_high)
n_snps <- unlist(n_snps)

output <- data.frame(window=as.numeric(window), chromosome=as.character(chromosome), start=as.numeric(start), end=as.numeric(end), gradient_estimate=as.numeric(gradient_estimate), gradient_low=as.numeric(gradient_low), gradient_high=as.numeric(gradient_high), center_estimate=as.numeric(center_estimate), center_low=as.numeric(center_low), center_high=as.numeric(center_high), n_snps=as.numeric(n_snps))

# edit window numbers (with gaps there are issues plotting)
output$window <- seq(from=1, to=nrow(output), by=1)


# initial plots
plot(0, 0, col="white", xlab="", ylab="", xlim=c(min(output$window), max(output$window)), ylim=c(min(output$gradient_low), max(output$gradient_high)))
segments(output$window, output$gradient_low, output$window, output$gradient_high, lwd=0.3, col="gray", lty=1)
abline(1,0)
points(output$window, output$gradient_estimate, pch=19, cex=0.1)

plot(0, 0, col="white", xlab="", ylab="", xlim=c(min(output$window), max(output$window)), ylim=c(min(output$center_low), max(output$center_high)))
segments(output$window, output$center_low, output$window, output$center_high, lwd=0.3, col="gray", lty=1)
abline(0.5,0)
points(output$window, output$center_estimate, pch=19, cex=0.1)


# separate outlier and non-outlier windows
output1a <- output[output$gradient_low <= 1,]
output1b <- output[output$gradient_low > 1,]

output2a <- output[output$center_low <= 0.5 & output$center_high >= 0.5, ]
output2b <- output[output$center_low > 0.5,]
output2c <- output[output$center_high < 0.5,]

# identify distinct regions of the original output
# 1. genome-wide random
# 2. chr 4 outlier region
# 3. chr Z outlier region
# make set 1 not overlap with other regions
set1 <- x[(x$chrom == "Ca_0001_Tg_2" & x$site >= 120000001 & x$site <= 150000000) == F & (x$chrom == "Ca_0004_Tg_Z" & x$site >= 12000001 & x$site <= 38000000) == F,]
set2 <- x[x$chrom == "Ca_0001_Tg_2" & x$site >= 120000001 & x$site <= 150000000,]
set3 <- x[x$chrom == "Ca_0004_Tg_Z" & x$site >= 12000001 & x$site <= 38000000,]

# count outliers in these subsets
nrow(set1)
nrow(set1[set1$gradient_low > 1 | set1$center_low > 0.5 | set1$center_high < 0.5,])
nrow(set2)
nrow(set2[set2$gradient_low > 1 | set2$center_low > 0.5 | set2$center_high < 0.5,])
nrow(set3)
nrow(set3[set3$gradient_low > 1 | set3$center_low > 0.5 | set3$center_high < 0.5,])

########################################################################
# plotting
########################################################################

# what are the unique chromosomes and their bounding areas for plotting?
total_windows <- max(output$window)
chr <- unique(output$chromosome)
chr_polygons <- list()
# make the plotting polygons
for(a in 1:length(chr)) {
	a1 <- output$window[output$chromosome == chr[a]]
	a2 <- max(a1)
	a1 <- a1[1]
	chr_polygons[[a]] <- rbind(c(a1, min(output$gradient_low)), c(a2, min(output$gradient_low)), c(a2, max(output$gradient_high)), c(a1, max(output$gradient_high)), c(a1, min(output$gradient_low)))
}


# set up plotting dimensions
par(mar=c(1,5,1,1))
layout(matrix(c(1,1,1,1,1,1,1,1,
				2,2,2,2,2,2,2,2,
				3,3,4,4,5,5,6,6,
				3,3,4,4,5,5,6,6), 
				4, 8, byrow=T))
par(lwd=1)


# plot windows
plot(0, 0, col="white", xlab="", xlim=c(min(output$window), max(output$window)), ylim=c(min(output$gradient_low), max(output$gradient_high)), xaxt="n",bty="n", ylab="Gradient Estimate")
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}

# plot non outliers
segments(output1a$window, output1a$gradient_low, output1a$window, output1a$gradient_high, lwd=0.6, col="gray", lty=1)
abline(1,0, lwd=0.5)
points(output1a$window, output1a$gradient_estimate, pch=19, cex=0.3)

# plot outliers
segments(output1b$window, output1b$gradient_low, output1b$window, output1b$gradient_high, lwd=0.8, col="blue", lty=1)
points(output1b$window, output1b$gradient_estimate, pch=19, cex=0.3, col="darkblue")


# plot windows
plot(0, 0, col="white", xlab="", xlim=c(min(output$window), max(output$window)), ylim=c(min(output$center_low), max(output$center_high)), xaxt="n",bty="n", ylab="Center Estimate")
odd <- 0
for(a in 1:length(chr_polygons)) {
	if(odd == 1) {
		polygon(chr_polygons[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}

# plot non outliers
segments(output2a$window, output2a$center_low, output2a$window, output2a$center_high, lwd=0.6, col="gray", lty=1)
abline(0.5,0, lwd=0.5)
points(output2a$window, output2a$center_estimate, pch=19, cex=0.3)

segments(output2b$window, output2b$center_low, output2b$window, output2b$center_high, lwd=0.8, col="orange", lty=1)
points(output2b$window, output2b$center_estimate, pch=19, cex=0.3, col="darkorange")

segments(output2c$window, output2c$center_low, output2c$window, output2c$center_high, lwd=0.8, col="orange", lty=1)
points(output2c$window, output2c$center_estimate, pch=19, cex=0.3, col="darkorange")


par(mar=c(5,5,3,1))


# plot per window genomic clines

non_outlier <- output[output$gradient_low <= 1 & output$center_low <= 0.5 & output$center_high >= 0.5,]
outlier <- output[output$gradient_low > 1 | output$center_low > 0.5 | output$center_high < 0.5,]
output_plot <- rbind(non_outlier, outlier)
gencline_plot(center=output_plot$center_estimate, v=output_plot$gradient_estimate, cvec=c(rep("gray", nrow(non_outlier)), rep("black", nrow(outlier))), pdf=FALSE)


# plot subset genomic clines (random 200 loci each)
set1_subset <- set1[sample(1:nrow(set1), 200),]
set2_subset <- set2[sample(1:nrow(set2), 200),]
set3_subset <- set3[sample(1:nrow(set3), 200),]

# pull out outliers and non-outliers and combine
set1_subset_non_outlier <- set1_subset[set1_subset$gradient_low <= 1 & set1_subset$center_low <= 0.5 & set1_subset$center_high >= 0.5,]
set1_subset_outlier <- set1_subset[set1_subset$gradient_low > 1 | set1_subset$center_low > 0.5 | set1_subset$center_high < 0.5,]
set1_subset <- rbind(set1_subset_non_outlier, set1_subset_outlier)

set2_subset_non_outlier <- set2_subset[set2_subset$gradient_low <= 1 & set2_subset$center_low <= 0.5 & set2_subset$center_high >= 0.5,]
set2_subset_outlier <- set2_subset[set2_subset$gradient_low > 1 | set2_subset$center_low > 0.5 | set2_subset$center_high < 0.5,]
set2_subset <- rbind(set2_subset_non_outlier, set2_subset_outlier)

set3_subset_non_outlier <- set3_subset[set3_subset$gradient_low <= 1 & set3_subset$center_low <= 0.5 & set3_subset$center_high >= 0.5,]
set3_subset_outlier <- set3_subset[set3_subset$gradient_low > 1 | set3_subset$center_low > 0.5 | set3_subset$center_high < 0.5,]
set3_subset <- rbind(set3_subset_non_outlier, set3_subset_outlier)


# plot
gencline_plot(center=set1_subset$center_estimate, v=set1_subset$gradient_estimate, cvec=c(rep("gray", nrow(set1_subset_non_outlier)), rep("black", nrow(set1_subset_outlier))), pdf=FALSE)

gencline_plot(center=set2_subset$center_estimate, v=set2_subset$gradient_estimate, cvec=c(rep("gray", nrow(set2_subset_non_outlier)), rep("black", nrow(set2_subset_outlier))), pdf=FALSE)

gencline_plot(center=set3_subset$center_estimate, v=set3_subset$gradient_estimate, cvec=c(rep("gray", nrow(set3_subset_non_outlier)), rep("black", nrow(set3_subset_outlier))), pdf=FALSE)











