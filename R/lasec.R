### LASEC: LANDMARK SAMPLING EVALUATION CURVE ###

# The function requires "geomorph" and "vegan" packages.
# To use LaSEC, run this entire script, then use function "lasec" with appropriate parameters: e.g., > lasec(coord.data=example.txt, n.dim=3, iter=1000)

lasec <- function(coord.data, n.dim, iter=1000, show.progress=T) {

# coord.data is a n x p matrix with n=specimens, p=shape variables (could be unaligned).
# n.dim is spatial dimensionality of coord.data (i.e., 2 or 3 for 2-D, 3-D)
# show.progress is for displaying the nth iteration on console.

	# LIBRARY #
	require(geomorph)
	require(vegan)
	
	# INITIALIZATION #
	matrix.pss.cs <- NULL
	matrix.pss <- NULL
	matrix.chosen <- NULL
	lm.maxpss <- NULL
	lm.minpss <- NULL
	output <- NULL
	if(show.progress==T) {
	  pb <- txtProgressBar(min=0, max=iter, char=" >", style=3)
	  }
	
	# DATA #
	n.lm <- ncol(coord.data)/n.dim   # no. landmarks in data
	list.lm <- seq(1:n.lm)   # make a list of landmarks
	coord.data <- arrayspecs(coord.data, n.lm, n.dim)
	gpa <- gpagen(coord.data, print.progress=FALSE)
	shape.data <- two.d.array(gpa$coords)
	cs.data <- gpa$Csize
	coord.data <- two.d.array(coord.data)
	
for(i in 1:iter) {
  if(show.progress==T) {
    setTxtProgressBar(pb, i)
    }
  new.lm.order <- sample(n.lm, replace=F)   # new sequence of landmarks
	list.pss <- NULL
	list.pss.cs <- NULL
	for(j in 3:n.lm) {
		subsample.lm <- new.lm.order[1:j]
		if(n.dim==3) {
			chosen.nvar <- sort(c(3*subsample.lm-2, 3*subsample.lm-1, 3*subsample.lm))
		}
		if(n.dim==2) {
			chosen.nvar <- sort(c(2*subsample.lm-1, 2*subsample.lm))
		}
		sampled.coord <- coord.data[,chosen.nvar]
		sampled.coord <- arrayspecs(sampled.coord, j, n.dim)
		sampled.gpa <- gpagen(sampled.coord, print.progress=FALSE)
		sampled.cs <- sampled.gpa$Csize
		sampled.shape <- two.d.array(sampled.gpa$coords)
		
		# CALCULATE 'FIT' BETWEEN FULL & SUBSAMPLED DATA #
		pss <- protest(shape.data, sampled.shape, permutations=0)$ss
		pss.cs <- protest(cs.data, sampled.cs, permutations=0)$ss
			
		# RECORD 'FIT' BETWEEN FULL & SUBSAMPLED DATA #
		list.pss <- append(list.pss, pss)
		list.pss.cs <- append(list.pss.cs, pss.cs)
	}
	matrix.chosen <- append(matrix.chosen, new.lm.order)
	matrix.pss <- append(matrix.pss, c(NA, NA, list.pss))
	matrix.pss.cs <- append(matrix.pss.cs, c(NA, NA, list.pss.cs))	
}

matrix.chosen <- matrix(matrix.chosen, nrow=iter, byrow=T)
matrix.fit <- 1-matrix(matrix.pss, nrow=iter, byrow=T)
median.fit <- apply(matrix.pss, 2, median)
matrix.fit.cs <- 1-matrix(matrix.pss.cs, nrow=iter, byrow=T)
median.fit.cs <- apply(matrix.pss.cs, 2, median)

# PLOT PROC SS CURVE #
pdf("LaSEC_SamplingCurve_Shape.pdf")
plot(x=NA, xlim=c(0, n.lm), ylim=c(0,1), xlab="No. landmarks", ylab="Fit")
for(i in 1:iter) {  
	par(new=T)
	plot(matrix.fit[i,], xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="grey", xlab="", ylab="", axes=F)
}
par(new=T)
plot(median.fit, xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="black", lwd=3, xlab="", ylab="", axes=F)
dev.off()

# PLOT PROC SS CENTROID SIZE CURVE #
pdf("LaSEC_SamplingCurve_CS.pdf")
plot(x=NA, xlim=c(0, n.lm), ylim=c(0,1), xlab="No. landmarks", ylab="Fit")
for(i in 1:iter) {  
	par(new=T)
	plot(matrix.fit.cs[i,], xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="grey", xlab="", ylab="", axes=F)
}
par(new=T)
plot(median.fit.cs, xlim=c(0, n.lm), ylim=c(0, 1), type="l", col="black", lwd=3, xlab="", ylab="", axes=F)
dev.off()

# OUTPUT #
output$fit <- matrix.fit
output$median.fit <- median.fit
output$maxpss.landmark <- table(lm.maxpss)
output$minpss.landmark <- table(lm.minpss)
output$pss.cs <- median.fit.cs
return(output)		
}
