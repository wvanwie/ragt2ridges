loglikLOOCVcontourVAR1fused <- function(lambdaAgrid, 
                                        lambdaFgrid, 
                                        Y, 
                                        id, 
                                        lambdaP, 
                                        figure=TRUE, 
                                        verbose = TRUE, 
					...){                                                                  

	if (as.character(class(Y)) != "array"){ stop("Input (Y) is of wrong class.") }
	if (length(dim(Y)) != 3){ stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") }
	if (as.character(class(id)) != "numeric" & as.character(class(id)) != "integer"){ stop("Input (id) is of wrong class.") }
	if (length(id) != dim(Y)[3]){ stop("Input (id) is of wrong length: should equal sample dimension of Y.") }	
	if (as.character(class(lambdaAgrid)) != "numeric"){ stop("Input (lambdaAgrid) is of wrong class.") }
	if (as.character(class(lambdaFgrid)) != "numeric"){ stop("Input (lambdaFgrid) is of wrong class.") }
	if (length(lambdaAgrid) < 2){ stop("Input (lambdaAgrid) is of wrong length.") }
	if (length(lambdaFgrid) < 2){ stop("Input (lambdaFgrid) is of wrong length.") }
	if (any(is.na(lambdaAgrid))){ stop("Input (lambdaAgrid) is not a vector of non-negative numbers.") }
	if (any(is.na(lambdaFgrid))){ stop("Input (lambdaFgrid) is not a vector of non-negative numbers.") }
	if (any(lambdaAgrid <= 0)){ stop("Input (lambdaAgrid) is not a vector of non-negative numbers.") }
	if (any(lambdaFgrid <= 0)){ stop("Input (lambdaFgrid) is not a vector of non-negative numbers.") }
	if (length(lambdaP) != 1){ stop("Input (lambdaP) is of wrong length.") }
	if (is.na(lambdaP)){ stop("Input (lambdaP) is not a vector of non-negative numbers.") }
	if (lambdaP <= 0){ stop("Input (lambdaP) is not a vector of non-negative numbers.") }
	if (as.character(class(figure)) != "logical"){ stop("Input (figure) is of wrong class.") }
	if (as.character(class(verbose)) != "logical"){ stop("Input (verbose) is of wrong class.") }

	lambdaAgrid <- sort(lambdaAgrid)
	lambdaFgrid <- sort(lambdaFgrid)
	llLOOCV <- matrix(NA, nrow = length(lambdaAgrid), ncol = length(lambdaFgrid))
	if (verbose) {
		cat("grid point:", "\n")
	}
	for (kA in 1:length(lambdaAgrid)) {
		for (kF in 1:length(lambdaFgrid)) {
			if (verbose) {
				cat(rep("\b", 100), sep = "")
				cat(paste("lambdaA=", lambdaAgrid[kA], "; lambdaF=",lambdaFgrid[kF], sep = ""))
			}
			llLOOCV[kA, kF] <- loglikLOOCVVAR1fused(c(lambdaAgrid[kA],lambdaFgrid[kF],lambdaP), Y, id, ...)
		}
	}
	if (figure) {
		contour(lambdaAgrid, lambdaFgrid, -llLOOCV, xlab="lambdaA", ylab="lambdaF", main="cross-validated log-likelihood")
	}
	return(list(lambdaA=lambdaAgrid, lambdaF=lambdaFgrid, llLOOCV=-llLOOCV))
}
