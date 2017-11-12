loglikLOOCVcontourVARX1 <- function(lambdaAgrid, 
                                    lambdaBgrid, 
                                    lambdaPgrid, 
                                    Y, 
                                    X, 
                                    lagX=0, 
                                    figure=TRUE, 
                                    verbose=TRUE, 
                                    ...){                                                                        

	# input checks
	if (as.character(class(Y)) != "array"){ 
		stop("Input (Y) is of wrong class.") 
	}
	if (length(dim(Y)) != 3){ 
		stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") 
	}
	if (as.character(class(X)) != "array"){ 
		stop("Input (X) is of wrong class.") 
	}
	if (length(dim(X)) != 3){ 
		stop("Input (X) is of wrong dimensions: either covariate, time or sample dimension is missing.") 
	}
	if (any(dim(Y)[2:3] != dim(X)[2:3])){ 
		stop("Input (X) do not have same dimensions as Y.") 
	}
	if (as.character(class(lambdaAgrid)) != "numeric"){ 
		stop("Input (lambdaAgrid) is of wrong class.") 
	}
	if (as.character(class(lambdaBgrid)) != "numeric"){ 
		stop("Input (lambdaBgrid) is of wrong class.") 
	}
	if (as.character(class(lambdaPgrid)) != "numeric"){ 
		stop("Input (lambdaPgrid) is of wrong class.") 
	}
	if (length(lambdaAgrid) < 1){ 
		stop("Input (lambdaAgrid) is of wrong length.") 
	}
	if (length(lambdaBgrid) < 1){ 
		stop("Input (lambdaBgrid) is of wrong length.") 
	}
	if (length(lambdaPgrid) < 1){ 
		stop("Input (lambdaPgrid) is of wrong length.") 
	}
	if (all(sort(c(length(lambdaAgrid), length(lambdaBgrid), length(lambdaPgrid)))[1:2] == 1)){ 
		stop("Input combination (lambdaA/B/Pgrid) does not forms a two-dimension grid.") 
	}
	if (length(lambdaAgrid) != length(unique(lambdaAgrid))){ 
		stop("Input (lambdaAgrid) contains non-unique values.") 
	}
	if (length(lambdaBgrid) != length(unique(lambdaBgrid))){ 
		stop("Input (lambdaBgrid) contains non-unique values.") 
	}
	if (length(lambdaPgrid) != length(unique(lambdaPgrid))){ 
		stop("Input (lambdaPgrid) contains non-unique values.") 
	}
	if (any(is.na(lambdaAgrid))){ 
		stop("Input (lambdaAgrid) is not a vector of non-negative numbers.") 
	}
	if (any(is.na(lambdaBgrid))){ 
		stop("Input (lambdaBgrid) is not a vector of non-negative numbers.") 
	}
	if (any(is.na(lambdaPgrid))){ 
		stop("Input (lambdaPgrid) is not a vector of non-negative numbers.") 
	}
	if (any(lambdaAgrid <= 0)){ 
		stop("Input (lambdaAgrid) is not a vector of non-negative numbers.") 
	}
	if (any(lambdaBgrid <= 0)){ 
		stop("Input (lambdaBgrid) is not a vector of non-negative numbers.") 
	}
	if (any(lambdaPgrid <= 0)){ 
		stop("Input (lambdaPgrid) is not a vector of non-negative numbers.") 
	}
	if (as.character(class(figure)) != "logical"){ 
		stop("Input (figure) is of wrong class.") 
	}
	if (as.character(class(verbose)) != "logical"){ 
		stop("Input (verbose) is of wrong class.") 
	}


	# case where a grid for both autoregression 
	# penalty parameters is specified
	if ((length(lambdaAgrid) > 1) && 
	    (length(lambdaBgrid) > 1) && 
	    (length(lambdaPgrid) > 1)){
		lambdaPgrid <- lambdaPgrid[1]
	}
	if ((length(lambdaAgrid) > 1) && 
	    (length(lambdaBgrid) > 1) && 
	    (length(lambdaPgrid) == 1)){
		lambdaAgrid <- sort(lambdaAgrid)
		lambdaBgrid <- sort(lambdaBgrid)
		llLOOCV <- matrix(NA, nrow=length(lambdaAgrid), 
		                      ncol=length(lambdaBgrid))
		if (verbose) {
			cat("grid point:", "\n")
		}
		for (kA in 1:length(lambdaAgrid)) {
			for (kB in 1:length(lambdaBgrid)) {
				if (verbose) {
					cat(rep("\b", 100), sep = "")
					cat(paste("lambdaA=", 
				                  lambdaAgrid[kA], 
					          "; lambdaB=",
					          lambdaBgrid[kB], 
					          sep = ""))
				}
				llLOOCV[kA, kB] <- loglikLOOCVVARX1(c(lambdaAgrid[kA],
				                                      lambdaBgrid[kB],
				                                      lambdaPgrid[1]), 
				                                      Y, 
				                                      X, 
				                                      lagX=lagX, ...)
			}
		}
		if (figure) {
			contour(lambdaAgrid, 
			        lambdaBgrid, 
			        -llLOOCV, 
				xlab="lambdaA", 
			        ylab="lambdaB", 
			        main="cross-validated log-likelihood")
		}
		return(list(lambdaA=lambdaAgrid, 
		            lambdaB=lambdaBgrid, 
		            lambdaP=lambdaPgrid[1], 
		            llLOOCV=-llLOOCV))
	}

	# case where a grid for endogeneous autoregression 
	# and precision penalty parameters is specified	
	if ((length(lambdaAgrid) > 1) && 
	    (length(lambdaBgrid) == 1) && 
	    (length(lambdaPgrid) > 1)){
		lambdaAgrid <- sort(lambdaAgrid)
		lambdaPgrid <- sort(lambdaPgrid)
		llLOOCV <- matrix(NA, nrow=length(lambdaAgrid), 
		                      ncol=length(lambdaPgrid))
		if (verbose) {
			cat("grid point:", "\n")
		}
		for (kA in 1:length(lambdaAgrid)) {
			for (kP in 1:length(lambdaPgrid)) {
				if (verbose) {
					cat(rep("\b", 100), sep = "")
					cat(paste("lambdaA=", 
					          lambdaAgrid[kA], 
					          "; lambdaP=",
					          lambdaPgrid[kP],
					          sep = ""))
				}
				llLOOCV[kA, kP] <- loglikLOOCVVARX1(c(lambdaAgrid[kA],
				                                      lambdaBgrid[1],
				                                      lambdaPgrid[kP]), 
				                                      Y, 
				                                      X, 
				                                      lagX=lagX, ...)
			}
		}
		if (figure) {
			contour(lambdaAgrid, 
			        lambdaPgrid, 
			        -llLOOCV, 
			        xlab="lambdaA", 
			        ylab="lambdaP", 
			        main="cross-validated log-likelihood")
		}
		return(list(lambdaA=lambdaAgrid, 
		            lambdaB=lambdaBgrid[1], 
		            lambdaP=lambdaPgrid, 
		            llLOOCV=-llLOOCV))
	}

	# case where a grid for exogeneous autoregression 
	# and precision penalty parameter is specified	
	if ((length(lambdaAgrid) == 1) && 
	    (length(lambdaBgrid) > 1) && 
	    (length(lambdaPgrid) > 1)){
		lambdaBgrid <- sort(lambdaBgrid)
		lambdaPgrid <- sort(lambdaPgrid)
		llLOOCV <- matrix(NA, nrow=length(lambdaBgrid), 
		                      ncol=length(lambdaPgrid))
		if (verbose) {
			cat("grid point:", "\n")
		}
		for (kB in 1:length(lambdaBgrid)) {
			for (kP in 1:length(lambdaPgrid)) {
				if (verbose) {
					cat(rep("\b", 100), sep = "")
					cat(paste("lambdaB=", 
					          lambdaBgrid[kB], 
					          "; lambdaP=",
					          lambdaPgrid[kP], 
					          sep = ""))
				}
				llLOOCV[kB, kP] <- loglikLOOCVVARX1(c(lambdaAgrid[1],
				                                      lambdaBgrid[kB],
				                                      lambdaPgrid[kP]), 
				                                      Y, 
				                                      X, 
				                                      lagX=lagX, ...)
			}
		}
		if (figure) {
			contour(lambdaBgrid, 
			        lambdaPgrid, 
			        -llLOOCV, 
			        xlab="lambdaB", 
			        ylab="lambdaP", 
			        main="cross-validated log-likelihood")
		}
		return(list(lambdaA=lambdaAgrid[1], 
		            lambdaB=lambdaBgrid, 
		            lambdaP=lambdaPgrid, 
		            llLOOCV=-llLOOCV))
	}




}
