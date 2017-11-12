loglikLOOCVcontourVAR2 <- function(lambdaA1grid, 
                                   lambdaA2grid, 
                                   lambdaPgrid, 
                                   Y, 
                                   figure = TRUE, 
                                   verbose = TRUE, 
                                   ...){

	# input checks
	if (as.character(class(Y)) != "array"){ 
		stop("Input (Y) is of wrong class.") 
	}
	if (length(dim(Y)) != 3){ 
		stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") 
	}
	if (as.character(class(lambdaA1grid)) != "numeric"){ 
		stop("Input (lambdaA1grid) is of wrong class.") 
	}
	if (as.character(class(lambdaA2grid)) != "numeric"){ 
		stop("Input (lambdaA2grid) is of wrong class.") 
	}
	if (as.character(class(lambdaPgrid)) != "numeric"){ 
		stop("Input (lambdaPgrid) is of wrong class.") 
	}
	if (length(lambdaA1grid) < 1){ 
		stop("Input (lambdaA1grid) is of wrong length.") 
	}
	if (length(lambdaA2grid) < 1){ 
		stop("Input (lambdaA2grid) is of wrong length.") 
	}
	if (length(lambdaPgrid) < 1){ 
		stop("Input (lambdaPgrid) is of wrong length.") 
	}
	if (all(sort(c(length(lambdaA1grid), length(lambdaA2grid), length(lambdaPgrid)))[1:2] == 1)){ 
		stop("Input combination (lambdaA1/A2/Pgrid) does not form a two-dimension grid.") 
	}
	if (length(lambdaA1grid) != length(unique(lambdaA1grid))){ 
		stop("Input (lambdaA1grid) contains non-unique values.") 
	}
	if (length(lambdaA2grid) != length(unique(lambdaA2grid))){ 
		stop("Input (lambdaA2grid) contains non-unique values.") 
	}
	if (length(lambdaPgrid) != length(unique(lambdaPgrid))){ 
		stop("Input (lambdaPgrid) contains non-unique values.") 
	}
	if (any(is.na(lambdaA1grid))){ 
		stop("Input (lambdaA1grid) is not a vector of non-negative numbers.") 
	}
	if (any(is.na(lambdaA2grid))){ 
		stop("Input (lambdaA2grid) is not a vector of non-negative numbers.") 
	}
	if (any(is.na(lambdaPgrid))){ 
		stop("Input (lambdaPgrid) is not a vector of non-negative numbers.") 
	}
	if (any(lambdaA1grid <= 0)){ 
		stop("Input (lambdaA1grid) is not a vector of non-negative numbers.") 
	}
	if (any(lambdaA2grid <= 0)){ 
		stop("Input (lambdaA2grid) is not a vector of non-negative numbers.") 
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

	# case where a grid for both autoregression penalty parameters is specified
	if ((length(lambdaA1grid) > 1) && 
	    (length(lambdaA2grid) > 1) && 
	    (length(lambdaPgrid) > 1)){
		lambdaPgrid <- lambdaPgrid[1]
	}
	if ((length(lambdaA1grid) > 1) && 
	    (length(lambdaA2grid) > 1) && 
	    (length(lambdaPgrid) == 1)){
		lambdaA1grid <- sort(lambdaA1grid)
		lambdaA2grid <- sort(lambdaA2grid)
		llLOOCV <- matrix(NA, nrow=length(lambdaA1grid), 
		                      ncol=length(lambdaA2grid))
		if (verbose) {
			cat("grid point:", "\n")
		}
		for (kA1 in 1:length(lambdaA1grid)) {
			for (kA2 in 1:length(lambdaA2grid)) {
				if (verbose) {
					cat(rep("\b", 100), sep = "")
					cat(paste("lambdaA1=", 
					          lambdaA1grid[kA1], 
					          "; lambdaA2=",
					          lambdaA2grid[kA2], 
					          sep=""))
				}
				llLOOCV[kA1, kA2] <- loglikLOOCVVAR2(c(lambdaA1grid[kA1],
				                                       lambdaA2grid[kA2],
				                                       lambdaPgrid[1]), 
				                                       Y, ...)
			}
		}
		if (figure) {
			contour(lambdaA1grid, 
			        lambdaA2grid, 
			        -llLOOCV, 
			        xlab="lambdaA1",
			        ylab="lambdaA2", 
				main="cross-validated log-likelihood")
		}
		return(list(lambdaA1=lambdaA1grid, 
		            lambdaA2=lambdaA2grid, 
		            lambdaPgrid=lambdaPgrid[1], 
		            llLOOCV=-llLOOCV))
	}

	# case where a grid for first autoregression 
	# and precision penalty parameters is specified
	if ((length(lambdaA1grid) > 1) && (length(lambdaA2grid) == 1) && (length(lambdaPgrid) > 1)){
		lambdaA1grid <- sort(lambdaA1grid)
		lambdaPgrid <- sort(lambdaPgrid)
		llLOOCV <- matrix(NA, nrow=length(lambdaA1grid), 
		                      ncol=length(lambdaPgrid))
		if (verbose) {
			cat("grid point:", "\n")
		}
		for (kA1 in 1:length(lambdaA1grid)) {
			for (kP in 1:length(lambdaPgrid)) {
				if (verbose) {
					cat(rep("\b", 100), sep = "")
					cat(paste("lambdaA1=", 
					          lambdaA1grid[kA1], 
					          "; lambdaP=",
					          lambdaPgrid[kP], sep=""))
				}
				llLOOCV[kA1, kP] <- loglikLOOCVVAR2(c(lambdaA1grid[kA1],
				                                      lambdaA2grid[1],
				                                      lambdaPgrid[kP]), 
				                                      Y, ...)
			}
		}
		if (figure) {
			contour(lambdaA1grid, 
			        lambdaPgrid, 
			        -llLOOCV, 
			        xlab="lambdaA1",
			        ylab="lambdaP", 
			        main="cross-validated log-likelihood")
		}
		return(list(lambdaA1=lambdaA1grid, 
		            lambdaA2=lambdaA2grid[1], 
		            lambdaP=lambdaPgrid, 
		            llLOOCV=-llLOOCV))
	}

	# case where a grid for second autoregression and 
	# precision penalty parameters is specified
	if ((length(lambdaA1grid) == 1) && 
            (length(lambdaA2grid) > 1) && 
	    (length(lambdaPgrid) > 1)){
		lambdaA2grid <- sort(lambdaA2grid)
		lambdaPgrid <- sort(lambdaPgrid)
		llLOOCV <- matrix(NA, nrow=length(lambdaA2grid), 
		                      ncol=length(lambdaPgrid))
		if (verbose) {
			cat("grid point:", "\n")
		}
		for (kA2 in 1:length(lambdaA2grid)) {
			for (kP in 1:length(lambdaPgrid)) {
				if (verbose) {
					cat(rep("\b", 100), sep="")
					cat(paste("lambdaA2=", 
					          lambdaA2grid[kA2], 
					          "; lambdaP=",
					          lambdaPgrid[kP], 
					          sep=""))
				}
				llLOOCV[kA2, kP] <- loglikLOOCVVAR2(c(lambdaA1grid[1],
				                                      lambdaA2grid[kA2],
				                                      lambdaPgrid[kP]), 
				                                      Y, ...)
			}
		}
		if (figure) {
			contour(lambdaA2grid, 
			        lambdaPgrid, 
			        -llLOOCV, 
			        xlab="lambdaA2",
			        ylab="lambdaP", 
			        main="cross-validated log-likelihood")
		}
		return(list(lambdaA1=lambdaA1grid[1], 
		            lambdaA2=lambdaA2grid, 
		            lambdaP=lambdaPgrid, 
		            llLOOCV=-llLOOCV))
	}
}

