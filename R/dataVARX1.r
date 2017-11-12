dataVARX1 <- function(X, A, B, SigmaE, lagX){
	########################################################################
	#
	# DESCRIPTION:
	# simulate data from a VARX(1) process
	#
	# ARGUMENTS:
	# -> X        : Matrix X of time-varying covariates.
	# -> A        : Matrix A of auto-regressive parameters.
	# -> B        : Matrix B of time-covarying parameters.
	# -> SigmaE   : Covariance matrix of the errors (innovations).
	# -> lagX     : Integer specifying how whether Xt of Xt-1 affects Yt.
	#
	# DEPENDENCIES:
	# ...
	#
	# NOTES:
	# ...
	#
	########################################################################

	# input check
	if (as.character(class(A)) != "matrix"){ 
		stop("Input (A) is of wrong class.") 
	}
	if (as.character(class(B)) != "matrix"){ 
		stop("Input (A) is of wrong class.") 
	}
	if (as.character(class(SigmaE)) != "matrix"){ 
		stop("Input (SigmaE) is of wrong class.") 
	}
	if (!isSymmetric(SigmaE)){ 
		stop("Non-symmetric covariance matrix is provided.") 
	}
	if (!all(eigen(SigmaE)$values > 0)){ 
		stop("Non positive-definite covariance matrix is provided.") 
	}
	if (nrow(A) != ncol(A)){ 
		stop("Matrix A is not square.") 
	}
	if (nrow(A) != nrow(SigmaE)){ 
		stop("Dimensions covariance matrix and A do not match.") 
	}
	if (nrow(B) != nrow(A)){ 
		stop("Dimensions of matrices A and B do not match.") 
	}
	if (length(lagX) != 1){ 
		stop("Input (lagX) is of wrong length.") 
	}
	if (is.na(lagX)){ 
		stop("Input (n) is an integer.") 
	}
	if (lagX != 0 & lagX != 1){ 
		stop("Input (lagX) does not equal zero or one.") 
}

	# warn user if the provided parameters do not correspond to a stationary process
	# .isStationary(A)

	# generate data
	Y <- array(NA, dim=c(nrow(A), dim(X)[2], dim(X)[3]))
	Y[,1,] <- matrix(rmvnorm(dim(X)[3], sigma=SigmaE), nrow=dim(Y)[1], byrow=TRUE)
	for (t in 2:dim(X)[2]){
        	Y[,t,] <- A %*% Y[,t-1,] + 
			  B %*% X[,t-lagX,] + 
			  matrix(rmvnorm(dim(X)[3], sigma=SigmaE), nrow=dim(A)[1], byrow=TRUE)
	}
	return(Y)
}

