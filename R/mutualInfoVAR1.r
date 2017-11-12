mutualInfoVAR1 <- function(A, 
                           SigmaE, 
	                   T){

	########################################################################
	#
	# DESCRIPTION:
	# Calculate the mutual informations for VARX(1) model.
	#
	# ARGUMENTS:
	# -> A        : Matrix A of regression parameters.
	# -> SigmaE   : Covariance matrix of the errors (innovations).
	# -> T        : Time points at which the mutual information is to 
	#               be evaluated.
	# 
	# DEPENDENCIES:
	# library(expm)	    # functions: %^%
	#
	# NOTES:
	# ....
	# 	
	########################################################################

	# input checks
	if (as.character(class(A)) != "matrix"){ 
		stop("Input (A) is of wrong class.") 
	}
	if (nrow(A) != ncol(A)){ 
		stop("Matrix A is not square.") 
	}
	if (as.character(class(T)) != "numeric"){ 
		stop("Input (Tmax) is of wrong class.") 
	}
	if (length(T) != 1){ 
		stop("Input (T) is of wrong length.") 
	}
	if (is.na(T)){ 
		stop("Input (T) is not a positive integer.") 
	}
	if (T < 1){ 
		stop("Input (T) should be an integer larger than one.") 
	}
	if (as.character(class(SigmaE)) != "matrix"){ 
		stop("Input (SigmaE) is of wrong class.") 
	}    
	if (!isSymmetric(SigmaE)){ 
		stop("Non-symmetrical matrix for the error covariance matrix provided") 
	} 
	if (nrow(A) != nrow(SigmaE)){ 
		stop("Dimensions of input (A) do not match that of other input (SigmaE).") 
	} 
	if (ncol(A) != ncol(SigmaE)){ 
		stop("Dimensions of input (A) do not match that of other input (SigmaE).") 
	} 

	# calculate impulse responses
	MI <- function(j, A, SigmaE, T){
		# covariance matrix with right zero
		SigmaEnull <- SigmaE
		SigmaEnull[j,] <- SigmaEnull[,j] <- 0
		SigmaEnull[-j,-j] <- SigmaEnull[-j,-j] - 
		                     SigmaE[-j, j, drop=FALSE] %*% 
                                            solve(SigmaE[j,j]) %*% 
                                            SigmaE[j,-j,drop=FALSE]

		# initial value of marginal and conditional covariances
		varMarg <- varCond <- SigmaE 
		Atau <- A %^% T
		varCond <- varMarg + Atau %*% SigmaEnull %*% t(Atau)
		varMarg <- varMarg + Atau %*% SigmaE %*% t(Atau)
		MIslh <- log(det(varMarg)) - log(det(varCond))
		return(MIslh)
	}    

	# reformat object
	MIs <- unlist(lapply(1:nrow(A), MI, A=A, SigmaE=SigmaE, T=T))
         
	# return object
	return(MIs)   
}


