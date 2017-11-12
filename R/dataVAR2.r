dataVAR2 <- function(n, T, A1, A2, SigmaE, TburnIn=1000){
	########################################################################
	#
	# DESCRIPTION:
	# simulate data from a VAR(2) process
	#
	# ARGUMENTS:
	# -> n        : Number of individuals to be sampled.
	# -> T        : Number of time points (per individual) to be sampled.
	# -> A1       : Matrix A1 with lag 1 auto-regression parameters.
	# -> A2       : Matrix A2 with lag 2 auto-regression parameters.	
	# -> SigmaE   : Covariance matrix of the errors (innovations).
	# -> TburnIn  : Number of time points to burn in the process.	
	#
	# DEPENDENCIES:
	# ...
	#
	# NOTES:
	# ...
	#
	########################################################################

	# input check
	if (as.character(class(A1)) != "matrix"){ 
		stop("Input (A1) is of wrong class.") 
	}
	if (as.character(class(A2)) != "matrix"){ 
		stop("Input (A2) is of wrong class.") 
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
	if (nrow(A1) != ncol(A1)){ 
		stop("Matrix A1 is not square.") 
	}
	if (nrow(A2) != ncol(A2)){ 
		stop("Matrix A2 is not square.") 
	}	
	if (nrow(A1) != nrow(SigmaE)){ 
		stop("Dimensions covariance matrix and A1 do not match.") 
	}
	if (as.character(class(n)) != "numeric"){ 
		stop("Input (n) is of wrong class.") 
	}
	if (length(n) != 1){ 
		stop("Input (n) is of wrong length.") 
	}
	if (is.na(n)){ 
		stop("Input (n) is not a positive integer.") 
	}
	if (n < 0){ 
		stop("Input (n) is not a positive integer.") 
	}
	if (as.character(class(T)) != "numeric"){ 
		stop("Input (T) is of wrong class.") 
	}
	if (length(T) != 1){ 
		stop("Input (T) is of wrong length.") 
	}
	if (is.na(T)){ 
		stop("Input (T) is not a positive integer.") 
	}
	if (T < 0){ 
		stop("Input (T) is not a positive integer.") 
	}
	if (as.character(class(TburnIn)) != "numeric"){ 
		stop("Input (TburnIn) is of wrong class.") 
	}
	if (length(TburnIn) != 1){ 
		stop("Input (TburnIn) is of wrong length.") 
	}
	if (is.na(TburnIn)){ 
		stop("Input (TburnIn) is not a positive integer.") 
	}
	if (TburnIn < 0){ 
		stop("Input (TburnIn) is not a positive integer.") 
	}

	# generate data
	# initiation
	Y        <- array(NA, dim=c(nrow(A1), T, n))
	Y[,1,]   <- matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A1), byrow=TRUE)
	Y[,2,]   <- A1 %*% Y[,1,] + 
	            matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A1), byrow=TRUE)
	Yupdate1 <- matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A1), byrow=FALSE)
	Yupdate2 <- A1 %*% Yupdate1 + 
		    matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A1), byrow=TRUE)	

	# burning in
	if (TburnIn > 0){
		for (t in 1:TburnIn){
			Yupdate3 <- A1 %*% Yupdate2 + 
			            A2 %*% Yupdate1 + 
				    matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A1), byrow=TRUE)
			Yupdate1 <- Yupdate2
			Yupdate2 <- Yupdate3
		}
	}

	# actual data generation
	Y[,1,] <- Yupdate1
	Y[,2,] <- Yupdate2
	for (t in 3:T){
        	Y[,t,] <- A1 %*% Y[,t-1,] + 
		          A2 %*% Y[,t-2,] + 
		          matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A1), byrow=TRUE)
	}
	return(Y)
}


