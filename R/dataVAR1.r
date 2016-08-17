dataVAR1 <- function(n, T, A, SigmaE, TburnIn=1000){
	#############################################################################
	#
	# DESCRIPTION:
	# simulate data from a VAR(1) process
	#
	# ARGUMENTS:
	# -> n        : Number of individuals to be sampled.
	# -> T        : Number of time points (per individual) to be sampled.
	# -> A        : Matrix A of regression parameters.
	# -> SigmaE   : Covariance matrix of the errors (innovations).
	# -> TburnIn  : Number of time points to burn in the process.
	#
	# DEPENDENCIES:
	# ...
	#
	# NOTES:
	# ...
	#
	#############################################################################

	# input check
	if (as.character(class(A)) != "matrix"){ stop("Input (A) is of wrong class.") }
	if (as.character(class(SigmaE)) != "matrix"){ stop("Input (SigmaE) is of wrong class.") }
	if (!isSymmetric(SigmaE)){ stop("Non-symmetric covariance matrix is provided.") }
	if (!all(eigen(SigmaE)$values > 0)){ stop("Non positive-definite covariance matrix is provided.") }
	if (nrow(A) != ncol(A)){ stop("Matrix A is not square.") }
	if (nrow(A) != nrow(SigmaE)){ stop("Dimensions covariance matrix and A do not match.") }
	if (as.character(class(n)) != "numeric"){ stop("Input (n) is of wrong class.") }
	if (length(n) != 1){ stop("Input (n) is of wrong length.") }
	if (is.na(n)){ stop("Input (n) is not a positive integer.") }
	if (n < 0){ stop("Input (n) is not a positive integer.") }
	if (as.character(class(T)) != "numeric"){ stop("Input (T) is of wrong class.") }
	if (length(T) != 1){ stop("Input (T) is of wrong length.") }
	if (is.na(T)){ stop("Input (T) is not a positive integer.") }
	if (T < 0){ stop("Input (T) is not a positive integer.") }
	if (as.character(class(TburnIn)) != "numeric"){ stop("Input (TburnIn) is of wrong class.") }
	if (length(TburnIn) != 1){ stop("Input (TburnIn) is of wrong length.") }
	if (is.na(TburnIn)){ stop("Input (TburnIn) is not a positive integer.") }
	if (TburnIn < 0){ stop("Input (TburnIn) is not a positive integer.") }

	# warn user if the provided parameters do not correspond to a stationary process
	# .isStationary(A)

	# generate data
	Y <- array(NA, dim=c(nrow(A), T, n))
	Yupdate <- matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A), byrow=FALSE)
	if (TburnIn > 0){
		for (t in 1:TburnIn){
	        	Yupdate <- A %*% Yupdate + matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A), byrow=TRUE)
		}
	}
	Y[,1,] <- Yupdate
	for (t in 2:T){
        	Y[,t,] <- A %*% Y[,t-1,] + matrix(rmvnorm(n, sigma=SigmaE), nrow=nrow(A), byrow=TRUE)
	}
	return(Y)
}


