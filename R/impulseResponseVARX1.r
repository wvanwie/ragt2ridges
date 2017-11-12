impulseResponseVARX1 <- function(A, 
				 B, 
				 T){

	########################################################################
	#
	# DESCRIPTION:
	# Calculate impulse responses for VARX(1) model. MA representation of 
	# the VARX(1) model dictates that these are given by powers of A 
	# (post-multiplied by B for the impulse response of X).  
	#
	# ARGUMENTS:
	# -> A        : Matrix A of auto-regressive parameters.
	# -> B        : Matrix B of time-covarying parameters.
	# -> T        : Time points at which the impulse response is to be 
	#               evaluated.
	#
	# DEPENDENCIES:
	# library(expm)	    # functions: %^% 
	#
	# NOTES:
	# ...
	# 
	########################################################################

	# input checks
	if (as.character(class(A)) != "matrix"){ 
		stop("Input (A) is of wrong class.") 
	}
	if (nrow(A) != ncol(A)){ 
		stop("Matrix A is not square.") 
	}
	if (nrow(B) != nrow(A)){ 
		stop("Dimensions of matrices A and B do not match.") 
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
	if (T < 1){ 
		stop("Input (T) should be an integer larger than one.") 
	}

	# calculate impulse responses for the innovations
	if (T==0){ 
		impResponse <- diag(nrow(A)) 
	}
	if (T > 0){ 
		impResponse <- abs(A %^% T) 
	}

	# calculate impulse responses for the time-varying covariate X
	impResponseX <- impResponse %*% B
    
	return(list(impResponse=impResponse, impResponseX=impResponseX))   
}

