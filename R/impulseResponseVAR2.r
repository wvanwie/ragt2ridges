impulseResponseVAR2 <- function(A1, 
				A2, 
				T){

	########################################################################
	#
	# DESCRIPTION:
	# Calculate impulse responses for VARX(1) model. MA representation of 
	# the VARX(1) model dictates that these are given by powers of A 
	# (post-multiplied by B for the impulse response of X).  
	#
	# ARGUMENTS:
	# -> A1       : Matrix A1 with lag 1 auto-regression parameters.
	# -> A2       : Matrix A2 with lag 2 auto-regression parameters.	
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
	if (as.character(class(A1)) != "matrix"){ 
		stop("Input (A1) is of wrong class.") 
	}
	if (as.character(class(A2)) != "matrix"){ 
		stop("Input (A2) is of wrong class.") 
	}
	if (nrow(A1) != ncol(A1)){ 
		stop("Matrix A1 is not square.") 
	}
	if (nrow(A2) != ncol(A2)){ 
		stop("Matrix A2 is not square.") 
	}	
	if (nrow(A1) != nrow(A2)){ 
		stop("Matrix A1 and A2 are not of same dimensions.") 
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
	impResponseTmin2 <- diag(nrow(A1)) 
	impResponseTmin1 <- A1 
	if (T==0){ impResponse <- impResponseTmin2 }
	if (T==1){ impResponse <- impResponseTmin1 }
	for (t in 2:T){
		impResponse <- A1 %*% impResponseTmin1 + A2 %*% impResponseTmin2
		impResponseTmin2 <- impResponseTmin1	    
		impResponseTmin1 <- impResponse	        
	}
	return(impResponse)   
}

