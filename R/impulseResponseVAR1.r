impulseResponseVAR1 <- function(A, 
				T){

	########################################################################
	#
	# DESCRIPTION:
	# Calculate impulse response for VAR(1) model. MA representation of the 
	# VAR(1) model dictates that these are given by powers of A.
	#
	# ARGUMENTS:
	# -> A       : Matrix A of regression parameters.
	# -> T       : Time point for which the impulse response is to be 
	#              evaluated.
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

	# calculate impulse responses
	return(A %^% T)   
}

