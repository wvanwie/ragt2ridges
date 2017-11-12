CIGofVAR2 <- function(sparseA1, 
                      sparseA2, 
                      sparseP, 
                      type="global"){

	#######################################################################
	#
	# DESCRIPTION:
	# Constructs the global or contemporaneous conditional independence 
	# graph (CIG) of the VAR(1) model.
	#
	# ARGUMENTS:
	# -> sparseA1                : Matrix A1 of regression parameters, 
	#                              which is assumed to be sparse.
	# -> sparseA2                : Matrix A2 of regression parameters, 
	#                              which is assumed to be sparse.
	# -> sparseP                 : Matrix P of precision of the error, 
	#                              which is assumed to be sparse.
	# -> type                    : A 'character' indicating whether 
	#                              the 'global' or 'contemp'
	#                              (contemporaneous)adjanceny matrix of 
	#                              the conditional independence graph is 
	#                              requested.
	# 
	# DEPENDENCIES:
	# ....
	#
	# NOTES:
	# ....
	#
	# REFERENCES: 
	# -> Dahlhaus (2000), "Graphical interaction models for multivariate 
	#    time series", Metrika, 51, 157-172.
	# -> Dahlhaus, Eichler (2003), "Causality and graphical models in 
	#    time series analysis", Oxford Statistical Science Series, 115-137. 
	#     	
	#######################################################################
    
	# input checks    
  	if (as.character(class(sparseA1)) != "matrix"){ 
		stop("Input (sparseA1) is of wrong class.") 
	}
	if (nrow(sparseA1) != ncol(sparseA1)){ 
		stop("Matrix sparseA1 is not square.") 
	}
  	if (as.character(class(sparseA2)) != "matrix"){ 
		stop("Input (sparseA2) is of wrong class.") 
	}
	if (nrow(sparseA2) != ncol(sparseA2)){ 
		stop("Matrix sparseA2 is not square.") 
	}
	if (as.character(class(sparseP)) != "matrix"){ 
		stop("Input (sparseP) is of wrong class.") 
	}
	if (nrow(sparseP) != ncol(sparseP)){ 
		stop("Matrix sparseP is not square.") 
	}
	if (nrow(sparseA1) != ncol(sparseP)){ 
		stop("Matrix sparseA1 and sparseP are not of equal dimensions.") 
	}
	if (nrow(sparseA2) != ncol(sparseP)){ 
		stop("Matrix sparseA2 and sparseP are not of equal dimensions.") 
	}
  	if (as.character(class(type)) != "character"){ 
		stop("Input (type) is of wrong class.") 
	}
  	if (!(type %in% c("global", "contemp"))){ 
		stop("Input (type) ill-specified.") 
	}

	# construct adjacency matrix of global markov (in)dependencies
	if (type=="global"){
		# restriction of Dahlhaus and Eichler, actually a slightly more restrictive version.
		CIG <- abs(t(sparseA1) %*% sparseP) + 
			abs(sparseP %*% sparseA1) + 
			abs(t(sparseA1) %*% sparseP %*% sparseA1) +
			abs(t(sparseA2) %*% sparseP) + 
			abs(sparseP %*% sparseA2) + 
			abs(t(sparseA2) %*% sparseP %*% sparseA2) +
			abs(t(sparseA1) %*% sparseP %*% sparseA2) + 
			abs(t(sparseA2) %*% sparseP %*% sparseA1)
		diag(sparseP) <- 0
		CIG <- CIG + abs(sparseP) 
		CIG[CIG != 0] <- 1
	}
    
	# construct adjacency matrix of global markov (in)dependencies
	if (type=="contemp"){
		diag(sparseP) <- 0
		CIG <- abs(sparseP) 
		CIG[CIG != 0] <- 1
	}
	return(CIG)
}

