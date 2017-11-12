loglikLOOCVVARX1 <- function(lambdas, 
                             Y, 
                             X, 
                             unbalanced=matrix(nrow=0, ncol=2), 
                             lagX=0, 
                             ...){

	########################################################################
	# 
	# DESCRIPTION:
	# Evaluation of the (minus) leave-one-out cross-validated 
	# log-likelihood of the VAR(1) model for given choices of the ridge
	# penalty parameters (lambdaA and lambdaO). The functions also works 
	# with a (possibly) unbalanced experimental set-up. The VAR(1)-process 
	# is assumed to have mean zero.
	#
	# ARGUMENTS:
	# -> lambdas      : Numeric of length three. It contains the ridge 
	#                   penalty parameters to be used in the estimation
	#                   of A, B and the Omega, the precision matrix of 
	#                   the errors (also called innovations). 
	# -> Y            : Three-dimensional array containing the data. The 
	#                   first, second and third dimensions correspond to 
	#                   covariates, time and samples, respectively. The 
	#                   data are assumed to centered covariate-wise. 
	# -> X            : Three-dimensional array containing the 
	#                   time-varying covariate data. The first, second 
	#                   and third dimensions correspond to covariates, 
	#                   time and samples, respectively. The data are 
	#                   assumed to be centered covariate-wise. 
	# -> unbalanced   : A matrix with two columns, indicating the 
	#                   unbalances in the design. Each row represents 
	#                   a missing design point in the 
	#                   (time x individual)-layout. The first and 
	#                   second column indicate the time and individual 
	#                   (respectively) specifics of the missing design 
	#                   point.
	# -> ...          : Other arguments to be passed to ridgeVAR1.
	#
	# DEPENDENCIES:
	# require("ragt2ridges")          # functions from package : 
	#                                   ridgeVAR1
	#
	# NOTES:
	# ...
	#
	########################################################################

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
	if (any(dim(X)[2:3] != dim(Y)[2:3])){ 
		stop("Inputs (X and Y) do not have same dimensions.") 
	}
	if (as.character(class(lambdas)) != "numeric"){ 
		stop("Input (lambdas) is of wrong class.") 
	}
	if (length(lambdas) != 3){ 
		stop("Input (lambdas) is of wrong length.") 
	}
	if (any(is.na(lambdas))){ 
		stop("Input (lambdas) is not a vector of positive numbers.") 
	}
	if (any(lambdas <= 0)){ 
		stop("Input (lambdas) is not a vector of positive numbers.") 
	}
	if (as.character(class(unbalanced)) != "matrix"){ 
		stop("Input (unbalanced) is of wrong class.") 
	}    
	if (ncol(unbalanced) != 2){ 
		stop("Wrong dimensions of the matrix unbalanced.") 
	} 

	# determine leave-one-out scheme
	LOOscheme <- cbind(rep(2:dim(Y)[2], dim(Y)[3]), 
				sort(rep(1:dim(Y)[3], dim(Y)[2]-1)))	
	if (nrow(unbalanced) > 0){
		LOO2unbalanced <- numeric()
		for (k in 1:nrow(unbalanced)){
			LOO2unbalanced <- c(LOO2unbalanced, 
				which(apply(LOOscheme, 1, 
					function(Y, Z){ all(Y == Z) }, 
					Z=unbalanced[k,])))
		}
		LOOscheme <- LOOscheme[-LOO2unbalanced,]
	}

	loglik <- 0
	for (k in 1:nrow(LOOscheme)){
		# obtain LOOCV estimates of VARX(1) parameters A, B and Se
		VARX1hat <- ridgeVARX1(Y, 
					X, 
					lambdas[1], 
					lambdas[2], 
					lambdas[3], 
					unbalanced=rbind(unbalanced, 
							LOOscheme[k,,drop=FALSE]), 
					lagX=lagX, ...)

		# evaluate LOOCV loglikelihood
		loglik <- loglik + 
				.armaVARX1_loglik_LOOCVinternal(
					Y[,LOOscheme[k,1], LOOscheme[k,2]], 
					Y[,LOOscheme[k,1]-1, LOOscheme[k,2]], 
					X[,LOOscheme[k,1]-lagX, LOOscheme[k,2]], 
					VARX1hat$A, 
					VARX1hat$B, 
					VARX1hat$P)
	}

	# return minus LOOCV loglikelihood
	return(-loglik)
}



