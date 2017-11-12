optPenaltyVAR1fused <- function (Y, 
				 id, 
				 lambdaMin, 
				 lambdaMax, 
				 lambdaInit=(lambdaMin+lambdaMax)/2, 
                                 optimizer="nlm",
				 ...){ 

	########################################################################
	# 
	# DESCRIPTION: 
	# Automatic selection of the ridge penalties for the estimation of the 
	# parameters of the VAR(2) model. The penalties are selected through 
	# maximization of the LOOCV log-likelihood. 
	# 
	# ARGUMENTS:
	# -> Y             : Three-dimensional array containing the data. The 
	#                    first, second and third dimensions correspond to 
	#                    covariates, time and samples, respectively. The 
	#                    data are assumed to centered covariate-wise. 
	# -> id            : A vector with groups indices comprising of integers 
	#                    only. First group is represented by '0', the next 
	#                    by '1', and so on until the last.	
	# -> lambdaMin     : Numeric of length three, containing the minimum 
	#                    values of ridge penalty parameters to be considered. 
	#                    The first element is the ridge parameter 
	#                    corresponding to the penalty on all A_g's, the 
	#                    matrices with lag one auto-regression coefficients, 
	#                    the second to the fused ridge parameter for these 
	#                    A_g' s, while the third parameter relates to the 
	#                    penalty on Omega, the precision matrix of the
	#                    errors.
	# -> lambdaMax     : Numeric of length three, containing the maximum 
	#                    values of ridge penalty parameters to be considered. 
	#                    The first element is the ridge parameter 
	#                    corresponding to the penalty on all A_g's, the 
	#                    matrices with lag one auto-regression coefficients, 
	#                    the second to the fused ridge parameter for these 
	#                    A_g's, while the third parameter relates to the 
	#                    penalty on Omega, the precision matrix of the 
	#                    errors.
	# -> lambdaInit    : Numeric of length three, containing the initial 
	#                    values of ridge penalty parameters to be considered. 
	#                    The first element is the ridge parameter 
	#                    corresponding to the penalty on all A_g's, the 
	#                    matrices with lag one auto-regression coefficients, 
	#                    the second to the fused ridge parameter for these 
	#                    A_g's, while the third parameter relates to the 
	#                    penalty on Omega, the precision matrix of the 
	#                    errors.
	# -> optimizer     : A character : which optimization function should be
	#                    used: "nlm" (default) or "optim"?
	# -> ...           : Additional arguments passed on to loglikLOOCVVAR2
	# 
	# DEPENDENCIES:
	# library(base)	        # functions: nlminb, constrOptim
	# library(ragt2ridges)  # functions: loglikLOOCVVAR2 and its dependencies.
	#
	# NOTES:
	# ....
	# 
	########################################################################

	# input checks
	if (as.character(class(Y)) != "array"){ 
		stop("Input (Y) is of wrong class.") 
	}
	if (length(dim(Y)) != 3){ 
		stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") 
	}
	if (as.character(class(id)) != "numeric" & as.character(class(id)) != "integer"){ 
		stop("Input (id) is of wrong class.") 
	}
	if (length(id) != dim(Y)[3]){ 
		stop("Input (id) is of wrong length: should equal sample dimension of Y.") 
	}		
	if (as.character(class(lambdaMin)) != "numeric"){ 
		stop("Input (lambdaMin) is of wrong class.") 
	}
	if (length(lambdaMin) != 3){ 
		stop("Input (lambdaMin) is of wrong length.") 
	}
	if (any(is.na(lambdaMin))){ 
		stop("Input (lambdaMin) does not comprise of positive numbers.") 
	}
	if (any(lambdaMin < 0)){ 
		stop("Input (lambdaMin) does not comprise of positive numbers.") 
	}
	if (as.character(class(lambdaMax)) != "numeric"){ 
		stop("Input (lambdaMax) is of wrong class.") 
	}
	if (length(lambdaMax) != 3){ 
		stop("Input (lambdaMax) is of wrong length.") 
	}
	if (any(is.na(lambdaMax))){ 
		stop("Input (lambdaMax) does not comprise of positive numbers.") 
	}
	if (any(lambdaMax < 0)){ 
		stop("Input (lambdaMax) does not comprise of positive numbers.") 
	}
	if (any(lambdaMax <= lambdaMin)){ 
		stop("Input (lambdaMax) must be larger (element-wise) than lambdaMin") 
	}
	if (as.character(class(lambdaInit)) != "numeric"){ 
		stop("Input (lambdaInit) is of wrong class.") 
	}
	if (length(lambdaInit) != 3){ 
		stop("Input (lambdaInit) is of wrong length.") 
	}
	if (any(is.na(lambdaInit))){ 
		stop("Input (lambdaInit) does not comprise of positive numbers.") 
	}
	if (any(lambdaInit < 0)){ 
		stop("Input (lambdaInit) does not comprise of positive numbers.") 
	}
	if (any(lambdaInit <= lambdaMin)){ 
		stop("Input (lambdaInit) must be larger (element-wise) than lambdaMin") 
	}
	if (any(lambdaInit >= lambdaMax)){ 
		stop("Input (lambdaInit) must be smaller (element-wise) than lambdaMax") 
	}

	# optimize LOOCV log-likelihood w.r.t. the penalty parameters
	if (optimizer=="optim"){    
		optLambdas <- constrOptim(lambdaInit, 
				          loglikLOOCVVAR1fused, 
					  grad=NULL, 
					  ui=rbind(diag(rep(1, 3)), diag(rep(-1, 3))), 
					  ci=c(lambdaMin, -lambdaMax), 
					  Y=Y, 
					  id=id,  
					  control=list(reltol=0.001), 
					  ...)$par 
	}
	if (optimizer=="nlm"){    
		optLambdas <- nlminb(lambdaInit, 
				     loglikLOOCVVAR1fused, 
				     gradient=NULL, 
				     lower=lambdaMin, 
				     upper=lambdaMax, 
				     Y=Y, 
		                     id=id, 
				     control=list(rel.tol=0.001), 
				     ...)$par 
	}
	return(optLambdas)
}


