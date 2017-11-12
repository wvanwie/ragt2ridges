optPenaltyVARX1 <- function (Y, 
                             X, 
                             lambdaMin, 
                             lambdaMax, 
                             lambdaInit=(lambdaMin+lambdaMax)/2,
                             optimizer="nlm", 
                             ...){ 
	########################################################################
	# 
	# DESCRIPTION: 
	# Automatic selection of the ridge penalties for the estimation of the 
	# parameters of the VARX(1) model. The penalties are selected through 
	# maximization of the LOOCV log-likelihood. 
	# 
	# ARGUMENTS:
	# -> Y             : Three-dimensional array containing the data. The 
	#                    first, second and third dimensions correspond to 
	#                    covariates, time and samples, respectively. The 
	#                    data are assumed to centered covariate-wise. 
	# -> X             : Three-dimensional array containing the time-varying 
	#                    covariate data. The first, second and third 
	#                    dimensions correspond to covariates, time and 
	#                    samples, respectively. The data are assumed to 
	#                    centered covariate-wise. 
	# -> lambdaMin     : Numeric of length three, containing the minimum 
	#                    values of ridge penalty parameters to be 
	#                    considered. The first element is the ridge 
	#                    parameter corresponding to the penalty on A, the 
	#                    matrix with auto-regression coefficients, the 
	#                    second to matrix B containing the regression 
	#                    coefficients of the time-varying covariates,
	#                    with with while the third parameter relates 
	#                    to the penalty on Omega, the precision matrix of 
	#                    the errors.
	# -> lambdaMax     : Numeric of length three, containing the maximum 
	#                    values of ridge penalty parameters to be considered. 
	#                    The first element is the ridge parameter 
	#                    corresponding to the penalty on A, the matrix with 
	#                    auto-regression coefficients, the second to matrix 
	#                    B containing the regression coefficients of the 
	#                    time-varying covariates, with with while the third 
	#                    parameter relates to the penalty on Omega, the 
	#                    precision matrix of the errors.
	# -> lambdaInit    : Numeric of length three, containing the initial 
	#                    values of ridge penalty parameters to be considered. 
	#                    The first element is the ridge parameter 
	#                    corresponding to the penalty on A, the matrix with 
	#                    auto-regression coefficients, the second to matrix 
	#                    B containing the regression coefficients of the 
	#                    time-varying covariates, while the third parameter 
	#                    relates to the penalty on Omega, the precision matrix 
	#                    of the errors.
	# -> optimizer     : A character : which optimization function should be 
	#                    used: "nlm" (default) or "optim"?
	# -> ...           : Additional arguments passed on to loglikLOOCVVARX1
	# 
	# DEPENDENCIES:
	# library(base)	        # functions: nlminb, constrOptim
	# library(ragt2ridges)  # functions: loglikLOOCVVARX1 and its dependencies.
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
	if (as.character(class(X)) != "array"){ 
		stop("Input (X) is of wrong class.") 
	}
	if (length(dim(X)) != 3){ 
		stop("Input (X) is of wrong dimensions: either covariate, time or sample dimension is missing.") 
	}
	if (any(dim(Y)[2:3] != dim(X)[2:3])){ 
		stop("Input (X) do not have same dimensions as Y.") 
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
	                                  loglikLOOCVVARX1, 
	                                  grad=NULL, 
	                                  ui=rbind(diag(rep(1, 3)), diag(rep(-1, 3))), 
	                                  ci=c(lambdaMin, -lambdaMax), 
	                                  Y=Y, 
	                                  X=X, 
	                                  control=list(reltol=0.001), 
	                                  ...)$par 
	}
	if (optimizer=="nlm"){    
		optLambdas <- nlminb(lambdaInit, 
	                             loglikLOOCVVARX1, 
	                             gradient=NULL, 
	                             lower=lambdaMin, 
	                             upper=lambdaMax, 
	                             Y=Y, 
	                             X=X, 
	                             control=list(rel.tol=0.001), 
	                             ...)$par 
	}
	return(optLambdas)
}


