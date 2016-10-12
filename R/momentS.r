momentS <- function(Sigma, shape, moment=1){
	############################################################################################################
	#
	# DESCRIPTION:
	# Returns the moments of a Wishart-distributed random variable.
	# Only those explicitly given in Lesac, Massam (2004) are implemented.
	#
	# ARGUMENTS:
	# -> Sigma	: Positive-definite 'matrix', the scale parameter of the Wishart distribution.
	# -> shape	: A 'numeric', the shape parameter of the Wishart distribution. Should exceed the number 
	# 		  of variates.
	# -> moment	: An 'integer'. Should be in the set {-4, -3, -2, -1, 0, 1, 2, 3, 4} (only those are explicitly 
	# 		  specified in Lesac, Massam, 2004).
	# 
	# DEPENDENCIES:
	# Currently, none.
	#
	# NOTES:
   	# ....
	#     	
	############################################################################################################

	# input checks
	if (shape < nrow(Sigma) + 1){ stop("shape parameter should exceed the number of variates") }
	if (all(moment != c(-c(4:1), 0:4))){ stop("moment not implemented") }

	# number of variates
	p <- nrow(Sigma)	

	# moment, case-wise
	if (moment==-4){ 
		Sinv <- solve(Sigma); 
		constant <- (shape-p+2) * (shape-p+1) * (shape-p-2) * (shape-p) * (shape-p-7) * (shape-p-5) * (shape-p-3) * (shape-p-1);
		c1 <- 5 * (shape-p) - 11;
		c2 <- 5 * (shape-p)^2 - 16*(shape-p) + 11;
		c3 <- (shape-p)^3 - 4*(shape-p)^2 + 7*(shape-p) - 4;
		c4 <- 2*(shape-p)^3 - 9*(shape-p)^2  + 24*(shape-p) + 19;
		c5 <- (shape-p)^4 - 5*(shape-p)^3 + 11*(shape-p)^2 - 11*(shape-p) + 4;
		ESr <- (c1 * Sinv * (sum(diag(Sinv)))^3  + 
	        	c2 * (Sinv * sum(diag(Sinv)) * sum(Sinv * Sinv) + Sinv %*% Sinv * (sum(diag(Sinv)))^2) +
	        	c3 * (Sinv * sum(diag(Sinv %*% Sinv %*% Sinv)) +  3 * Sinv %*% Sinv %*% Sinv * sum(diag(Sinv))) + 
		        c4 * Sinv %*% Sinv * sum(Sinv * Sinv) + c5 * Sinv %*% Sinv %*% Sinv %*% Sinv) / constant 
	}	
	if (moment==-3){ 
		Sinv <- solve(Sigma); 
		constant <- ((shape - p) * (shape - p - 1) * (shape - p - 3) * (shape - p + 1) * (shape - p - 5))
		ESr <- ( (2 * (.trace(Sinv))^2 * Sinv + (shape - p - 1) * (.trace(Sinv %*% Sinv) * Sinv + 2 * .trace(Sinv) * Sinv %*% Sinv) + 
			(shape - p - 1)^2 * Sinv %*% Sinv %*% Sinv) / constant )
	}
	if (moment==-2){ 
		Sinv <- solve(Sigma); 
		ESr <- (((shape - p - 1) * Sinv %*% Sinv + .trace(Sinv) * Sinv) / ((shape - p) * (shape - p - 1) * (shape - p - 3)))
	}
	if (moment==-1){ ESr <- solve(Sigma) / (shape - p - 1) }
	if (moment==0){ ESr <- diag(p) }
	if (moment==1){ ESr <- shape * Sigma }
	if (moment==2){ ESr <- (shape * (shape + 1) * Sigma %*% Sigma + shape * .trace(Sigma) * Sigma) }
	if (moment==3){ ESr <- (shape * (.trace(Sigma))^2 * Sigma + shape * (shape + 1) * (.trace(Sigma %*% Sigma) * Sigma + 
                            2 * .trace(Sigma) * Sigma %*% Sigma) + shape * (shape^2 + 3*shape + 4) * Sigma %*% Sigma %*% Sigma) }
	if (moment==4){ ESr <- (shape * .trace(Sigma %*% Sigma %*% Sigma) * Sigma + 3 * shape * (shape + 1) * (.trace(Sigma) * .trace(Sigma %*% Sigma) * Sigma + 
                            (.trace(Sigma))^2 * Sigma %*% Sigma) + 
                            shape * (shape^2 + 3 * shape + 4) * (3* .trace(Sigma) * Sigma %*% Sigma %*% Sigma + .trace(Sigma %*% Sigma %*% Sigma) * Sigma) +                   
                            shape * (2 * shape^2 + 5 * shape + 5) * .trace(Sigma %*% Sigma) * Sigma %*% Sigma + 
                            shape * (shape^3 + 6 * shape^2 + 21 * shape + 20) * Sigma %*% Sigma %*% Sigma %*% Sigma
                            ) }
    return(ESr / shape^moment)
}

