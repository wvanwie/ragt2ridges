mutualInfoVAR2 <- function(A1, 
                           A2, 
                           SigmaE, 
                           T){

	########################################################################
	#
	# DESCRIPTION:
	# Calculate the mutual informations for VAR(2) model.
	#
	# ARGUMENTS:
	# -> A1       : Matrix A1 with lag 1 auto-regression parameters.
	# -> A2       : Matrix A2 with lag 2 auto-regression parameters.	
	# -> SigmaE   : Covariance matrix of the errors (innovations).
	# -> T        : Time point at which the mutual information is to 
	#               be evaluated.
	# 
	# DEPENDENCIES:
	# library(expm)	    # functions: %^%
	#
	# NOTES:
	# ....
	# 	
	########################################################################

	# input checks
	if (as.character(class(A1)) != "matrix"){ 
		stop("Input (A1) is of wrong class.") 
	}
	if (nrow(A1) != ncol(A1)){ 
		stop("Matrix A1 is not square.") 
	}
	if (as.character(class(A2)) != "matrix"){ 
		stop("Input (A2) is of wrong class.") 
	}
	if (nrow(A2) != ncol(A2)){ 
		stop("Matrix A2 is not square.") 
	}
	if (nrow(A1) != nrow(A2)){ 
		stop("Matrices A1 and A2 do not have same dimensions.") 
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
	if (as.character(class(SigmaE)) != "matrix"){ 
		stop("Input (SigmaE) is of wrong class.") 
	}    
	if (!isSymmetric(SigmaE)){ 
		stop("Non-symmetrical matrix for the error covariance matrix provided") 
	} 
	if (nrow(A1) != nrow(SigmaE)){ 
		stop("Dimensions of input (A1) do not match that of other input (SigmaE).") 
	} 
	if (ncol(A1) != ncol(SigmaE)){ 
		stop("Dimensions of input (A1) do not match that of other input (SigmaE).") 
	} 
       
	# calculate impulse responses
	MI <- function(j, A1, A2, SigmaE, T){
		# covariance matrix with right zero
		SigmaEnull <- SigmaE
		SigmaEnull[j,] <- SigmaEnull[,j] <- 0
		SigmaEnull[-j,-j] <- SigmaEnull[-j,-j] - 
		                     SigmaE[-j, j, drop=FALSE] %*% 
                                     solve(SigmaE[j,j]) %*% 
		                     SigmaE[j,-j,drop=FALSE]

		# initial variances
		varMarg1 <- SigmaE
		varMarg2 <- SigmaE + 
		            A1 %*% SigmaE %*% t(A1)
		varMarg3 <- SigmaE + 
		            A1 %*% SigmaE %*% t(A1) + 
		            A2 %*% SigmaE %*% t(A2) + 
		            A1 %*% A1 %*% SigmaE %*% t(A1) %*% t(A1)         
		varCond1 <- SigmaEnull
		varCond2 <- SigmaE + 
		            A1 %*% SigmaEnull %*% t(A1)
		varCond3 <- SigmaE + 
		            A1 %*% SigmaE %*% t(A1) + 
		            A2 %*% SigmaEnull %*% t(A2) + 
		            A1 %*% A1 %*% SigmaEnull %*% t(A1) %*% t(A1)

		# special cases
		if (T == 1){ 
			MIslh <- determinant(varMarg1)$modulus - 
			         determinant(varCond1)$modulus 
		}
		if (T == 2){ 
			MIslh <- determinant(varMarg2)$modulus - 
			         determinant(varCond2)$modulus 
		}
		if (T == 3){
			MIslh <- determinant(varMarg3)$modulus - 
			         determinant(varCond3)$modulus 
		}
		if (T > 3){
			# initiating covariances
			covMarg10 <- matrix(0, nrow(A1), ncol(A1))
			covMarg21 <- A1 %*% SigmaE
			covCond32 <- A1 %*% SigmaE + 
			             A2 %*% SigmaE %*% t(A1) + 
			             A1 %*% A1 %*% SigmaE %*% t(A1)         
			covCond10 <- matrix(0, nrow(A1), ncol(A1))
			covCond21 <- A1 %*% SigmaEnull
			covCond32 <- A1 %*% SigmaE + 
			             A2 %*% SigmaEnull %*% t(A1) + 
			             A1 %*% A1 %*% SigmaEnull %*% t(A1)         

			covMargT_2andT_3 <- A1 %*% SigmaE
			covMargT_1andT_2 <- A1 %*% SigmaE + 
			                    A2 %*% SigmaE %*% t(A1) + 
			                    A1 %*% A1 %*% SigmaE %*% t(A1)         
			varMargT_1 <- varMarg3
			varMargT_2 <- varMarg2
			for (tau in 4:T){        
				varMargT <- SigmaE + 
				            A1 %*% varMargT_1 %*% t(A1) + 
				            A2 %*% varMargT_2 %*% t(A2) + 
				            A1 %*% covMargT_1andT_2 %*% t(A2) + 
				            t(A1 %*% covMargT_1andT_2 %*% t(A2))
				covMargTandT_1 <- A1 %*% varMargT_1 + 
				                  A2 %*% varMargT_1 %*% t(A1) + 
				                  A2 %*% covMargT_2andT_3 %*% t(A2)
				varMargT_2 <- varMargT_1
				varMargT_1 <- varMargT
				covMargT_2andT_3 <- covMargT_1andT_2
				covMargT_1andT_2 <- covMargTandT_1
			}        

			covCondT_2andT_3 <- A1 %*% SigmaEnull
            		covCondT_1andT_2 <- A1 %*% SigmaE + 
			                    A2 %*% SigmaEnull %*% t(A1) + 
			                    A1 %*% A1 %*% SigmaEnull %*% t(A1)        
			varCondT_1 <- varCond3
			varCondT_2 <- varCond2
			for (tau in 4:T){        
				varCondT <- SigmaE + 
				            A1 %*% varCondT_1 %*% t(A1) + 
				            A2 %*% varCondT_2 %*% t(A2) + 
				            A1 %*% covCondT_1andT_2 %*% t(A2) + 
				            t(A1 %*% covCondT_1andT_2 %*% t(A2))
				covCondTandT_1 <- A1 %*% varCondT_1 + 
				                  A2 %*% varCondT_1 %*% t(A1) + 
				                  A2 %*% covCondT_2andT_3 %*% t(A2)
				varCondT_2 <- varCondT_1
				varCondT_1 <- varCondT
				covCondT_2andT_3 <- covCondT_1andT_2
				covCondT_1andT_2 <- covCondTandT_1
			}        
			MIslh <- as.numeric(determinant(varMargT)$modulus) - 
				 as.numeric(determinant(varCondT)$modulus) 
		}
		return(MIslh)
	}    
	MIs <- unlist(lapply(1:nrow(A1), MI, A1=A1, A2=A2, SigmaE=SigmaE, T=T))
        
	return(MIs)   
}


