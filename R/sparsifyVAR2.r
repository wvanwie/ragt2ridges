sparsifyVAR2 <- function(A1, 
			 A2, 
			 SigmaE, 
			 threshold=c("absValue", "localFDR", "top"), 
			 absValueCut=c(0.25, 0.25), 
			 FDRcut=c(0.8, 0.8), 
			 top=c(10,10), 
			 zerosA1=matrix(nrow=0, ncol=2), 
			 zerosA2=matrix(nrow=0, ncol=2), 
			 statistics=FALSE, 
			 verbose=TRUE){

	########################################################################
	#
	# DESCRIPTION: 
	# Function that sparsifies/determines support of the matrix A of a 
	# VAR(1) model. Support can be determined by absolute value thresholding 
	# or by local FDRs thresholding. One can also choose to threshold based
	# on the top X of absolute value of the elements of A. Function is to 
	# some extent a wrapper around certain 'fdrtool' functions
	#
	# ARGUMENTS:
	# -> A           : (Possibly shrunken) matrix with regression 
	#                  coefficients
	# -> SigmaE      : (Possibly shrunken) error covariance matrix
	# -> threshold   : Signifies type of thresholding
	# -> absValueCut : Cut-off for elements selection based on absolute 
	#                  value thresholding (of A). Only when 
	#                  threshold='absValue'. Default = .25
	# -> FDRcut      : Cut-off for element selection based on local FDR
	#                  thresholding on test statistics.
	#                  Only when threshold = 'localFDR'. Default = .9
	# -> top         : element selection based on retainment 'top' number 
	#                  of absolute values of A.
	#                  Only when threshold = 'top'. Default = 10
	# -> statistics  : logical indicating if test statistics should be 
	#                  returned. Only when threshold = "localFDR". 
	#                  Default = FALSE
	# -> verbose     : Logical indicating if intermediate output should be 
	#                  printed on screen. Only when threshold='localFDR'. 
	#                  Default = TRUE
	#
	# DEPENDENCIES:
	# require("fdrtool")          # functions from package : fdrtool
	#
	# NOTES:
	# For the sparsification of the (inverse of the) SigmaE matrix, 
	# confer 'sparsifyS'
	# 
	########################################################################

	# input check
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
	if (nrow(A1) != ncol(A2)){ 
		stop("Dimensions of matrices A1 and A2 do not match.") 
	}	
	if (as.character(class(threshold)) != "character"){ 
		stop("Input (threshold) is of wrong class.") 
	}
	if (missing(threshold)){ 
		stop("Need to specify type of sparsification ('absValue' or 'localFDR' or 'top')") 
	}
	if (!(threshold %in% c("absValue", "localFDR", "top"))){ 
		stop("Input (threshold) should be one of {'absValue', 'localFDR', 'top'}") 
	}
	if (as.character(class(verbose)) != "logical"){ 
		stop("Input (verbose) is of wrong class.") 
	}

	# top elements of A
	if (threshold == "top"){
		# extra input checks for this option
		if (class(top) != "numeric"){ 
			stop("Input (top) is of wrong class") 
		} 
		if (length(top) != 2){ 
			stop("Input (top) must be of a numeric of length two") 
		}
		if (any(!.is.int(top))){ 
			stop("Ech element of input (top) should be a numeric integer") 
		}
		if (any(top <= 0)){ 
			stop("Each element of input (top) must be strictly positive") 
		} 
		if (any(top >= (ncol(A1))^2)){ 
			stop("Input (top) must be smaller than the number of elements of the input matrices A1 and A2") 
		}
		
		# determine threshold for top
		absValueCut <- c(sort(abs(A1), decreasing = TRUE)[ceiling(top[1])], 
					absValueCut1 <- sort(abs(A2), decreasing = TRUE)[ceiling(top[2])])
		threshold   <- "absValue"
	}

	# absolute value criterion
	if (threshold == "absValue"){
		# extra input checks for this option
		if (class(absValueCut) != "numeric"){ 
			stop("Input (absValueCut) is of wrong class") 
		} 
		if (length(absValueCut) != 2){ 
			stop("Input (absValueCut) must be a numeric of length two") 
		}  
		if (any(absValueCut <= 0)){ 
			stop("Input (absValueCut) must be positive") 
		}
        
		# selection
		nonzerosA1 <- which(abs(A1) >= absValueCut[1], arr.ind=TRUE)
		zerosA1 <- which(abs(A1) < absValueCut[1], arr.ind=TRUE)
		H0statA1 <- NULL
		nonzerosA2 <- which(abs(A2) >= absValueCut[2], arr.ind=TRUE)
		zerosA2 <- which(abs(A2) < absValueCut[2], arr.ind=TRUE)
		H0statA2 <- NULL
	}

	# local FDR values 
	if (threshold=="localFDR"){
		# extra input checks for this option
		if (as.character(class(SigmaE)) != "matrix"){ 
			stop("Input (SigmaE) is of wrong class.") 
		}
		if (!isSymmetric(SigmaE)){ 
			stop("Non-symmetric covariance matrix is provided.") 
		}
		if (!all(eigen(SigmaE, only.values=TRUE)$values > 0)){ 
			stop("Non positive-definite covariance matrix is provided.") 
		}
		if (nrow(A1) != nrow(SigmaE)){ 
			stop("Dimensions covariance matrix and A1 do not match.") 
		}
		if (as.character(class(FDRcut)) != "numeric"){ 
			stop("Input (FDRcut) is of wrong class.") 
		}
		if (length(FDRcut) != 2){ 
			stop("Input (FDRcut) is of wrong length.") 
		}
		if (any(is.na(FDRcut))){ 
			stop("Input (FDRcut) is not a positive number.") 
		}
		if (any(FDRcut <= 0)){ 
			stop("Input (FDRcut) is not a positive number.") 
		}
		if (any(FDRcut >= 1)){ 
			stop("Input (FDRcut) should be smaller than one.") 
		}
		if (as.character(class(statistics)) != "logical"){ 
			stop("Input (testStat) is of wrong class.") 
		}
	
		# calculate the variance of Y	
		# initiating co- and variances
		varMarg1 <- SigmaE
		varMarg2 <- SigmaE + A1 %*% SigmaE %*% t(A1)
		varMarg3 <- SigmaE + A1 %*% SigmaE %*% t(A1) + 
					A2 %*% SigmaE %*% t(A2) + 
					A1 %*% A1 %*% SigmaE %*% t(A1) %*% t(A1)         
		covMarg10 <- matrix(0, nrow(A1), ncol(A1))
		covMarg21 <- A1 %*% SigmaE
		covMargT_2andT_3 <- A1 %*% SigmaE
		covMargT_1andT_2 <- A1 %*% SigmaE + 
					A2 %*% SigmaE %*% t(A1) + 
					A1 %*% A1 %*% SigmaE %*% t(A1)         
		varMargT_1 <- varMarg3
		varMargT_2 <- varMarg2
		for (tau in 4:1000){        
			varMargT <- SigmaE + A1 %*% varMargT_1 %*% t(A1) + 
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
		        if (sum(abs(varMargT - varMargT_1)) < 10^(-10)){ break }
		}
		 
		# sparsify A1      
		# calculate "H0 statistics" for entries of A1: beta-hat / s.e.(beta-hat)
		H0statA1 <- as.numeric(A1 / sqrt(t(outer(diag(solve(varMargT)), diag(SigmaE)))))

		# local FDR procedure for A2
    		if (nrow(zerosA1) > 0){ 
			# local FDR procedure excluding already known zero entries of A1
			H0statA1[zerosA1] <- 0 
			lFDRs <- 1 - fdrtool(as.numeric(H0statA1[H0statA1 != 0]), 
						"normal", 
						cutoff.method="locfdr", 
						plot=verbose, 
						verbose=verbose)$lfdr
			nonzerosA1 <- H0statA1
			nonzerosA1[H0statA1 != 0] <- lFDRs
			zerosA1 <- which(nonzerosA1 <= FDRcut[1], arr.ind=TRUE)
			nonzerosA1 <- which(nonzerosA1 > FDRcut[1], arr.ind=TRUE)
		} else {
			# local FDR procedure without already known zero entries of A1
			lFDRs <- 1 - fdrtool(as.numeric(H0statA1), 
						"normal", 
						cutoff.method="locfdr", 
						plot=verbose, 
						verbose=verbose)$lfdr
			nonzerosA1 <- which(matrix(lFDRs, ncol=ncol(A1), byrow=FALSE) > FDRcut[1], arr.ind=TRUE)
			zerosA1 <- which(matrix(lFDRs, ncol=ncol(A1), byrow=FALSE) <= FDRcut[1], arr.ind=TRUE)
		}

		# sparsify A2
		# calculate "H0 statistics" for entries of A2: beta-hat / s.e.(beta-hat)
		H0statA2 <- as.numeric(A2 / sqrt(t(outer(diag(solve(varMargT)), diag(SigmaE)))))

		# local FDR procedure for A2
    		if (nrow(zerosA2) > 0){ 
			# local FDR procedure excluding already known zero entries of A2
			H0statA2[zerosA2] <- 0 
			lFDRs <- 1 - fdrtool(as.numeric(H0statA2[H0statA2 != 0]), 
						"normal", 
						cutoff.method="locfdr", 
						plot=verbose, 
						verbose=verbose)$lfdr
			nonzerosA2 <- H0statA2
			nonzerosA2[H0statA2 != 0] <- lFDRs
			zerosA2 <- which(nonzerosA2 <= FDRcut[2], arr.ind=TRUE)
			nonzerosA2 <- which(nonzerosA2 > FDRcut[2], arr.ind=TRUE)
		} else {
			# local FDR procedure without already known zero entries of A2
			lFDRs <- 1 - fdrtool(as.numeric(H0statA2), 
						"normal", 
						cutoff.method="locfdr", 
						plot=verbose, 
						verbose=verbose)$lfdr
			nonzerosA2 <- which(matrix(lFDRs, ncol=ncol(A2), byrow=FALSE) > FDRcut[2], arr.ind=TRUE)
			zerosA2 <- which(matrix(lFDRs, ncol=ncol(A2), byrow=FALSE) <= FDRcut[2], arr.ind=TRUE)
		}
	}    

	# report results of sparsification on screen?	
	if (verbose){
	        cat("-> Retained elements of A1: ", nrow(nonzerosA1), "\n")
	        cat("-> Corresponding to", round(nrow(nonzerosA1)/(nrow(nonzerosA1) +
			nrow(zerosA1)), 4) * 100, "% of possible A1 elements \n")
	        cat("-> Retained elements of A2: ", nrow(nonzerosA2), "\n")
	        cat("-> Corresponding to", round(nrow(nonzerosA2)/(nrow(nonzerosA2) + 
			nrow(zerosA2)), 4) * 100, "% of possible A2 elements \n")
	}
    
	# should the "H0 statistics" be exported too?
	if (!statistics){ 
		H0statA1 <- NULL 
		H0statA2 <- NULL 
	}

	# return the zero and nonzero elements of VAR(2) parameters A1 and A2
	return(list(nonzerosA1=nonzerosA1, 
			zerosA1=zerosA1, 
			statisticsA1=H0statA1, 
			nonzerosA2=nonzerosA2, 
			zerosA2=zerosA2, 
			statisticsA2=H0statA2))	    
}


