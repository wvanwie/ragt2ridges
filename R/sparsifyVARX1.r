sparsifyVARX1 <- function(X, 
                          A, 
                          B, 
                          SigmaE, 
                          threshold=c("absValue", "localFDR", "top"), 
                          absValueCut=rep(0.25, 2), 
                          FDRcut=rep(0.8, 2), 
                          top=rep(10, 2), 
                          zerosA=matrix(nrow=0, ncol=2), 
                          zerosB=matrix(nrow=0, ncol=2), 
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
	# -> B           : (Possibly shrunken) matrix with regression 
	#                  coefficients
	# -> SigmaE      : (Possibly shrunken) error covariance matrix
	# -> threshold   : Signifies type of thresholding
	# -> absValueCut : Cut-off for elements selection based on 
	#                  absolute value thresholding (of A). Only when 
	#                  threshold='absValue'. Default=0.25
	# -> FDRcut      : Cut-off for element selection based on local 
	#                  FDR thresholding on test statistics. Only when 
	#                  threshold='localFDR'. Default=0.9
	# -> top         : Element selection based on retainment 'top' 
	#                  number of absolute values of A. Only when 
	#                  threshold='top'. Default=10
	# -> statistics  : logical indicating if test statistics should be 
	#                  returned. Only when threshold = 'localFDR'. 
	#                  Default = FALSE.
	# -> verbose     : Logical indicating if intermediate output should 
	#                  be printed on screen. Only when 
	#                  threshold='localFDR'. Default = TRUE
	#
	# DEPENDENCIES:
	# require("fdrtool")          # functions from package : fdrtool
	# require("expm")             # functions from package : expm
	#
	# NOTES:
	# For the sparsification of the (inverse of the) SigmaE matrix, 
	# confer 'sparsifyS'
	# 
	########################################################################

	# input check
	if (as.character(class(A)) != "matrix"){ 
		stop("Input (A) is of wrong class.") 
	}
	if (nrow(A) != ncol(A)){ 
			stop("Matrix A is not square.") 
	}	
	if (as.character(class(B)) != "matrix"){ 
		stop("Input (B) is of wrong class.") 
	}
	if (nrow(A) != nrow(B)){ 
		stop("Dimensiona mismatch between matrices A and B.") 
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

	# top elements of A and B
	if (threshold == "top"){
		# extra input checks for this option
		if (class(top) != "numeric"){ 
			stop("Input (top) is of wrong class") 
		} 
		if (length(top) != 2){ 
			stop("Input (top) must be a scalar") 
		}
		if (any(!.is.int(top))){ 
			stop("Input (top) should be a numeric integer") 
		}
		if (any(top <= 0)){ 
			stop("Input (top) must be strictly positive") 
		} 
		if (top[1] >= prod(dim(A))){ 
			stop("Input (top) must be smaller than the number of elements of the input matrix A") 
		}
		if (top[2] >= prod(dim(B))){ 
			stop("Input (top) must be smaller than the number of elements of the input matrix B") 
		}
		
		# determine threshold for top
		A[zerosA] <- 0
		B[zerosB] <- 0
		absValueCut <- c(sort(abs(A), decreasing = TRUE)[ceiling(top[1])], 
					sort(abs(B), decreasing = TRUE)[ceiling(top[2])])
		threshold   <- "absValue"
	}

	# absolute value criterion
	if (threshold == "absValue"){
		# extra input checks for this option
		if (class(absValueCut) != "numeric"){ 
			stop("Input (absValueCut) is of wrong class") 
		} 
		if (length(absValueCut) != 2){ 
			stop("Input (absValueCut) must be a scalar") 
		}  
		if (any(absValueCut <= 0)){ 
			stop("Input (absValueCut) must be positive") 
		}
        
		# selection
		A[zerosA] <- 0
		B[zerosB] <- 0
		nonzerosA <- which(abs(A) >= absValueCut[1], arr.ind=TRUE)
		zerosA <- which(abs(A) < absValueCut[1], arr.ind=TRUE)
		H0statA <- NULL
		nonzerosB <- which(abs(B) >= absValueCut[2], arr.ind=TRUE)
		zerosB <- which(abs(B) < absValueCut[2], arr.ind=TRUE)
		H0statB <- NULL
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
		if (nrow(A) != nrow(SigmaE)){ 
			stop("Dimensions error covariance matrix and A do not match.") 
		}
		if (nrow(B) != nrow(SigmaE)){ 
			stop("Dimensions error covariance matrix and B do not match.") 
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
		Syy <- SigmaE
		for (tau in 1:1000) {
			Atau <- A %^% tau
			Syy <- Syy + Atau %*% SigmaE %*% t(Atau)
			if (max(abs(Atau)) < 10^(-20)){ break }
		}

		# sparsify A
		# calculate "H0 statistics" for entries of A
		H0statA <- A / sqrt(t(outer(diag(MASS::ginv(Syy)), diag(SigmaE))))

		# local FDR procedure for A
		if (nrow(zerosA) > 0){ 
			# local FDR procedure excluding already known zero entries of A
			H0statA[zerosA] <- 0 
			lFDRs <- 1 - fdrtool(as.numeric(H0statA[H0statA != 0]), 
						"normal", 
						cutoff.method="locfdr", 
						plot=verbose, 
						verbose=verbose)$lfdr
			nonzerosA <- H0statA
			nonzerosA[H0statA != 0] <- lFDRs
			zerosA <- which(nonzerosA <= FDRcut[1], arr.ind=TRUE)
			nonzerosA <- which(nonzerosA > FDRcut[1], arr.ind=TRUE)
		} else {
			# local FDR procedure without already known zero entries of B
			lFDRs <- 1 - fdrtool(as.numeric(H0statA), 
						"normal", 
						cutoff.method="locfdr", 
						plot=verbose, 
						verbose=verbose)$lfdr
			nonzerosA <- which(matrix(lFDRs, ncol=ncol(A), byrow=FALSE) > FDRcut[1], arr.ind=TRUE)
			zerosA <- which(matrix(lFDRs, ncol=ncol(A), byrow=FALSE) <= FDRcut[1], arr.ind=TRUE)
		}

		# sparsify B
		# calculate "H0 statistics" for entries of B
		Sxx <- .armaVAR1_VARYhat(X, TRUE, unbalanced=matrix(nrow=0, ncol=2))
		H0statB <- B / sqrt(t(outer(diag(MASS::ginv(Sxx)), diag(SigmaE))))

		# local FDR procedure for B
		if (nrow(zerosB) > 0){
			# local FDR procedure excluding already known zero entries of B
			H0statB[zerosB] <- 0 
			lFDRs <- 1 - fdrtool(as.numeric(H0statB[H0statB != 0]), 
						"normal", 
						cutoff.method="locfdr", 
						plot=verbose, 
						verbose=verbose)$lfdr
			nonzerosB <- H0statB
			nonzerosB[H0statB != 0] <- lFDRs
			zerosB <- which(nonzerosB <= FDRcut[2], arr.ind=TRUE)
			nonzerosB <- which(nonzerosB > FDRcut[2], arr.ind=TRUE)
		} else {
			# local FDR procedure without already known zero entries of B
			lFDRs <- 1 - fdrtool(as.numeric(H0statB), 
						"normal", 
						cutoff.method="locfdr", 
						plot=verbose, 
						verbose=verbose)$lfdr
			nonzerosB <- which(matrix(lFDRs, ncol=ncol(B), byrow=FALSE) > FDRcut[2], arr.ind=TRUE)
			zerosB <- which(matrix(lFDRs, ncol=ncol(B), byrow=FALSE) <= FDRcut[2], arr.ind=TRUE)
		}
	}    
	
	# report results of sparsification on screen?
	if (verbose){
		cat("-> Retained elements in A: ", nrow(nonzerosA), "\n")
		cat("-> Corresponding to", round(nrow(nonzerosA)/(nrow(nonzerosA) + 
			nrow(zerosA)), 4) * 100, "% of possible elements of A. \n")
		cat("-> Retained elements in B: ", nrow(nonzerosB), "\n")
		cat("-> Corresponding to", round(nrow(nonzerosB)/(nrow(nonzerosB) + 
			nrow(zerosB)), 4) * 100, "% of possible elements of B. \n")
	}
    
	# should the "H0 statistics" be exported too?
	if (!statistics){ 
		H0statA <- NULL
		H0statB <- NULL 
	}

	# return the zero and nonzero elements of VARX1 parameters A and B
	return(list(nonzerosA=nonzerosA, 
			zerosA=zerosA, 
			statisticsA=H0statA, 
			nonzerosB=nonzerosB, 
			zerosB=zerosB, 
			statisticsB=H0statB))	    
}


