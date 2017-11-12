ridgeVAR2 <- function(Y, 
                      lambdaA1=-1, 
                      lambdaA2=-1, 
                      lambdaP=-1, 
                      targetA1=matrix(0, dim(Y)[1], dim(Y)[1]), 
                      targetA2=matrix(0, dim(Y)[1], dim(Y)[1]), 
                      targetP=matrix(0, dim(Y)[1], dim(Y)[1]), 
                      targetPtype="none", fitA12="ml", 
                      zerosA1=matrix(nrow=0, ncol=2), 
                      zerosA2=matrix(nrow=0, ncol=2), 
                      zerosA1fit="sparse", 
                      zerosA2fit="sparse", 
                      zerosP=matrix(nrow=0, ncol=2), 
                      cliquesP=list(), 
                      separatorsP=list(), 
                      unbalanced=matrix(nrow=0, ncol=2), 
                      diagP=FALSE, 
                      efficient=TRUE, 
                      nInit=100, 
                      minSuccDiff=0.001){

	########################################################################
	# 
	# DESCRIPTION: 
	# Ridge estimation of the parameters of the VARX(1) model. The 
	# log-likelihood is augmented with a ridge penalty for all three 
	# parameters, A, the matrix of auto-regression coefficients, B, the 
	# matrix with regression coefficient of the time-varying covariates, 
	# and OmegaE, the inverse of the error variance. 
	# 
	# ARGUMENTS:
	# -> Y             : Three-dimensional array containing the data. The 
	#                    first, second and third dimensions correspond to 
	#                    covariates, time and samples, respectively. The 
	#                    data are assumed to centered covariate-wise.
	# -> lambdaA1      : Ridge penalty parameter to be used in the 
	#                    estimation of A1, the matrix with autro-regressive 
	#                    coefficients.
	# -> lambdaA2      : Ridge penalty parameter to be used in the 
	#                    estimation of A2, the matrix with regression 
	#                    coefficients.
	# -> lambdaP       : Ridge penalty parameter to be used in the
	#                    estimation of Omega, the precision matrix of the 
	#                    errors.
	# -> targetA1      : Target matrix to which the matrix A1 is to be 
	#                    shrunken.
	# -> targetA2      : Target matrix to which the matrix A2 is to be 
	#                    shrunken.
	# -> targetP       : Target matrix to which the precision matrix Omega
	#                    is to be shrunken.
	# -> zerosA1       : Matrix with indices of entries of A1 that are 
	#                    constrained to zero. The matrix comprises two 
	#                    columns, each row corresponding to an entry of A1.
	#                    The first column contains the row indices and the 
	#                    second the column indices.
	# -> zerosA2       : Matrix with indices of entries of A2 that are 
	#                    constrained to zero. The matrix comprises two 
	#                    columns, each row corresponding to an entry of A2.
	#                    The first column contains the row indices and the 
	#                    second the column indices.
	# -> zerosA1fit    : Character, either "sparse" or "dense". With 
	#                    "sparse", the matrix A1 is assumed to contain many
	#                    zeros and a computational efficient implementation 
	#                    of its estimation is employed. If "dense", it is 
	#                    assumed that A1 contains only few zeros and the 
	#                    estimation method is optimized computationally 
	#                    accordingly.
	# -> zerosA2fit    : Character, either "sparse" or "dense". With 
	#                    "sparse", the matrix A2 is assumed to contain many
	#                    zeros and a computational efficient implementation 
	#                    of its estimation is employed. If "dense", it is 
	#                    assumed that A2 contains only few zeros and the 
	#                    estimation method is optimized computationally 
	#                    accordingly.
	# -> zerosP        : A matrix with indices of entries of the precision 
	#                    matrix that are constrained to zero. The matrix 
	#                    comprises two columns, each row corresponding to an
	#                    entry of the adjacency matrix. The first column
	#                    contains the row indices and the second the column 
	#                    indices. The specified graph should be undirected 
	#                    and decomposable. If not, it is symmetrized and 
	#                    triangulated (unless cliquesP and seperatorsP are 
	#                    supplied). Hence, the employed zero structure may 
	#                    differ from the input 'zerosP'.
	# -> cliquesP      : A 'list'-object containing the node indices per 
	#                    clique as object from the 'rip-function.
	# -> separatorsP   : A 'list'-object containing the node indices per 
	#                    clique as object from the 'rip-function.
	# -> unbalanced    : A matrix with two columns, indicating the 
	#                    unbalances in the design. Each row represents a 
	#                    missing design point in the (time x individual)-
	#                    layout. The first and second column indicate the 
	#                    time and individual (respectively) specifics of 
	#                    the missing design point.
	# -> diagP         : Logical, indicates whether the error covariance 
	#                    matrix is assumed to be diagonal.
	# -> efficient     : Logical, affects estimation of A. Details below.
	# -> nInit         : Maximum number of iterations to used in maximum 
	#                    likelihood estimation.
	# -> minSuccDiff   : Minimum distance between estimates of two 
	#                    successive iterations to be achieved.
	# 
	# DEPENDENCIES:
	# library(rags2ridges)	    # functions: default.target, ridgeP, 
	#                                        ridgePchordal. Former two may 
	#                                        be called on the C++-side only.
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
	if (as.character(class(lambdaA1)) != "numeric"){ 
		stop("Input (lambdaA1) is of wrong class.") 
	}
	if (length(lambdaA1) != 1){ 
		stop("Input (lambdaA1) is of wrong length.") 
	}
	if (is.na(lambdaA1)){ 
		stop("Input (lambdaA1) is not a non-negative number.") 
	}
	if (lambdaA1 < 0){ 
		stop("Input (lambdaA1) is not a non-negative number.") 
	}
	if (as.character(class(lambdaA2)) != "numeric"){ 
		stop("Input (lambdaA2) is of wrong class.") 
	}
	if (length(lambdaA2) != 1){ 
		stop("Input (lambdaA2) is of wrong length.") 
	}
	if (is.na(lambdaA2)){ 
		stop("Input (lambdaA2) is not a non-negative number.") 
	}
	if (lambdaA2 < 0){ 
		stop("Input (lambdaA2) is not a non-negative number.") 
	}
	if (as.character(class(lambdaP)) != "numeric"){ 
		stop("Input (lambdaP) is of wrong class.") 
	}
	if (length(lambdaP) != 1){ 
		stop("Input (lambdaP) is of wrong length.") 
	}
	if (is.na(lambdaP)){ 
		stop("Input (lambdaP) is not a non-negative number.") 
	}
	if (lambdaP < 0){ 
		stop("Input (lambdaP) is not a non-negative number.") 
	}
	if (!is.null(unbalanced) & as.character(class(unbalanced)) != "matrix"){ 
		stop("Input (unbalanced) is of wrong class.") 
	}    
	if (!is.null(unbalanced)){ 
		if(ncol(unbalanced) != 2){ 
			stop("Wrong dimensions of the matrix unbalanced.") 
		} 
	} 
	if (as.character(class(zerosA1fit)) != "character"){ 
		stop("Input (zerosA1fit) is of wrong class.") 
	}
	if (as.character(class(zerosA1fit)) == "character"){ 
		if (!(zerosA1fit %in% c("dense", "sparse"))){ 
			stop("Input (zerosA1fit) ill-specified.") 
		} 
	}
	if (as.character(class(zerosA2fit)) != "character"){ 
		stop("Input (zerosA2fit) is of wrong class.") 
	}
	if (as.character(class(zerosA2fit)) == "character"){ 
		if (!(zerosA2fit %in% c("dense", "sparse"))){ 
			stop("Input (zerosA2fit) ill-specified.") 
		} 
	}
	if (as.character(class(diagP)) != "logical"){ 
		stop("Input (diagP) is of wrong class.") 
	}
	if (as.character(class(efficient)) != "logical"){ 
		stop("Input (efficient) is of wrong class.") 
	}
	if (as.character(class(nInit)) != "numeric" & as.character(class(nInit)) != "logical"){ 
		stop("Input (nInit) is of wrong class.") 
	}
	if (length(nInit) != 1){ 
		stop("Input (nInit) is of wrong length.") 
	}
	if (is.na(nInit)){ 
		stop("Input (nInit) is not a positive integer.") 
	}
	if (nInit < 0){ 
		stop("Input (nInit) is not a positive integer.") 
	}
	if (as.character(class(minSuccDiff)) != "numeric"){ 
		stop("Input (minSuccDiff) is of wrong class.") 
	}
	if (length(minSuccDiff) != 1){ 
		stop("Input (minSuccDiff) is of wrong length.") 
	}
	if (is.na(minSuccDiff)){ 
		stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (minSuccDiff <= 0){ 
		stop("Input (minSuccDiff) is not a positive number.") 
	}
	if (as.character(class(targetA1)) != "matrix"){ 
		stop("Input (targetA1) is of wrong class.") 
	}
	if (!is.null(targetA1)){ 
		if (dim(Y)[1] != nrow(targetA1)){ 
			stop("Dimensions of input (# rows of targetA1) do not match that of other input (# variates of Y).") 
		} 
	}
	if (!is.null(targetA1)){ 
		if (dim(Y)[1] != ncol(targetA1)){ 
			stop("Dimensions of input (# columns of targetA1) do not match that of other input (# variates of Y).") 
		} 
	}
	if (as.character(class(targetA2)) != "matrix"){ 
		stop("Input (targetA2) is of wrong class.") 
	}
	if (!is.null(targetA2)){ 
		if (dim(Y)[1] != nrow(targetA2)){ 
			stop("Dimensions of input (# rows of targetA2) do not match that of other input (# variates of Y).") 
		} 
	}
	if (!is.null(targetA2)){ 
		if (dim(Y)[1] != ncol(targetA2)){ 
			stop("Dimensions of input (# columns of targetA2) do not match that of other input (# variates of Y).") 
		} 
	}
	if (is.null(targetP)){ 
		targetP <- "Null" 
	}    
	if (as.character(class(targetP)) != "matrix"){ 
		stop("Input (targetP) is of wrong class.") 
	}    
	if (as.character(class(targetP)) == "matrix"){ 
		if(!isSymmetric(targetP)){ 
			stop("Non-symmetrical target for the precision matrix provided") 
		} 
	} 
	if (diagP & as.character(class(targetP)) == "matrix"){ 
		if(max(abs(targetP[upper.tri(targetP)])) != 0){ 
			stop("Inconsistent input (targetP v. diagP) provided") 
		} 
	}
 	if (as.character(class(targetPtype)) != "character"){ 
		stop("Input (targetPtype) of wrong class.") 
	} 
	if (targetPtype != "none"){ 
		if( length(intersect(targetPtype, c("DAIE", "DIAES", "DUPV", "DAPV", "DCPV", "DEPV", "Null"))) != 1 ){ 
			stop("Wrong default target for the precision matrix provided: see default.target for the options.") 
		} 
	} 
	if (as.character(class(targetP)) == "matrix"){ 
		if (dim(Y)[1] != nrow(targetP)){ 
			stop("Dimensions of input (targetP) do not match that of other input (Y).") 
		} 
	}
	if (!is.null(zerosA1) & as.character(class(zerosA1)) != "matrix"){ 
		stop("Input (zerosA1) is of wrong class.") 
	}    
	if (!is.null(zerosA1)){ 
		if(ncol(zerosA1) != 2){ 
			stop("Wrong dimensions of the (zerosA1) matrix.") 
		} 
	} 
	if (!is.null(zerosA1)){ 
		zerosA1 <- zerosA1[order(zerosA1[,2], zerosA1[,1]),] 
	}
	if (!is.null(zerosA2) & as.character(class(zerosA2)) != "matrix"){ 
		stop("Input (zerosA2) is of wrong class.") 
	}    
	if (!is.null(zerosA2)){ 
		if(ncol(zerosA2) != 2){ 
			stop("Wrong dimensions of the (zerosA2) matrix.") 
		} 
	} 
	if (!is.null(zerosA2)){ 
		zerosA2 <- zerosA2[order(zerosA2[,2], zerosA2[,1]),] 
	}
	if (!is.null(zerosP) & as.character(class(zerosP)) != "matrix"){ 
		stop("Input (zerosP) is of wrong class.") 
	}    
	if (!is.null(zerosP)){ 
		if(ncol(zerosP) != 2){ 
			stop("Wrong dimensions of the (zerosP).") 
		} 
	} 
	if (!is.null(zerosP)){ 
		zerosP <- zerosP[order(zerosP[,2], zerosP[,1]),] 
	}

	# targets only appears in a product with lambdas. 
	# moreover, the multiplication of a matrix times a scaler is faster in R.
	targetA1 <- lambdaA1 * targetA1;
	targetA2 <- lambdaA2 * targetA2;

	# estimation without support information neither on A nor on P
	if (nrow(zerosA1) == 0 && nrow(zerosA2) == 0 && nrow(zerosP) == 0){
		VAR2hat <- .armaVAR2_ridgeML(Y, 
	                                     lambdaA1, 
	                                     lambdaA2, 
	                                     lambdaP, 
	                                     targetA1, 
	                                     targetA2, 
	                                     targetP, 
	                                     targetPtype, 
	                                     fitA12, 
	                                     unbalanced, 
	                                     diagP, 
	                                     efficient, 
	                                     nInit, 
	                                     minSuccDiff);
	 	Phat <- VAR2hat$P; 
		A1hat <- VAR2hat$A[,1:dim(Y)[1]]; 
		A2hat <- VAR2hat$A[,-c(1:dim(Y)[1])];
	}

	# estimation with support information on A but not on P
	if ((nrow(zerosA1) > 0 | nrow(zerosA2) > 0) && nrow(zerosP) == 0){
		VAR2hat <- .armaVAR2_ridgeML_zerosA(Y, 
                                                    lambdaA1, 
                                                    lambdaA2, 
                                                    lambdaP, 
                                                    targetA1, 
                                                    targetA2, 
                                                    targetP, 
                                                    targetPtype, 
                                                    fitA12,
                                                    unbalanced, 
                                                    diagP, 
                                                    efficient, 
                                                    nInit, 
                                                    minSuccDiff, 
                                                    zerosA1[,1], 
                                                    zerosA1[,2], 
                                                    zerosA2[,1], 
                                                    zerosA2[,2], 
                                                    zerosA1fit, 
                                                    zerosA2fit);
		Phat <- VAR2hat$P;
		A1hat <- VAR2hat$A[,1:dim(Y)[1]]; 
		A2hat <- VAR2hat$A[,-c(1:dim(Y)[1])];
	}

	# estimation with support information both on A and on P
	if (nrow(zerosP) > 0){
		if (fitA12 == "ss"){
			# set profiles of missing (time, sample)-points to missing
			if (!is.null(unbalanced)){ 
				Y <- .armaVAR_array2cube_withMissing(Y, 
                                                                     unbalanced[,1], 
                                                                     unbalanced[,2]); 
			}

			# estimate A by SS minimization
			VARY <- .armaVAR2_VARYhat(Y, efficient);
			COVY <- .armaVAR2_COVYhat(Y);
                
			# estimate A
			if (nrow(zerosA1) == 0 && nrow(zerosA2) == 0){
				Ahat <- .armaVAR2_Ahat_ridgeSS(COVY, 
                                                               VARY, 
                                                               lambdaA1, 
                                                               lambdaA2, 
                                                               targetA1, 
                                                               targetA2)
			} else {
				# eigen-decomposition of VARY
				VARY <- .armaEigenDecomp(VARY)
				Ahat <- .armaVAR2_Ahat_zeros(diag(ncol(targetA1)), 
                                                             COVY, 
                                                             VARY$vectors, 
                                                             VARY$values, 
                                                             lambdaA1, 
                                                             lambdaA2, 
                                                             targetA1, 
                                                             targetA2, 
                                                             fitA12, 
                                                             zerosA1[,1], 
                                                             zerosA1[,2], 
                                                             zerosA2[,1], 
                                                             zerosA2[,2], 
                                                             zerosA1fit, 
                                                             zerosA2fit)
			}

			# calculate Se
			Se <- .armaVAR2_Shat_ML(Y, 
                                                Ahat[,1:dim(Y)[1]], 
                                                Ahat[,-c(1:dim(Y)[1])]);
	
			# if cliques and separators of support of P are not provided:
			if (length(cliquesP)==0){
				supportPinfo <- support4ridgeP(zeros=zerosP, nNodes=dim(Y)[1]);
				cliquesP     <- supportPinfo$cliques; 
				separatorsP  <- supportPinfo$separators; 
				zerosP       <- supportPinfo$zeros;
			}
	
			# ridge ML estimation of Se
			if (is.character(targetP)){ 
				target <- .armaP_defaultTarget(Se, 
                                                               targetType=targetPtype, 
                                                               fraction=0.0001, 
                                                               multiplier=0) 
			} else { 
				target <- targetP 
			}
			Phat <- ridgePchordal(Se, 
                                              lambda=lambdaP, 
                                              target=target, 
                                              zeros=zerosP, 
                                              cliques=cliquesP, 
                                              separators=separatorsP, 
                                              type="Alt", 
                                              verbose=FALSE)

		}
		if (fitA12 == "ml"){
            		# set profiles of missing (time, sample)-points to missing
			if (!is.null(unbalanced)){ 
				Y <- .armaVAR_array2cube_withMissing(Y, 
                                                                     unbalanced[,1], 
                                                                     unbalanced[,2]); 
			}

			# estimate A by SS minimization
			VARY <- .armaVAR2_VARYhat(Y, efficient);
			COVY <- .armaVAR2_COVYhat(Y);
			Ahat <- .armaVAR2_Ahat_ridgeSS(COVY, 
                                                       VARY, 
                                                       lambdaA1, 
                                                       lambdaA2, 
                                                       targetA1, 
                                                       targetA2)

			# calculate Se
			Se <- .armaVAR2_Shat_ML(Y, 
                                                Ahat[,1:dim(Y)[1]], 
                                                Ahat[,-c(1:dim(Y)[1])]);	
            
			# if cliques and separators of support of P are not provided:
			if (length(cliquesP)==0){
				supportPinfo <- support4ridgeP(zeros=zerosP, 
                                                               nNodes=dim(Y)[1]);
				cliquesP     <- supportPinfo$cliques; 
				separatorsP  <- supportPinfo$separators; 
				zerosP       <- supportPinfo$zeros;
			}
	
			# ridge ML estimation of Se
			if (is.character(targetP)){ 
				target <- .armaP_defaultTarget(Se, 
                                                               targetType=targetPtype, 
                                                               fraction=0.0001, 
                                                               multiplier=0) 
			} else { 
				target <- targetP 
			}
			Phat <- ridgePchordal(Se, 
                                              lambda=lambdaP, 
                                              target=target, 
                                              zeros=zerosP, 
                                              cliques=cliquesP, 
                                              separators=separatorsP, 
                                              type="Alt", 
                                              verbose=FALSE)

			########################################################
			# estimate parameters by ML, using the SS estimates as initials
			########################################################
            
			# eigen-decomposition of VARY               
			VARY <- .armaEigenDecomp(VARY)

			for (u in 1:nInit){
				# store latest estimates
				Aprev <- Ahat; 
                                Pprev <- Phat;

				# estimate A
				if (nrow(zerosA1) == 0 && nrow(zerosA2) == 0){
					Ahat <- .armaVAR2_Ahat_ridgeML(Phat, 
                                                                       COVY, 
                                                                       VARY$vectors, 
                                                                       VARY$values, 
                                                                       lambdaA1, 
                                                                       lambdaA2, 
                                                                       targetA1, 
                                                                       targetA2)
				} else {
					Ahat <- .armaVAR2_Ahat_zeros(Phat, 
                                                                     COVY, 
                                                                     VARY$vectors, 
                                                                     VARY$values, 
                                                                     lambdaA1, 
                                                                     lambdaA2, 
                                                                     targetA2, 
                                                                     targetA2, 
                                                                     fitA12, 
                                                                     zerosA1[,1], 
                                                                     zerosA1[,2], 
                                                                     zerosA2[,1], 
                                                                     zerosA2[,2], 
                                                                     zerosA1fit, 
                                                                     zerosA2fit)
				}

				# calculate Se
				Se <- .armaVAR2_Shat_ML(Y, 
                                                        Ahat[,1:dim(Y)[1]], 
                                                        Ahat[,-c(1:dim(Y)[1])]);
                
				# ridge ML estimation of Se
				if (is.character(targetP)){ 
					target <- .armaP_defaultTarget(Se, 
                                                                       targetType=targetPtype, 
                                                                       fraction=0.0001, 
                                                                       multiplier=0) 
				} else { 
					target <- targetP 
				}
				Phat <- ridgePchordal(Se, 
                                                      lambda=lambdaP, 
                                                      target=target, 
                                                      zeros=zerosP, 
                                                      cliques=cliquesP, 
                                                      separators=separatorsP, 
                                                      type="Alt", 
                                                      verbose=FALSE)
		
				# assess convergence
				if (.armaVAR2_convergenceEvaluation(Ahat, Aprev, Phat, Pprev) < minSuccDiff){ 
					break 
				}
			}
		}
		A1hat <- Ahat[,1:dim(Y)[1]]; 
                A2hat <- Ahat[,-c(1:dim(Y)[1])];    	
	}
	return(list(A1=A1hat, 
                    A2=A2hat, 
                    P=Phat, 
                    lambdaA1=lambdaA1, 
                    lambdaA2=lambdaA2, 
                    lambdaP=lambdaP))
}


