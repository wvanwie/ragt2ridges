nodeStatsVAR2 <- function(sparseA1, 
                          sparseA2, 
                          sparseP, 
                          as.table=FALSE){

	#######################################################################
	#
	# DESCRIPTION:
	# -> Function that calculates various network statistics from a 
	#    sparse VAR(1) model.
	#
	# ARGUMENTS:
	# -> sparseA  : sparse regression coefficient matrix
	# -> sparseP  : sparse precision/partial correlation matrix
	# -> as.table : logical indicating if output should be returned 
	#               as table; default = FALSE.
	# 
	# NOTES (network statistics produced):
	# -> degreeA1in        	 : number of one-lag (temporal) edges 
	#                          pointing to each node ('in'-degree).
	# -> degreeA1out       	 : number of one-lag (temporal) edges 
	#                          leaving each node ('out'-degree).
	# -> nNegA1in          	 : number of negative one-lag (temporal) 
	#                          edges pointing to each node.
	# -> nPosA1in          	 : number of positive one-lag (temporal)
	#                          edges pointing to each node ('in'-degree)
	# -> nNegA1out         	 : number of negative one-lag (temporal)
	#                          edges leaving each node ('out'-degree)
	# -> nPosA1out         	 : number of positive one-lag (temporal)
	#                          edges leaving each node ('out'-degree)
	# -> degreeA2in        	 : number of two-lag (temporal) edges 
	#                          pointing to each node ('in'-degree).
	# -> degreeA2out       	 : number of two-lag (temporal) edges 
	#                          leaving each node ('out'-degree).
	# -> nNegA2in          	 : number of negative two-lag (temporal) 
	#                          edges pointing to each node.
	# -> nPosA2in          	 : number of positive two-lag (temporal)
	#                          edges pointing to each node ('in'-degree)
	# -> nNegA2out         	 : number of negative two-lag (temporal)
	#                          edges leaving each node ('out'-degree)
	# -> nPosA2out         	 : number of positive two-lag (temporal)
	#                          edges leaving each node ('out'-degree)
	# -> degreePe          	 : number of contemporaneous edges of 
	#                          each node (as implied by the error 
	#                          precision matrix)
	# -> betweennessPe     	 : vector representing the contemporaneous 
	#                          betweenness centrality for each node.
	# -> closenessPe       	 : vector representing the contemporaneous 
	#                          closeness centrality for each node.
	# -> eigenCentralityPe 	 : vector representing the contemporaneous 
	#                          eigen centrality for each node.
	# -> nNegPe            	 : vector representing the number of negative 
	#                          contemporaneous edges for each node.
	# -> nPosPe            	 : vector representing the number of positive 
	#                          contemporaneous edges for each node.
	# -> variancePe        	 : vector representing the error variance of 
	#                          each node.
	# -> partialVarPe      	 : vector representing the partial error 
	#                          variance of each node.
	# -> varianceY		 : vector representing the variance of 
	#                          each node.
	# -> degreePy		 : number of edges of each node in the 
	#                          global Markov graph.
	# -> betweennessPy	 : vector representing the betweenness 
	#                          centrality for each node in the global 
	#                          Markov graph.
	# -> closenessPy	 : vector representing the closeness 
	#                          centrality for each node in the global 
	#                          Markov graph.
	# -> eigenCentralityPy	 : vector representing the eigen centrality 
	#                          for each node in the global Markov graph.
	# -> mutualInfo_Tplus1	 : vector with for each node its mutual 
	#                          information with all other nodes at the 
	#                          next (t+1) time point.
	# -> mutualInfo_Tplus2	 : vector with for each node its mutual 
	#                          information with all other nodes at the 
	#                          (t+2)-th time point.
	# -> itemResponse_Tplus1 : vector with for each node its  mean absolute 
	#                          impulse response on all other nodes at the 
	#                          next (t+1) time point.
	# -> itemResponse_Tplus2 : vector with for each node its  mean 
	#                          absolute impulse response on all other 
	#                          nodes at the (t+2)-th time point.
	# - Future versions of this function may include additional statistics
	# 
	# DEPENDENCIES:
	# require("igraph")      : functions from package : 
	#                          graph.adjacency, degree, closeness, 
	#                          betweenness, evcent
	#
	# REFERENCES:
	# -> Newman, M.E.J. (2010), "Networks: an introduction", 
	#                            Oxford University Press
	# 
	#######################################################################

	# Dependencies
	# require("base")
	# require("igraph")

	if (!is.matrix(sparseA1)){
		stop("Input (sparseA1) should be a matrix")
	}
	if (!is.matrix(sparseA2)){
		stop("Input (sparseA2) should be a matrix")
	}
	if (!is.matrix(sparseP)){
		stop("Input (sparseP) should be a matrix")
	}
	if (!isSymmetric(sparseP)){
		stop("Input (sparseP) should be a symmetric matrix")
	}
	if (!evaluateS(sparseP, verbose = FALSE)$posEigen){
		stop("Input (sparseP) is expected to be positive definite")
	}
	if (class(as.table) != "logical"){
		stop("Input (as.table) is of wrong class")
	} 
	
	# some warnings
	if (all(sparseA1 != 0)){
		warning("Given input (sparseA1) implies a saturated conditional independence graph")
	}
	if (all(sparseA2 != 0)){
		warning("Given input (sparseA2) implies a saturated conditional independence graph")
	}
	if (all(sparseA1 == 0) & all(sparseA2 == 0)){
		warning("Given inputs (sparseA1, sparseA2) imply an empty conditional independence graph")
	}
	if (all(sparseP != 0)){
		warning("Given input (sparseP) implies a saturated contemporaneous conditional independence graph")
	}
	if (all(sparseP[!diag(nrow(sparseP))] == 0)){
		warning("Given input (sparseP) implies an empty contemporaneous conditional independence graph")
	}

	###############################################
        # statistics from A1
        ###############################################
	
       	# in and out degree        
       	degreeA1out <- ncol(sparseA1) - colSums(sparseA1==0)
       	degreeA1in  <- nrow(sparseA1) - rowSums(sparseA1==0)
        
	# signs of edges
	nPosA1out <- apply(sign(sparseA1), 2, function(Z){ sum(Z == 1) }) 
	nNegA1out <- apply(sign(sparseA1), 2, function(Z){ sum(Z == -1) })

	# signs of edges
	nPosA1in <- apply(sign(sparseA1), 1, function(Z){ sum(Z == 1) }) 
	nNegA1in <- apply(sign(sparseA1), 1, function(Z){ sum(Z == -1) })

	# centrality measures of A1
    	# to be included


	###############################################
       	# statistics from A2
       	###############################################
	
       	# in and out degree        
       	degreeA2out <- ncol(sparseA2) - colSums(sparseA2==0)
       	degreeA2in  <- nrow(sparseA2) - rowSums(sparseA2==0)
        
	# signs of edges
	nPosA2out <- apply(sign(sparseA2), 2, function(Z){ sum(Z == 1) }) 
	nNegA2out <- apply(sign(sparseA2), 2, function(Z){ sum(Z == -1) })

	# signs of edges
	nPosA2in <- apply(sign(sparseA2), 1, function(Z){ sum(Z == 1) }) 
	nNegA2in <- apply(sign(sparseA2), 1, function(Z){ sum(Z == -1) })

	# centrality measures of A2
    	# to be included


    	############################################### 
    	# statistics from Pe
       	###############################################

       	# (partial) variance of the error
       	pvarsPe <- 1/diag(sparseP)
       	Se      <- solve(sparseP)
       	varsPe  <- diag(Se)

	# signs of edges of P
	slh <- diag(sparseP)
	diag(sparseP) <- 0 
	nPosPe <- apply(sign(sparseP), 2, function(Z){ sum(Z == 1) }) 
	nNegPe <- apply(sign(sparseP), 2, function(Z){ sum(Z == -1) })
       	diag(sparseP) <- slh

    	# adjacency to graphical object
	adjMat  <- adjacentMat(sparseP)
	CIGerror <- graph.adjacency(adjMat, mode = "undirected")
	    
	# centrality measures of P
	degreePe          <- degree(CIGerror)
	betweennessPe     <- betweenness(CIGerror)
	closenessPe       <- closeness(CIGerror)
	eigenCentralityPe <- evcent(CIGerror)$vector


	###############################################
	# statistics for Y
	###############################################

	# centrality measures of Py
	CIGy <- graph.adjacency(CIGofVAR2(sparseA1, 
                                          sparseA2, 
                                          sparseP, 
                                          "global"), 
                                mode="undirected")
	degreePy          <- degree(CIGy)
	betweennessPy     <- betweenness(CIGy)
	closenessPy       <- closeness(CIGy)
	eigenCentralityPy <- evcent(CIGy)$vector
        		

	# calculate the variance of Y	
	# initiating co- and variances
	varMarg1 <- Se
	varMarg2 <- Se + sparseA1 %*% Se %*% t(sparseA1)
	varMarg3 <- Se + sparseA1 %*% Se %*% t(sparseA1) + 
	            sparseA2 %*% Se %*% t(sparseA2) + 
	            sparseA1 %*% sparseA1 %*% Se %*% t(sparseA1) %*% t(sparseA1)         
	covMarg10 <- matrix(0, nrow(sparseA1), ncol(sparseA1))
	covMarg21 <- sparseA1 %*% Se
	covMargT_2andT_3 <- sparseA1 %*% Se
	covMargT_1andT_2 <- sparseA1 %*% Se + 
	                    sparseA2 %*% Se %*% t(sparseA1) + 
	                    sparseA1 %*% sparseA1 %*% Se %*% t(sparseA1)         
	varMargT_1 <- varMarg3
	varMargT_2 <- varMarg2
	for (tau in 4:1000){        
		varMargT <- Se + sparseA1 %*% varMargT_1 %*% t(sparseA1) + 
		            sparseA2 %*% varMargT_2 %*% t(sparseA2) + 
		            sparseA1 %*% covMargT_1andT_2 %*% t(sparseA2) + 
	                    t(sparseA1 %*% covMargT_1andT_2 %*% t(sparseA2))
		covMargTandT_1 <- sparseA1 %*% varMargT_1 + 
		                  sparseA2 %*% varMargT_1 %*% t(sparseA1) + 
		                  sparseA2 %*% covMargT_2andT_3 %*% t(sparseA2)
		varMargT_2 <- varMargT_1
		varMargT_1 <- varMargT
		covMargT_2andT_3 <- covMargT_1andT_2
		covMargT_1andT_2 <- covMargTandT_1
	        if (sum(abs(varMargT - varMargT_1)) < 10^(-10)){ break }
	}
	varsY <- diag(varMargT)

	# Calculate nodes' mutual information
	MIatTplus2         <- mutualInfoVAR2(sparseA1, sparseA2, sparseP, 2) 
	names(MIatTplus2)  <- colnames(sparseP)
	MIatTplus3         <- mutualInfoVAR2(sparseA1, sparseA2, sparseP, 3)  
	names(MIatTplus3)  <- colnames(sparseP)

	# Calculate nodes' influence response
	IRFatTplus2        <- apply(abs(impulseResponseVAR2(sparseA1, sparseA2, 2)), 2, mean)
	names(MIatTplus2)  <- colnames(sparseP)
	IRFatTplus3        <- apply(abs(impulseResponseVAR2(sparseA1, sparseA2, 3)), 2, mean)
	names(IRFatTplus3) <- colnames(sparseP)

	# return
	if (as.table){
		networkStats <- cbind(degreeA1in, 
					degreeA1out, 
					nNegA1in, 
					nPosA1in, 
					nNegA1out, 
					nPosA1out, 
					degreeA2in, 
					degreeA2out, 
					nNegA2in, 
					nPosA2in, 
					nNegA2out, 
					nPosA2out, 
					degreePe, 
					betweennessPe, 
					closenessPe, 
					eigenCentralityPe,
					nNegPe, 
					nPosPe, 
					varsPe, 
					pvarsPe, 
					varsY, 
					degreePy, 
					betweennessPe, 
					closenessPy, 
					eigenCentralityPy,
					MIatTplus2, 
					MIatTplus3, 
					IRFatTplus2, 
					IRFatTplus3)
		colnames(networkStats) <- c("degreeA1in", 
						"degreeA1out", 
						"nNegA1in", 
						"nPosA1in", 
						"nNegA1out", 
						"nPosA1out", 
						"degreeA2in", 
						"degreeA2out", 
						"nNegA2in", 
						"nPosA2in", 
						"nNegA2out", 
						"nPosA2out",
						"degreePe", 
						"betweennessPe", 
						"closenessPe", 
						"eigenCentralityPe", 
						"nNegPe", 
						"nPosPe", 
						"variancePe", 
						"partialVarPe", 
						"varianceY", 
						"degreePy", 
						"betweennessPy", 
						"closenessPy", 
						"eigenCentralityPy",
						"mutualInfo_Tplus2", 
						"mutualInfo_Tplus3", 
						"itemResponse_Tplus2", 
						"itemResponse_Tplus3")
		return(networkStats)
	} 
	if (!as.table){
		return(list(degreeA1in=degreeA1in, 
				degreeA1out=degreeA1out, 
				nNegA1in=nNegA1in, 
				nPosA1in=nPosA1in, 
				nNegA1out=nNegA1out, 
				nPosA1out=nPosA1out, 
				degreeA2in=degreeA2in, 
				degreeA2out=degreeA2out, 
				nNegA2in=nNegA2in, 
				nPosA2in=nPosA2in, 
				nNegA2out=nNegA2out, 
				nPosA2out=nPosA2out, 
				degreePe=degreePe, 
				betweennessPe=betweennessPe, 
				closenessPe=closenessPe, 
				eigenCentralityPe=eigenCentralityPe, 
				nNegPe=nNegPe, 
				nPosPe=nPosPe, 
				variancePe=varsPe, 
				partialVarPe=pvarsPe, 
				varianceY=varsY, 
				degreePy=degreePy, 
				betweennessPy=betweennessPy, 
				closenessPy=closenessPy, 
				eigenCentralityPy=eigenCentralityPy,
				mutualInfo_Tplus2=MIatTplus2, 
				mutualInfo_Tplus3=MIatTplus3, 
				itemResponse_Tplus2=IRFatTplus2, 
				itemResponse_Tplus3=IRFatTplus3))
	}
}




