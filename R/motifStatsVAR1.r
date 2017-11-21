motifStatsVAR1 <- function(sparseA, verbose=TRUE){
	#######################################################################
	#
	# DESCRIPTION:
	# -> Function that detects standard network motifs from the lag-one 
	#    autoregression relationships.
	#
	# ARGUMENTS:
	# -> sparseA  : sparse autoregression matrix
	# -> verbose  : logical specifying whether summary output 
	#               should be displayed.
	# 
	# DEPENDENCIES:
	# ...         : ...
	#
	# REFERENCES:
	# -> Alon, U. (2007), "Network motifs: theory and experimental approaches", 
	#                     Nature Reviews Genetics, 8, 450-461.
	# 
	#######################################################################

	# Dependencies
	# require("base")
	# require("igraph")

	if (!is.matrix(sparseA)){
		stop("Input (sparseA) should be a matrix.")
	}
	if (nrow(sparseA) != ncol(sparseA)){
		stop("Input (sparseA) should be a square matrix.")
	}
  	if (all(sparseA != 0)){
		warning("Given input (sparseA) implies a saturated time-series chain graph.")
	}
	if (all(sparseA == 0)){
		warning("Given input (sparseA) implies an empty time-series chain graph.")
	}
	if (!is.logical(verbose)){
		stop("Input (verbose) should be a logical.")
	}

	# detection of self-regulators
	selfregulators <- list()
	idSelf <- which(diag(sparseA) != 0)
	if (length(idSelf) > 0){
		for (m in 1:length(idSelf)){
			selfregulator           <- matrix(c(idSelf[m], idSelf[m], 
                                                sign(diag(sparseA)[idSelf[m]])), nrow=1)
			colnames(selfregulator) <- c("t", "t+1", "sign")
			rownames(selfregulator) <- "path"
			selfregulators[[m]]     <- selfregulator
		}
	}

	# detection of feedback pairs
	feedbackpairs <- list()
	idSelfLag2 <- which(diag(adjacentMat(sparseA) %*% adjacentMat(sparseA)) != 0, arr.ind=TRUE)
	if (length(idSelfLag2) > 2){
		for (j1 in 1:(length(idSelfLag2)-1)){
			for (j2 in (j1+1):(length(idSelfLag2))){
				if (sparseA[idSelfLag2[j1], idSelfLag2[j2]] * sparseA[idSelfLag2[j2], idSelfLag2[j1]] != 0){
					feedbackpair <- rbind(c(idSelfLag2[j1], idSelfLag2[j2], 
                                            sign(sparseA[idSelfLag2[j2], idSelfLag2[j1]])),
						                  c(idSelfLag2[j2], idSelfLag2[j1], 
                                            sign(sparseA[idSelfLag2[j1], idSelfLag2[j2]])))
					colnames(feedbackpair) <- c("t", "t+1", "sign")
					rownames(feedbackpair) <- c("path1", "path2")
					feedbackpairs[[length(feedbackpairs) + 1]] <- feedbackpair
				}
			}
		}
	}

	# detection of bi-fans
	bifans <- list()
	idRegulators <- which(colSums(adjacentMat(sparseA)) != 0)
	if (length(idRegulators) > 1){
		for (j1 in 1:(length(idRegulators)-1)){
			for (j2 in (j1+1):(length(idRegulators))){
                		idSharedRegulatees <- intersect(which(sparseA[,j1] != 0), 
				                               which(sparseA[,j2] != 0))
				if (length(idSharedRegulatees) > 2){
                    			for (j3 in 1:(length(idSharedRegulatees)-1)){
                        			for (j4 in (j3+1):length(idSharedRegulatees)){
        						bifan <- rbind(c(idRegulators[j1], idSharedRegulatees[j3], 
                                                                         sign(sparseA[idSharedRegulatees[j3], 
										      idRegulators[j1]])),
						                       c(idRegulators[j1], idSharedRegulatees[j4], 
									  sign(sparseA[idSharedRegulatees[j4], 
                                                                                       idRegulators[j1]])),
							               c(idRegulators[j2], idSharedRegulatees[j3], 
                                             				 sign(sparseA[idSharedRegulatees[j3], 
                                                                                      idRegulators[j2]])),
								       c(idRegulators[j3], idSharedRegulatees[j4], 
						                         sign(sparseA[idSharedRegulatees[j4], 
								                      idRegulators[j2]])))       
							colnames(bifan) <- c("t", "t+1", "sign")
							rownames(bifan) <- c("path1", "path2", "path3", "path4")
							bifans[[length(bifans) + 1]] <- bifan
						}
					}
				}
			}
		}
	}

	# detection of diamonds
	diamonds <- list()
	idStartAndFinish <- which(adjacentMat(sparseA) %*% adjacentMat(sparseA) != 0, arr.ind=TRUE)
	if (nrow(idStartAndFinish) > 0){
		idIntermediates  <- lapply(1:nrow(idStartAndFinish), 
        	                   	   function(m, Z){ which(adjacentMat(sparseA)[,Z[m,2], drop=TRUE] * 
        	                                         adjacentMat(sparseA)[Z[m,1],, drop=TRUE] != 0) }, 
        	                   	   Z=idStartAndFinish)
		nPaths <- sapply(idIntermediates, function(Z){ length(Z) }, simplify=TRUE)
		if (any(nPaths > 1)){
			idDiamonds <- which(nPaths > 1)
			for (d in idDiamonds){
				for (j1 in 1:(length(idIntermediates[[d]])-1)){
					for (j2 in (j1+1):(length(idIntermediates[[d]]))){
						diamond <- rbind(c(idStartAndFinish[d,2], idIntermediates[[d]][j1], idStartAndFinish[d,1],
        	                               sign(sparseA[idIntermediates[[d]][j1], idStartAndFinish[d,2]] * 
        	                                    sparseA[idStartAndFinish[d,1], idIntermediates[[d]][j1]])),
								         c(idStartAndFinish[d,2], idIntermediates[[d]][j2], idStartAndFinish[d,1],
        	                               sign(sparseA[idIntermediates[[d]][j2], idStartAndFinish[d,2]] * 
        	                                    sparseA[idStartAndFinish[d,1], idIntermediates[[d]][j2]])))
						colnames(diamond) <- c("t", "t+1", "t+2", "sign")
						rownames(diamond) <- c("path1", "path2")
						diamonds[[length(diamonds)+1]] <- diamond
					}
				}
			}
		}
	}

	# detection of feedforwardloops
	feedforwardloops <- list()
	idFFloopStartAndFinish <- which((adjacentMat(sparseA) %*% adjacentMat(sparseA)) * adjacentMat(sparseA) != 0, arr.ind=TRUE)
	if (nrow(idFFloopStartAndFinish) > 0){
		for (d in 1:nrow(idFFloopStartAndFinish)){
			idIntermediates <- which(sparseA[,idFFloopStartAndFinish[d,1]] * sparseA[,idFFloopStartAndFinish[d,2]] != 0)
			# ensure the intermediate node does not coincide with the end node
			idSlh1 <- which(idIntermediates == idFFloopStartAndFinish[d,1])
			if (length(idSlh1) > 0){
				idIntermediates <- idIntermediates[-idSlh1]
			}
			# ensure the intermediate node does not coincide with the start node
			idSlh2 <- which(idIntermediates == idFFloopStartAndFinish[d,2])                
			if (length(idSlh1) > 0){
				idIntermediates <- idIntermediates[-idSlh2]
			}
			if (length(idIntermediates) > 0){
				for (j in 1:length(idIntermediates)){
					feedforwardloop <- rbind(c(idFFloopStartAndFinish[d,2], 
								   idFFloopStartAndFinish[d,1], 
                                                                   sign(sparseA[idFFloopStartAndFinish[d,1], 
                                                                                idFFloopStartAndFinish[d,2]])),
                                                                 c(idFFloopStartAndFinish[d,2], 
                                                                   idIntermediates[j],           
                                                                   sign(sparseA[idIntermediates[j], 
                                                                   idFFloopStartAndFinish[d,2]])),
                                                                 c(idIntermediates[j], 
                                                                   idFFloopStartAndFinish[d,1], 
                                                                   sign(sparseA[idFFloopStartAndFinish[d,1], 
                                                                   idIntermediates[j]])))                
					colnames(feedforwardloop) <- c("t", "t+1", "sign")
					rownames(feedforwardloop) <- c("path1", "path2", "path3")
					feedforwardloops[[length(feedforwardloops)+1]] <- feedforwardloop
				}
			}
		}
	}

	# detection of feedbackloops
	feedbackloops <- list()
	idFBloopStartAndFinish <- which((adjacentMat(sparseA) %*% adjacentMat(sparseA)) * t(adjacentMat(sparseA)) != 0, arr.ind=TRUE)
	if (nrow(idFBloopStartAndFinish) > 0){
		for (d in 1:nrow(idFBloopStartAndFinish)){
			idIntermediates <- which(sparseA[idFBloopStartAndFinish[d,1],] * sparseA[,idFBloopStartAndFinish[d,2]] != 0)
			# ensure the intermediate node does not coincide with the end node
			idSlh1 <- which(idIntermediates == idFBloopStartAndFinish[d,1])
			if (length(idSlh1) > 0){
				idIntermediates <- idIntermediates[-idSlh1]
			}
			# ensure the intermediate node does not coincide with the start node
			idSlh2 <- which(idIntermediates == idFBloopStartAndFinish[d,2])                
			if (length(idSlh1) > 0){
				idIntermediates <- idIntermediates[-idSlh2]
			}
			if (length(idIntermediates) > 0){
				for (j in 1:length(idIntermediates)){
					feedbackloop <- rbind(c(idFBloopStartAndFinish[d,1], 
						                idFBloopStartAndFinish[d,2], 
						                sign(sparseA[idFBloopStartAndFinish[d,2], 
						                             idFBloopStartAndFinish[d,1]])),
							       c(idFBloopStartAndFinish[d,2], 
						                 idIntermediates[j],          
						                 sign(sparseA[idIntermediates[j], 
						                              idFBloopStartAndFinish[d,2]])),
						                c(idIntermediates[j], 
						                  idFBloopStartAndFinish[d,1], 
						                  sign(sparseA[idFBloopStartAndFinish[d,1], 
						                               idIntermediates[j]])))                
					colnames(feedbackloop) <- c("t", "t+1", "sign")
					rownames(feedbackloop) <- c("path1", "path2", "path3")
					feedbackloops[[length(feedbackloops)+1]] <- feedbackloop
				}
			}
		}
	}

	# make motif table
	if (verbose){
		cat("motif summary table:", "\n")
		print(data.frame(motif=c("-> selfregulator    : ", 
                                 "-> feedback pair    : ", 
                                 "-> feedforward loop : ", 
                                 "-> feedback loop    : ", 
                                 "-> bifan            : ", 
                                 "-> diamond          : "), 
		                  frequency=c(length(selfregulators), 
                                      length(feedbackpairs), 
                                      length(feedforwardloops), 
                                      length(feedbackloops), 
                                      length(bifans), 
                                      length(diamonds))), row.names=FALSE)

		# specify grid
		hBetweenNodeWi  <- 4
		hBetweenMotifWi <- 3 * hBetweenNodeWi
		hCoord <- c(1, 
        	            1+1*hBetweenNodeWi, 
        	            1+1*hBetweenNodeWi+1*hBetweenMotifWi,
        	            1+2*hBetweenNodeWi+1*hBetweenMotifWi,
        	            1+2*hBetweenNodeWi+2*hBetweenMotifWi,
        	            1+3*hBetweenNodeWi+2*hBetweenMotifWi,
        	            1+4*hBetweenNodeWi+2*hBetweenMotifWi)
		vBetweenNodeWi  <- 3
		vBetweenMotifWi <- 3 * vBetweenNodeWi
		vCoord <- c(1, 
        	            1+1*vBetweenNodeWi, 
        	            1+2*vBetweenNodeWi,
        	            1+3*vBetweenNodeWi,  
        	            1+3*vBetweenNodeWi+1*vBetweenMotifWi,
        	            1+4*vBetweenNodeWi+1*vBetweenMotifWi,
        	            1+5*vBetweenNodeWi+1*vBetweenMotifWi)
		grid <- cbind(rep(hCoord, length(vCoord)), sort(rep(vCoord, length(hCoord))))
		idRemove <- c(intersect(which(grid[,2] == 1), which(grid[,1] < 10)),
			      intersect(which(grid[,2] > 15), which(grid[,1] > 40)),
			      intersect(intersect(which(grid[,2] > 15), which(grid[,1] < 30)), which(grid[,2] < 20)),
			      intersect(intersect(which(grid[,2] > 20), which(grid[,1] < 10)), which(grid[,2] < 24)))
		grid <- grid[-idRemove,]
	
		# specify adjacency matrix
		adjMat <- matrix(0, nrow=nrow(grid), ncol=nrow(grid))
		adjMat[33,34] <- 1; adjMat[35,30] <- 1;	adjMat[29,36] <- 1;
		adjMat[37,32] <- 1; adjMat[31,28] <- 1; adjMat[37,28] <- 1;
		adjMat[20,14] <- 1; adjMat[13, 7] <- 1;	adjMat[ 6,21] <- 1;
		adjMat[15,23] <- 1; adjMat[15, 2] <- 1;	adjMat[ 8,23] <- 1;
		adjMat[ 8, 2] <- 1; adjMat[24,18] <- 1; adjMat[24,11] <- 1;
		adjMat[18, 5] <- 1; adjMat[11, 5] <- 1;

		# specify node names
		nNames <- c(rep("D", 5), rep("C", 7), rep("B", 7), 
        	            rep("A", 7), rep("C", 2), rep("B", 4), rep("A", 6))


		# convert adjacency matrix to graph-object
		gObj <- graph.adjacency(adjMat, mode="directed")
		plot(gObj, layout=grid, vertex.color="black", vertex.label.color="white", 
		     vertex.label=nNames, edge.width=1, vertex.label.family="sans",
		     vertex.label.font=1.5, vertex.size=10, edge.arrow.size=0.5, edge.width=1)
		text(-1.0,  0.85, expression(italic(t)),   cex=1, font=3)
		text(-0.8,  0.85, expression(italic(t)+1), cex=1, font=3)
		text(-0.2,  0.60, expression(italic(t)),   cex=1, font=3)
		text( 0.0,  0.60, expression(italic(t)+1), cex=1, font=3)
		text( 0.6,  0.35, expression(italic(t)),   cex=1, font=3)
		text( 0.8,  0.35, expression(italic(t)+1), cex=1, font=3)
		text(-1.0, -0.90, expression(italic(t)),   cex=1, font=3)
		text(-0.8, -0.90, expression(italic(t)+1), cex=1, font=3)
		text(-0.2, -1.15, expression(italic(t)),   cex=1, font=3)
		text( 0.0, -1.15, expression(italic(t)+1), cex=1, font=3)
		text( 0.6, -1.15, expression(italic(t)),   cex=1, font=3)
		text( 0.8, -1.15, expression(italic(t)+1), cex=1, font=3)
		text(   1, -1.15, expression(italic(t)+2), cex=1, font=3)
		text(-1.0,  1.15, paste("# self-regulators: ",   length(selfregulators), sep=""))
		text(-0.1,  1.15, paste("# feeback-pairs: ",     length(feedbackpairs), sep=""))
		text( 0.7,  1.15, paste("# feedforward loops: ", length(feedforwardloops), sep=""))
		text(-1.0, -0.1,  paste("# feedback loops: ",    length(feedbackloops), sep=""))
		text(-0.1, -0.1,  paste("# bi-fans: ",           length(bifans), sep=""))
		text( 0.8, -0.1,  paste("# diamonds: ",          length(diamonds), sep=""))
	}

	# return the list of the various detected motif types
	return(list(selfregulators  =selfregulators, 
	            feedbackpairs   =feedbackpairs,
	            feedforwardloops=feedforwardloops,
	            feedbackloops   =feedbackloops,
	            bifans          =bifans,
                    diamonds        =diamonds))
}

