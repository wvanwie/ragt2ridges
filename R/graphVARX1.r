graphVARX1 <- function(sparseA, 
			sparseB, 
			sparseP, 
			type="TSCG", 
			side="right", 
			prune=TRUE, 
			nNamesY=NULL, 
			nNamesX=NULL, 
			main=NULL, 
			vertex.color.X="lightcyan2", 
			vertex.color.T0="lightcyan2", 
			vertex.color.T1="lightcyan2", 
			vertex.frame.color="steelblue", 
			vertex.label.cex=-1, 
			vertex.label.color.X="black", 
			vertex.label.color.T0="black", 
			vertex.label.color.T1="black", 
			vertex.label.font=1.5, 
			vertex.size=-1, 
			edge.arrow.size=-1, 
			edge.width=-1, 
                        ...){

	#######################################################################
	#
	# DESCRIPTION:
	# Plot temporal relations as implied by A and B (the matrices with lag 
	# one and two regression coefficient of the VAR(2)) model.
	#
	# ARGUMENTS:
	#
	# -> sparseA                 : Matrix A1 of lag one regression 
	#                              parameters, which is assumed to be sparse.
	# -> sparseB                 : Matrix A2 of lag two regression 
	#                              parameters, which is assumed to be sparse.
	# -> sparseP                 : Matrix P of precision of the error, which 
	#                              is assumed to be sparse.
	# -> type                    : A 'character' indicating what should be 
	#                              plotted. If type='TSCG'. the time series
	#                              chain graph is plotted, while type 
	#                              type='Aonly' limits this graph
	#                              to the temporal relations. If 
	#                              type='globalPC' or type='contempPC', the 
	#                              global or contemporaneous (respectively) 
	#                              partial correlation graph is plotted.
	# -> side                    : A 'character' indicating wether the 
	#                              contemporaneous dependencies should be 
	#                              plotted on the 'left' (time t) or the 
	#                              'right' (time t+1) side. Only active 
	#                              when type='TSCG'.
	# -> prune                   : A 'logical' indicating whether to remove 
	#                              covariates without any temporal relations 
	#                              (as implied by 'sparseA', 'sparseB', 
	#                              'sparseP').
	# -> nNamesY                  : A 'character' containing the Y variate 
	#                              names to written inside the nodes.
	# -> nNamesX                  : A 'character' containing the X covariate 
	#                              names to written inside the nodes.
	# -> main                    : The 'character' to be plotted as title
	#                              above the graph.
	# -> vertex.color.X          : Color of nodes at time point t.
	# -> vertex.color.T0         : Color of nodes at time point t+1. Ignored 
	#                              when type='globalPC' or type='contempPC'.
	# -> vertex.color.T1         : Color of nodes at time point t+2. Ignored 
	#                              when type='globalPC' or type='contempPC'.
	# -> vertex.frame.color      : Refer to 'plot.igraph'.
	# -> vertex.label.cex        : Refer to 'plot.igraph'.
	# -> vertex.label.color.X    : Color of the node label at time point t.
	# -> vertex.label.color.T0   : Color of the node label at time point t+1. 
	#                              Ignored when type='globalPC' or 
	#                              type='contempPC'.
	# -> vertex.label.color.T1   : Color of the node label at time point t+2. 
	#                              Ignored when type='globalPC' or
	#                              type='contempPC'.
	# -> vertex.label.font       : Refer to 'plot.igraph'.
	# -> vertex.size             : Refer to 'plot.igraph'.
	# -> edge.arrow.size         : Refer to 'plot.igraph'.
	# -> edge.width              : Refer to 'plot.igraph'.
	# -> ...                     : Other arguments to be passed to 
	#                              'plot.igraph'.
	# 
	# DEPENDENCIES:
	# library(igraph)	    # functions: delete.edges, E, ecount, 
	#                             plot.igraph, graph.adjacency, 
	#                             maximum.cardinality.search
	#
	# NOTES:
   	# - igraph does not support visualization of mixed graphs, including both 
	#   directed and undirected edges. As a consequence the time-series chain 
	#   graph is undirected: the edges from t to t+1 should be thought of 
	#   as directed. 
	# - if "white" labels, color combination "navy" and "blue" looks nice. 
	# - if "black" labels, color combination "tomato1" and "red" looks nice. 
	#     	
	#######################################################################

	# input check
	if (as.character(class(sparseA)) != "matrix"){ 
		stop("Input (sparseA) is of wrong class.") 
	}
	if (nrow(sparseA) != ncol(sparseA)){ 
		stop("Matrix sparseA is not square.") 
	}
	if (as.character(class(sparseB)) != "matrix"){ 
		stop("Input (sparseB) is of wrong class.") 
	}
	if (nrow(sparseB) != ncol(sparseA)){ 
		stop("Number of rows of matrix sparseB does not match with that of sparseA.") 
	}
	if (as.character(class(sparseP)) != "matrix"){ 
		stop("Input (sparseP) is of wrong class.") 
	}
	if (nrow(sparseP) != ncol(sparseP)){ 
		stop("Matrix sparseP is not square.") 
	}
	if (as.character(class(prune)) != "logical"){ 
		stop("Input (prune) is of wrong class.") 
	}
	if (as.character(class(vertex.size)) != "numeric"){ 
		stop("Input (vertex.size) is of wrong class.") 
	}
	if (as.character(class(vertex.color.X)) != "character"){ 
		stop("Input (vertex.color.X) is of wrong class.") 
	}
	if (as.character(class(vertex.color.T0)) != "character"){ 
		stop("Input (vertex.color.T0) is of wrong class.") 
	}
	if (as.character(class(vertex.label.color.X)) != "character"){ 
		stop("Input (vertex.label.color.X) is of wrong class.") 
	}
	if (as.character(class(vertex.label.color.T0)) != "character"){ 
		stop("Input (vertex.label.color.T0) is of wrong class.") 
	}

	# if no covariate names are specified the columns and 
	# row names of A and B and P are given names 1, 2, et cetera
	if (is.null(nNamesY)){ 
		nNamesY <- as.character(1:nrow(sparseA)) 
	}
	if (is.null(nNamesX)){ 
		nNamesX <- as.character(1:ncol(sparseB)) 
	}

	if (type=="TSCG"){     
		# prune covariate without connections (according to sparseA, sparseB and sparseP)
		if (prune){ 
			idRemoveA <- intersect(which(apply(sparseA, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA, 2, function(Z){ all(Z == 0) })))
			idRemoveBy <- which(apply(sparseB, 1, function(Z){ all(Z == 0) }))
			idRemoveX <- which(apply(sparseB, 2, function(Z){ all(Z == 0) }))
			diag(sparseP) <- 0 
			idRemoveP <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			idRemoveY <- intersect(intersect(idRemoveA, idRemoveBy), idRemoveP)
			if (length(idRemoveY) > 0){ 
				sparseA <- sparseA[-idRemoveY, -idRemoveY]
				sparseB <- sparseB[-idRemoveY, ]
				sparseP <- sparseP[-idRemoveY, -idRemoveY] 
				nNamesY <- nNamesY[-idRemoveY]
			}
			if (length(idRemoveX) > 0){ 
				sparseB <- sparseB[,-idRemoveX]
				nNamesX <- nNamesX[-idRemoveX]
			}
		}
		
		# store number of nodes of A and B
		nNodesY <- nrow(sparseA)
		nNodesX <- ncol(sparseB)

		# default node size and font size if not provided
		if (vertex.label.cex <= 0){ 
			vertex.label.cex <- max(6*(nrow(sparseA))^(-0.8), 0.1) 
		}
		if (vertex.size <= 0){ 
			vertex.size <- max(75*(nrow(sparseA))^(-0.7), 1) 
		}
		
		# default plot title if not provided
		if (is.null(main)){ 
			main <- "VARX(1) time-series chain graph" 
		}

		# specify plot layout
		if (nrow(sparseA) %% 2 == 1){
			grid <- rbind(cbind(seq(0, -1, , length.out=nrow(sparseA)), 
						seq(0.5, floor(nrow(sparseA)/2) + 
						1/2, length.out=nrow(sparseA))), 
						cbind(seq(-1, 0, , length.out=ncol(sparseB)), 
						seq(ceiling(nrow(sparseA)/2) + 1/2, nrow(sparseA) + 
						1/2, length.out=ncol(sparseB))), 
						cbind(1, 1:nrow(sparseA)))
		}
		if (nrow(sparseA) %% 2 == 0){
			grid <- rbind(cbind(seq(0, -1, , length.out=nrow(sparseA)), 
						seq(0.5, floor(nrow(sparseA)/2), 
						length.out=nrow(sparseA))), 
						cbind(seq(-1, 0, , length.out=ncol(sparseB)), 
						seq(ceiling(nrow(sparseA)/2 + 1), nrow(sparseA) + 
						1/2, length.out=ncol(sparseB))), 
						cbind(1, 1:nrow(sparseA)))
		}


		# contempory independencies
		contCIG <- sparseP 
		contCIG[sparseP != 0] <- 1
		diag(contCIG) <- 0
		contCIG[lower.tri(contCIG)] <- 0

		if (side == "right"){
			# size of the edges       
			edge.widthA <- abs(sparseA)[which(sparseA != 0, arr.ind=TRUE)]
			edge.widthB <- abs(sparseB)[which(sparseB != 0, arr.ind=TRUE)] 
			edge.widthP <- abs(sparseP[which(contCIG != 0, arr.ind=TRUE)])
			edge.widthA  <- edge.widthA  / max(edge.widthA)
			edge.widthB  <- edge.widthB  / max(edge.widthB)
			edge.widthP  <- edge.widthP  / max(edge.widthP)
			if (edge.width <= 0){ 
				edge.width <- 40 * (nrow(sparseA)^(-0.8)) * c(edge.widthA, edge.widthB, edge.widthP)
			} else { 
				edge.width <- edge.width * c(edge.widthA, edge.widthB, edge.widthP)  
			}

			# create adjacency matrix			
			adjMat <- rbind(cbind(0*sparseA, 0*sparseB, t(sparseA)), 
					cbind(0*t(sparseB), matrix(0, ncol(sparseB), ncol(sparseB)), t(sparseB)), 
					cbind(0*sparseA, 0*sparseB, contCIG))
			adjMat[adjMat != 0] <- 1
        
			# convert adjacency matrix to graph-object
			gObj <- graph.adjacency(adjMat, mode="undirected")

			# arrow heads of curved edges should be zero width
			curved <- rep(0, sum(adjMat))
			curved[-c(1:sum(adjMat[1:(nNodesX + nNodesY),]))] <- -1      

			# specify negative edges  
			negEdgesA <- which(sparseA[which(sparseA != 0)] < 0)
			negEdgesB <- length(which(sparseA != 0)) + 
					which(sparseB[which(sparseB != 0)] < 0)
			negEdges <- c(negEdgesA, negEdgesB, length(which(sparseA != 0)) + 
					length(which(sparseB != 0)) + 
					which(sparseP[which(contCIG != 0, arr.ind=TRUE)] < 0))

			# specify positive edges
			posEdgesA <- which(sparseA[which(sparseA != 0)] > 0)
			posEdgesB <- length(which(sparseA != 0)) + 
					which(sparseB[which(sparseB != 0)] > 0)
			posEdges <- c(posEdgesA, posEdgesB, length(which(sparseA != 0)) + 
					length(which(sparseB != 0)) + 
					which(sparseP[which(contCIG != 0, arr.ind=TRUE)] > 0))

			# code sign of edges as solid/dashed
			igraph::E(gObj)[negEdges]$style <- 2
			igraph::E(gObj)[posEdges]$style <- 1
    
			# make plot
			plot(gObj, 
			     edge.lty=0, 
			     edge.arrow.size=0, 
			     layout=grid, 
			     vertex.shape="none", 
			     vertex.color="white", 
			     main=main, 
			     margin=c(0, -0.5, 0.1, 0.3), 
			     vertex.label.color="white", ...)
			plot(gObj, 
			     add=TRUE, 
			     layout=grid, 
			     vertex.size=rep(vertex.size, 2*nNodesY + nNodesX), 
			     vertex.label.font=vertex.label.font, 
			     edge.width=edge.width,
			     vertex.color=c(rep(vertex.color.T0, nNodesY), 
				            rep(vertex.color.X, nNodesX), 
					    rep(vertex.color.T1, nNodesY)),
			     vertex.frame.color=vertex.frame.color, 
			     vertex.label=c(nNamesY, nNamesX, nNamesY), 
			     vertex.label.cex=vertex.label.cex, 
			     vertex.label.color=c(rep(vertex.label.color.T0, nNodesY), 
						  rep(vertex.label.color.X, nNodesX), 
						  rep(vertex.label.color.T1, nNodesY)),
			     edge.color="black", 
			     edge.width=edge.width, 
			     vertex.label.family="sans", 
			     edge.curved=curved, 
			     edge.lty=igraph::E(gObj)$style, margin=c(0, -0.5, 0.1, 0.3), ...)
			text(0, -1.2, expression(paste(Y[italic(t)], "", sep="")), cex=1.2, font=3)
			text(0,  1.2, expression(paste(X[italic(t)], "", sep="")), cex=1.2)
			text(1, -1.2, expression(paste(Y[italic(t)+1], "", sep="")), cex=1.2)
		}

		if (side == "bottomleft"){
			# create adjacency matrix			
			adjMat <- rbind(cbind(contCIG, 0*sparseB, t(sparseA)), 
					cbind(0*t(sparseB), matrix(0, ncol(sparseB), ncol(sparseB)), t(sparseB)), 
					cbind(0*sparseA, 0*sparseB, 0*contCIG))
			adjMat[adjMat != 0] <- 1

			# convert adjacency matrix to graph-object
			gObj <- graph.adjacency(adjMat, mode="undirected")
            
			# width of the edges       
			Pslh <-sparseP
			Pslh[lower.tri(Pslh, diag=TRUE)] <- 0
			maxA <- max(abs(sparseA))
			maxB <- max(abs(sparseB))
			maxP <- max(abs(sparseP))
			adjMatSlh <- rbind(cbind(Pslh/maxP, 0*sparseB, t(sparseA)/maxA), 
						cbind(0*t(sparseB), matrix(0, ncol(sparseB), ncol(sparseB)), t(sparseB)/maxB), 
						cbind(0*sparseA, 0*sparseB, 0*contCIG))	

			# default node size and font size if not provided
			if (edge.width <= 0){ 
				edge.width <- 40 * (nrow(sparseA)^(-0.8)) * abs(t(adjMatSlh)[which(t(adjMatSlh) != 0)])
			} else { 
				edge.width <- edge.width * abs(t(adjMatSlh)[which(t(adjMatSlh) != 0)])
			}

			# arrow heads of curved edges should be zero width
			curved <- rep(0, sum(adjMat))
			slh <- which(adjMat[1, c(1:nrow(sparseP))]==1)
			if (length(slh) > 0){ idCurved <- 1:length(slh) } else { idCurved <- numeric() }
			for (j in 2:nrow(sparseP)){
				slh <- which(adjMat[j, c(1:nrow(sparseP))]==1)
				if (length(slh) > 0){
					idCurved <- c(idCurved, 1:length(slh) + sum(adjMat[1:(j-1),]))
				}
			}
			curved[idCurved] <- 1

			# identify negative and positive edges
			Pslh <-sparseP
			Pslh[lower.tri(Pslh, diag=TRUE)] <- 0
			adjMatSlh <- rbind(cbind(Pslh, 0*sparseB, t(sparseA)), 
						cbind(0*t(sparseB), matrix(0, ncol(sparseB), ncol(sparseB)), t(sparseB)), 
						cbind(0*sparseA, 0*sparseB, 0*contCIG))
			diag(adjMatSlh) <- 0
			edgeSigns <- sign(t(adjMatSlh)[which(t(adjMatSlh) != 0)])
			negEdges <- which(edgeSigns < 0)
			posEdges <- which(edgeSigns > 0)

			# code sign of edges as solid/dashed
			igraph::E(gObj)[negEdges]$style <- 2
			igraph::E(gObj)[posEdges]$style <- 1
    
			# make plot
			plot(gObj, 
			     edge.lty=0, 
			     edge.arrow.size=0, 
			     layout=grid, 
			     vertex.shape="none", 
			     vertex.color="white", 
			     main=main, 
			     margin=c(0, -0.5, 0.1, 0.3), 
			     vertex.label.color="white", ...)
			plot(gObj, 
			     add=TRUE, 
			     layout=grid, 
			     vertex.size=rep(vertex.size, 2*nNodesY + nNodesX), 
			     vertex.label.font=vertex.label.font, edge.width=edge.width,
			     vertex.color=c(rep(vertex.color.T0, nNodesY), 
					    rep(vertex.color.X, nNodesX), 
					    rep(vertex.color.T1, nNodesY)), 
			     vertex.frame.color=vertex.frame.color, 
			     vertex.label=c(nNamesY, nNamesX, nNamesY), 
			     vertex.label.cex=vertex.label.cex, 
			     vertex.label.color=c(rep(vertex.label.color.T0, nNodesY), 
			     			  rep(vertex.label.color.X, nNodesX), 
						  rep(vertex.label.color.T1, nNodesY)),
			     edge.color="black", 
			     edge.width=edge.width, 
			     vertex.label.family="sans", 
			     edge.curved=curved, 
			     edge.lty=igraph::E(gObj)$style, 
			     margin=c(0, -0.5, 0.1, 0.3), ...)
			text(0, -1.2, expression(paste(Y[italic(t)], "", sep="")), cex=1.2, font=3)
			text(0,  1.2, expression(paste(X[italic(t)], "", sep="")), cex=1.2)
			text(1, -1.2, expression(paste(Y[italic(t)+1], "", sep="")), cex=1.2)
		}
	}
    
	if (type=="ABonly"){
		# prune covariate without connections (according to sparseA, sparseB and sparseP)
		if (prune){ 
			idRemoveA <- intersect(which(apply(sparseA, 1, function(Z){ all(Z == 0) })), 
					       which(apply(sparseA, 2, function(Z){ all(Z == 0) })))
			idRemoveBy <- which(apply(sparseB, 1, function(Z){ all(Z == 0) }))
			idRemoveX <- which(apply(sparseB, 2, function(Z){ all(Z == 0) }))
			diag(sparseP) <- 0 
			idRemoveP <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), 
					       which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			idRemoveY <- intersect(intersect(idRemoveA, idRemoveBy), idRemoveP)
			if (length(idRemoveY) > 0){ 
				sparseA <- sparseA[-idRemoveY, -idRemoveY]
				sparseB <- sparseB[-idRemoveY, ]
				sparseP <- sparseP[-idRemoveY, -idRemoveY] 
				nNamesY <- nNamesY[-idRemoveY]
			}
			if (length(idRemoveX) > 0){ 
				sparseB <- sparseB[,-idRemoveX]
				nNamesX <- nNamesX[-idRemoveX]
			}
		}

		# store number of nodes of A and B
		nNodesY <- nrow(sparseA)
		nNodesX <- ncol(sparseB)

		# default node size and font size if not provided
		if (vertex.label.cex <= 0){ 
			vertex.label.cex <- max(6*(nrow(sparseA))^(-0.8), 0.1) 
		}
		if (vertex.size <= 0){ 
			vertex.size <- max(75*(nrow(sparseA))^(-0.7), 1) 
		}

		# default plot title if not provided
		if (is.null(main)){ 
			main <- "Cross-time and covariate relations of VARX(1) model" 
		}

		# specify plot layout
		if (nrow(sparseA) %% 2 == 1){
			grid <- rbind(cbind(seq(0, -1, , length.out=nrow(sparseA)), 
					seq(0.5, floor(nrow(sparseA)/2) + 1/2, length.out=nrow(sparseA))), 
					cbind(seq(-1, 0, , length.out=ncol(sparseB)), 
					seq(ceiling(nrow(sparseA)/2) + 1/2, nrow(sparseA) + 1/2, length.out=ncol(sparseB))), 
					cbind(1, 1:nrow(sparseA)))
		}
		if (nrow(sparseA) %% 2 == 0){
			grid <- rbind(cbind(seq(0, -1, , length.out=nrow(sparseA)), 
					seq(0.5, floor(nrow(sparseA)/2), length.out=nrow(sparseA))), 
					cbind(seq(-1, 0, , length.out=ncol(sparseB)), 
					seq(ceiling(nrow(sparseA)/2 + 1), nrow(sparseA) + 1/2, length.out=ncol(sparseB))), 
					cbind(1, 1:nrow(sparseA)))
		}
    
		# create adjacency matrix			
		adjMat <- rbind(cbind(0*contCIG, 0*sparseB, t(sparseA)), 
				cbind(0*t(sparseB), matrix(0, ncol(sparseB), ncol(sparseB)), t(sparseB)), 
				cbind(0*sparseA, 0*sparseB, 0*contCIG))
		adjMat[adjMat != 0] <- 1

		# convert adjacency matrix to graph-object
		gObj <- graph.adjacency(adjMat, mode="undirected")
            
		# width of the edges       
		Pslh <- 0 * sparseP
		Pslh[lower.tri(Pslh, diag=TRUE)] <- 0
		maxA <- max(abs(sparseA))
		maxB <- max(abs(sparseB))
		maxP <- max(abs(sparseP))
		adjMatSlh <- rbind(cbind(Pslh/maxP, 0*sparseB, t(sparseA)/maxA), 
				   cbind(0*t(sparseB), matrix(0, ncol(sparseB), ncol(sparseB)), t(sparseB)/maxB), 
				   cbind(0*sparseA, 0*sparseB, 0*contCIG))	

		# size of the edges       
		if (edge.width <= 0){ 
			edge.width <- 40 * (nrow(sparseA)^(-0.8)) * abs(t(adjMatSlh)[which(t(adjMatSlh) != 0)])
		} else { 
			edge.width <- edge.width * abs(t(adjMatSlh)[which(t(adjMatSlh) != 0)])
		}

		# identify negative and positive edges
		Pslh <- 0*sparseP
		Pslh[lower.tri(Pslh, diag=TRUE)] <- 0
		adjMatSlh <- rbind(cbind(0*Pslh, 0*sparseB, t(sparseA)), 
				   cbind(0*t(sparseB), matrix(0, ncol(sparseB), ncol(sparseB)), t(sparseB)), 
				   cbind(0*sparseA, 0*sparseB, 0*contCIG))
		diag(adjMatSlh) <- 0
		edgeSigns <- sign(t(adjMatSlh)[which(t(adjMatSlh) != 0)])
		negEdges <- which(edgeSigns < 0)
		posEdges <- which(edgeSigns > 0)

		# code sign of edges as solid/dashed
		igraph::E(gObj)[negEdges]$style <- 2
		igraph::E(gObj)[posEdges]$style <- 1
    
		# make plot
		plot(gObj, 
		     edge.lty=0, 
		     edge.arrow.size=0, 
		     layout=grid, 
		     vertex.shape="none", 
		     vertex.color="white", 
		     main=main, 
		     margin=c(0, -0.5, 0.1, 0.3), 
		     vertex.label.color="white", ...)
		plot(gObj, 
		     add=TRUE, 
		     layout=grid, 
		     vertex.size=rep(vertex.size, 2*nNodesY + nNodesX), 
		     vertex.label.font=vertex.label.font, edge.width=edge.width,
		     vertex.color=c(rep(vertex.color.T0, nNodesY), 
		                    rep(vertex.color.X, nNodesX), rep(vertex.color.T1, nNodesY)), 
		     vertex.frame.color=vertex.frame.color, 
		     vertex.label=c(nNamesY, nNamesX, nNamesY), 
		     vertex.label.cex=vertex.label.cex, 
		     vertex.label.color=c(rep(vertex.label.color.T0, nNodesY), 
			                  rep(vertex.label.color.X, nNodesX), 
					  rep(vertex.label.color.T1, nNodesY)),
		     edge.color="black", 
		     edge.width=edge.width, 
		     vertex.label.family="sans", 
		     edge.lty=igraph::E(gObj)$style, 
		     margin=c(0, -0.5, 0.1, 0.3), ...)
		text(0, -1.2, expression(paste(Y[italic(t)], "", sep="")), cex=1.2, font=3)
		text(0,  1.2, expression(paste(X[italic(t)], "", sep="")), cex=1.2)
		text(1, -1.2, expression(paste(Y[italic(t)+1], "", sep="")), cex=1.2)
	}
    
	if (type=="contempPC"){
		# prune covariate without connections (according to sparseA, sparseB and sparseP)
		if (prune){ 
			idRemoveA <- intersect(which(apply(sparseA, 1, function(Z){ all(Z == 0) })), 
					       which(apply(sparseA, 2, function(Z){ all(Z == 0) })))
			idRemoveBy <- which(apply(sparseB, 1, function(Z){ all(Z == 0) }))
			idRemoveX <- which(apply(sparseB, 2, function(Z){ all(Z == 0) }))
			diag(sparseP) <- 0 
			idRemoveP <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), 
					       which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			idRemoveY <- intersect(intersect(idRemoveA, idRemoveBy), idRemoveP)
			if (length(idRemoveY) > 0){ 
				sparseA <- sparseA[-idRemoveY, -idRemoveY]
				sparseB <- sparseB[-idRemoveY, ]
				sparseP <- sparseP[-idRemoveY, -idRemoveY] 
				nNamesY <- nNamesY[-idRemoveY]
			}
			if (length(idRemoveX) > 0){ 
				sparseB <- sparseB[,-idRemoveX]
				nNamesX <- nNamesX[-idRemoveX]
			}
		}

		# store number of nodes
		nNodes <- nrow(sparseP)

		# default node size and font size if not provided
		if (vertex.label.cex <= 0){ 
			vertex.label.cex <- max(8*(nrow(sparseP))^(-0.6), 0.1) 
		}
		if (vertex.size <= 0){ 
			vertex.size <- max(75*(nrow(sparseP))^(-0.6), 1) 
		}

		# default plot title if not provided
		if (is.null(main)){ 
			main <- "Contemp. cond. independence graph of VARX(1) model" 
		}
 
		# get global CIG
		CIG <- CIGofVAR1(sparseA, sparseP, type="contemp")
		diag(CIG) <- 0
		gObj <- graph.adjacency(CIG, "undirected")

		# size of the edges       
		if (edge.width <= 0){ edge.width <- 10 * (nrow(sparseP)^(-0.8)) }
		if (edge.width > 0){ edge.width <- edge.width }
        
		# actual plotting
		plot.igraph(gObj, 
			    layout=layout.circle, 
			    vertex.size=rep(vertex.size, nNodes), 
			    vertex.label.font=vertex.label.font, 
                	    vertex.color=rep(vertex.color.X, nNodes), 
			    vertex.frame.color=vertex.frame.color,
			    edge.color="black",
                	    vertex.label=nNamesY, 
			    vertex.label.cex=vertex.label.cex, 
			    vertex.label.color=rep(vertex.label.color.X, nNodes), 
                	    edge.width=edge.width, 
			    vertex.label.family="sans", 
			    main=main, ...)               
	}
}

