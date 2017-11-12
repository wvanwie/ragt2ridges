graphVAR2 <- function(sparseA1, 
			sparseA2, 
			sparseP, 
			type="TSCG", 
			side="right", 
			prune=TRUE, 
			nNames=NULL, 
			main=NULL, 
			vertex.color.T0="lightcyan2", 
			vertex.color.T1="lightcyan2", 
			vertex.color.T2="lightcyan2", 
			vertex.frame.color="steelblue", 
			vertex.label.cex=-1, 
			vertex.label.color.T0="black", 
			vertex.label.color.T1="black", 
			vertex.label.color.T2="black", 
			vertex.label.font=1.5, 
			vertex.size=-1, 
			edge.arrow.size=-1, 
			edge.width=-1, 
                        ...){

	#######################################################################
	#
	# DESCRIPTION:
	# Plot temporal relations as implied by A1 and A2 (the matrice with 
	# lag one and two regression coefficient of the VAR(2)) model.
	#
	# ARGUMENTS:
	#
	# -> sparseA1                : Matrix A1 of lag one regression 
	#                              parameters, which is assumed to be sparse.
	# -> sparseA2                : Matrix A2 of lag two regression 
	#                              parameters, which is assumed to be sparse.
	# -> sparseP                 : Matrix P of precision of the error, 
	#                              which is assumed to be sparse.
	# -> type                    : A 'character' indicating what should 
	#                              be plotted. If type='TSCG'. the time
	#                              series chain graph is plotted, while 
	#                              type type='Aonly' limits this graph
	#                              to the temporal relations. If 
	#                              type='globalPC' or type='contempPC', the 
	#                              global or contemporaneous (respectively) 
	#                              partial correlation graph is plotted.
	# -> side                    : A 'character' indicating wether the 
	#                              contemporaneous dependencies should be 
	#                              plotted on the 'left' (time t) or the 
	#                              'right' (time t+1) side. Only active 
	#                              when type='TSCG'.
	# -> prune                   : A 'logical' indicating whether to 
	#                              remove covariates without any 
	#                              temporal relations (as implied 
	#                              by 'sparseA1', 'sparseA2', 'sparseP').
	# -> nNames                  : A 'character' containing the 
	#                              covariate names to written inside 
	#                              the nodes.
	# -> main                    : The 'character' to be plotted as 
	#                              title above the graph.
	# -> vertex.color.T0         : Color of nodes at time point t.
	# -> vertex.color.T1         : Color of nodes at time point t+1. 
	#                              Ignored when type='globalPC' 
	#                              or type='contempPC'.
	# -> vertex.color.T2         : Color of nodes at time point t+2. 
	#                              Ignored when type='globalPC' 
	#                              or type='contempPC'.
	# -> vertex.frame.color      : Refer to 'plot.igraph'.
	# -> vertex.label.cex        : Refer to 'plot.igraph'.
	# -> vertex.label.color.T0   : Color of the node label at time 
	#                              point t.
	# -> vertex.label.color.T1   : Color of the node label at time 
	#                              point t+1. Ignored when 
	#                              type='globalPC' or type='contempPC'.
	# -> vertex.label.color.T2   : Color of the node label at time 
	#                              point t+2. Ignored when 
	#                              type='globalPC' or type='contempPC'.
	# -> vertex.label.font       : Refer to 'plot.igraph'.
	# -> vertex.size             : Refer to 'plot.igraph'.
	# -> edge.arrow.size         : Refer to 'plot.igraph'.
	# -> edge.width              : Refer to 'plot.igraph'.
	# -> ...                     : Other arguments to be passed to 
	#                              'plot.igraph'.
	# 
	# DEPENDENCIES:
	# library(igraph)	     : functions: delete.edges, E, ecount, 
	#                              plot.igraph, graph.adjacency, 
	#                              maximum.cardinality.search
	#
	# NOTES:
   	# - igraph does not support visualization of mixed graphs, including 
	#   both directed and undirected edges. As a consequence the 
	#   time-series chain graph is undirected: the edges from t to t+1 or
	#   those from t to t+2 should be thought of as directed. 
	# - if "white" labels, color combination "navy" and "blue" looks nice. 
	# - if "black" labels, color combination "tomato1" and "red" looks nice. 
	#     	
	#######################################################################
	# input check
	if (as.character(class(sparseA1)) != "matrix"){ 
		stop("Input (sparseA1) is of wrong class.") 
	}
	if (nrow(sparseA1) != ncol(sparseA1)){ 
		stop("Matrix sparseA1 is not square.") 
	}
	if (as.character(class(sparseA2)) != "matrix"){ 
		stop("Input (sparseA2) is of wrong class.") 
	}
	if (nrow(sparseA2) != ncol(sparseA2)){ 
		stop("Matrix sparseA2 is not square.") 
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
	if (as.character(class(vertex.color.T0)) != "character"){ 
		stop("Input (vertex.color.T0) is of wrong class.") 
	}
	if (as.character(class(vertex.color.T1)) != "character"){ 
		stop("Input (vertex.color.T1) is of wrong class.") 
	}
	if (as.character(class(vertex.label.color.T0)) != "character"){ 
		stop("Input (vertex.label.color.T0) is of wrong class.") 
	}
	if (as.character(class(vertex.label.color.T1)) != "character"){ 
		stop("Input (vertex.label.color.T1) is of wrong class.") 
	}

	# if no covariate names are specified the columns and 
	# row names of A1 and A2 and P are given names 1, 2, et cetera
	if (is.null(nNames)){ nNames <- as.character(1:nrow(sparseA1)) }

	if (type=="TSCG"){     
		# prune covariate without connections (according to sparseA1, sparseA2 and sparseP)
		if (prune){ 
			idRemoveA1 <- intersect(which(apply(sparseA1, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA1, 2, function(Z){ all(Z == 0) })))
			idRemoveA2 <- intersect(which(apply(sparseA2, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA2, 2, function(Z){ all(Z == 0) })))
			diag(sparseP) <- 0 
			idRemoveP <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			idRemove <- intersect(intersect(idRemoveA1, idRemoveA2), idRemoveP)
			if (length(idRemove) > 0){ 
				sparseA1 <- sparseA1[-idRemove, -idRemove]
				sparseA2 <- sparseA2[-idRemove, -idRemove]
				sparseP <- sparseP[-idRemove, -idRemove] 
				nNames <- nNames[-idRemove]
			}
		}
	
		# store number of nodes
		nNodes <- nrow(sparseA1)
	
		# default node size and font size if not provided
		if (vertex.label.cex <= 0){ 
			vertex.label.cex <- max(6*(nrow(sparseA1))^(-0.8), 0.1) 
		}
		if (vertex.size <= 0){ 
			vertex.size <- max(75*(nrow(sparseA1))^(-0.7), 1) 
		}

		# default plot title if not provided
		if (is.null(main)){ 
			main <- "VAR(2) time-series chain graph" 
		}

		# specify plot layout
		if (nrow(sparseA1) %% 2 == 1){
			grid <- rbind(cbind(seq(0, -1, , length.out=nrow(sparseA1)), 
					seq(0.5, floor(nrow(sparseA1)/2) + 1/2, length.out=nrow(sparseA1))), 
					cbind(seq(-1, 0, , length.out=nrow(sparseA1)), 
					seq(ceiling(nrow(sparseA1)/2) + 1/2, nrow(sparseA1) + 1/2, 
					length.out=nrow(sparseA1))), cbind(1, 1:nrow(sparseA1)))
		}
		if (nrow(sparseA1) %% 2 == 0){
			grid <- rbind(cbind(seq(0, -1, , length.out=nrow(sparseA1)), 
					seq(0.5, floor(nrow(sparseA1)/2), length.out=nrow(sparseA1))), 
					cbind(seq(-1, 0, , length.out=nrow(sparseA1)), 
					seq(ceiling(nrow(sparseA1)/2 + 1), nrow(sparseA1) + 1/2, 
					length.out=nrow(sparseA1))), cbind(1, 1:nrow(sparseA1)))
		}

		# contempory independencies
		contCIG <- sparseP 
		contCIG[sparseP != 0] <- 1
		diag(contCIG) <- 0
		contCIG[lower.tri(contCIG)] <- 0

		if (side == "right"){
			# adjacency matrix 
			adjMat <- rbind(cbind(0*sparseA1, 0*sparseA2, t(sparseA1)), 	
					cbind(0*sparseA1, 0*sparseA2, t(sparseA2)), 
					cbind(0*sparseA1, 0*sparseA2, contCIG))
			adjMat[adjMat != 0] <- 1
			
			# names 
			rownames(adjMat) <- colnames(adjMat) <- 
						c(paste(nNames, "l1", sep="_"), 
							paste(nNames, "l2", sep="_"), 
							paste(nNames, "p", sep="_"))
         
			# size of the edges       
			edge.widthA1 <- abs(sparseA1)[which(sparseA1 != 0, arr.ind=TRUE)]
			edge.widthA2 <- abs(sparseA2)[which(sparseA2 != 0, arr.ind=TRUE)] 
			edge.widthP <- abs(sparseP[which(contCIG != 0, arr.ind=TRUE)])
			edge.widthA1  <- edge.widthA1  / max(c(edge.widthA1, edge.widthA2))
			edge.widthA2  <- edge.widthA2  / max(c(edge.widthA2, edge.widthA1))
			edge.widthP  <- edge.widthP  / max(edge.widthP)
			if (edge.width <= 0){ 
				edge.width <- 40 * (nrow(sparseA1)^(-0.8)) * c(edge.widthA1, edge.widthA2, edge.widthP)
			} else { 
				edge.width <- edge.width * c(edge.widthA1, edge.widthA2, edge.widthP)  
			}
			
			# convert adjacency matrix to graph-object
			gObj <- graph.adjacency(adjMat, mode="undirected")

			# arrow heads of curved edges should be zero width
			curved <- rep(0, sum(adjMat))
			curved[-c(1:sum(adjMat[1:(2*nNodes),]))] <- -1  

			# sort out negative edges
			negEdgesA1 <- which(sparseA1[which(sparseA1 != 0)] < 0)
			negEdgesA2 <- length(which(sparseA1 != 0)) + which(sparseA2[which(sparseA2 != 0)] < 0)
			negEdges <- c(negEdgesA1, negEdgesA2, length(which(sparseA1 != 0)) + 
								length(which(sparseA2 != 0)) + 
								which(sparseP[which(contCIG != 0, arr.ind=TRUE)] < 0))

			# sort out positive edges
			posEdgesA1 <- which(sparseA1[which(sparseA1 != 0)] > 0)
			posEdgesA2 <- length(which(sparseA1 != 0)) + which(sparseA2[which(sparseA2 != 0)] > 0)
			posEdges <- c(posEdgesA1, posEdgesA2, length(which(sparseA1 != 0)) + 
								length(which(sparseA2 != 0)) + 
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
				vertex.size=rep(vertex.size, 3*nNodes), 
				vertex.label.font=vertex.label.font, 
				edge.width=edge.width,
				vertex.color=c(rep(vertex.color.T1, nNodes), 
						rep(vertex.color.T0, nNodes), 
						rep(vertex.color.T2, nNodes)), 
				vertex.frame.color=vertex.frame.color, 
				vertex.label=rep(nNames, 3), 
				vertex.label.cex=vertex.label.cex, 
				vertex.label.color=c(rep(vertex.label.color.T1, nNodes), 
							rep(vertex.label.color.T0, nNodes), 
							rep(vertex.label.color.T2, nNodes)),  
				edge.color="black", 
				edge.width=edge.width, 
				vertex.label.family="sans", 
				edge.curved=curved, 
				edge.lty=igraph::E(gObj)$style, 
				margin=c(0, -0.5, 0.1, 0.3), ...)
			text(0, -1.2, "t", cex=1.2, font=3)
			text(0,  1.2, expression(paste(italic(t), "-1 ", sep="")), cex=1.2)
			text(1, -1.2, expression(paste(italic(t), "+1 ", sep="")), cex=1.2)
		}

		if (side == "bottomleft"){
			# adjacency matrix 
			adjMat <- rbind(cbind(contCIG, 0*sparseA2, t(sparseA1)), 
					cbind(0*sparseA1, 0*contCIG, t(sparseA2)), 
					cbind(0*sparseA1, 0*sparseA2, 0*contCIG))
			adjMat[adjMat != 0] <- 1

			# convert adjacency matrix to graph-object
			gObj <- graph.adjacency(adjMat, mode="undirected")
            
			# size of the edges       
			Pslh <-sparseP
			Pslh[lower.tri(Pslh, diag=TRUE)] <- 0
			maxAs <- max(max(abs(t(sparseA1))), max(abs(t(sparseA2))))
			adjMatSlh <- rbind(cbind(abs(Pslh), 0*sparseA2, abs(t(sparseA1))/maxAs), 
						cbind(0*sparseA1, 0*contCIG, abs(t(sparseA2))/maxAs), 
						cbind(0*sparseA1, 0*sparseA2, 0*contCIG))
			if (edge.width <= 0){ 
				edge.width <- 40 * (nrow(sparseA1)^(-0.8)) * abs(t(adjMatSlh)[which(t(adjMatSlh) != 0)])
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
			adjMatSlh <- rbind(cbind(Pslh, 0*sparseA2, t(sparseA1)), 
						cbind(0*sparseA1, 0*contCIG, t(sparseA2)), 
						cbind(0*sparseA1, 0*sparseA2, 0*contCIG))
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
				vertex.size=rep(vertex.size, 3*nNodes), 
				vertex.label.font=vertex.label.font, 
				edge.width=edge.width,
				vertex.color=c(rep(vertex.color.T1, nNodes), 
						rep(vertex.color.T0, nNodes), 
						rep(vertex.color.T2, nNodes)), 
				vertex.frame.color=vertex.frame.color, 
				vertex.label=c(nNames, nNames, nNames), 
				vertex.label.cex=vertex.label.cex, 
				vertex.label.color=c(rep(vertex.label.color.T1, nNodes), 
							rep(vertex.label.color.T0, nNodes), 
							rep(vertex.label.color.T2, nNodes)), 
				edge.color="black", 
				edge.width=edge.width, 
				vertex.label.family="sans", 
				edge.curved=curved, 
				edge.lty=igraph::E(gObj)$style, 
				margin=c(0, -0.5, 0.1, 0.3), ...)
			text(0, -1.2, expression(paste(Y[italic(t)+1], "", sep="")), cex=1.2, font=3)
			text(0,  1.2, expression(paste(Y[italic(t)], "", sep="")), cex=1.2)
			text(1, -1.2, expression(paste(Y[italic(t)+2], "", sep="")), cex=1.2)
		}


		if (side == "topleft"){
			# adjacency matrix 
			adjMat <- rbind(cbind(0*contCIG, 0*sparseA2, t(sparseA1)), 
					cbind(0*sparseA1, contCIG, t(sparseA2)), 
					cbind(0*sparseA1, 0*sparseA2, 0*contCIG))
			adjMat[adjMat != 0] <- 1

			# convert adjacency matrix to graph-object
			gObj <- graph.adjacency(adjMat, mode="undirected")
            
			# size of the edges       
			Pslh <-sparseP
			Pslh[lower.tri(Pslh, diag=TRUE)] <- 0
			maxAs <- max(max(abs(t(sparseA1))), max(abs(t(sparseA2))))
			adjMatSlh <- rbind(cbind(0*abs(Pslh), 0*sparseA2, abs(t(sparseA1))/maxAs), 
						cbind(0*sparseA1, abs(Pslh), abs(t(sparseA2))/maxAs), 
						cbind(0*sparseA1, 0*sparseA2, 0*contCIG))
			if (edge.width <= 0){ 
				edge.width <- 40 * (nrow(sparseA1)^(-0.8)) * t(adjMatSlh)[which(t(adjMatSlh) != 0)]
			} else { 
				edge.width <- edge.width * t(adjMatSlh)[which(t(adjMatSlh) != 0)]
			}

			# arrow heads of curved edges should be zero width
			curved <- rep(0, sum(adjMat))  
			Pslh <-sparseP
			Pslh[lower.tri(Pslh, diag=TRUE)] <- 0
			adjMatSlh <- rbind(cbind(0*abs(Pslh), 0*sparseA2, abs(t(sparseA1))), 
						cbind(0*sparseA1, -abs(Pslh), abs(t(sparseA2))), 
						cbind(0*sparseA1, 0*sparseA2, 0*contCIG))
			curved[which(t(adjMatSlh)[which(t(adjMatSlh) != 0)] < 0)] <- 1

			# identify negative and positive edges
			Pslh <-sparseP
			Pslh[lower.tri(Pslh, diag=TRUE)] <- 0
			adjMatSlh <- rbind(cbind(0*abs(Pslh), 0*sparseA2, t(sparseA1)), 
						cbind(0*sparseA1, abs(Pslh), t(sparseA2)), 
						cbind(0*sparseA1, 0*sparseA2, 0*contCIG))
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
				vertex.size=rep(vertex.size, 3*nNodes), 
				vertex.label.font=vertex.label.font, 
				edge.width=edge.width,
				vertex.color=c(rep(vertex.color.T1, nNodes), 
						rep(vertex.color.T0, nNodes), 
						rep(vertex.color.T2, nNodes)), 
				vertex.frame.color=vertex.frame.color, 
				vertex.label=c(nNames, nNames, nNames), 
				vertex.label.cex=vertex.label.cex, 
				vertex.label.color=c(rep(vertex.label.color.T1, nNodes), 
							rep(vertex.label.color.T0, nNodes), 
							rep(vertex.label.color.T2, nNodes)),
				edge.color="black", 
				edge.width=edge.width, 
				vertex.label.family="sans", 
				edge.curved=curved, 
				edge.lty=igraph::E(gObj)$style, 
				margin=c(0, -0.5, 0.1, 0.3), ...)
			text(0, -1.2, expression(paste(Y[italic(t)+1], "", sep="")), cex=1.2, font=3)
			text(0,  1.2, expression(paste(Y[italic(t)], "", sep="")), cex=1.2)
			text(1, -1.2, expression(paste(Y[italic(t)+2], "", sep="")), cex=1.2)
		}



	}
    
	if (type=="A12only"){
		# prune covariate without connections (according to sparseA1, sparseA2 and sparseP)
		if (prune){ 
			idRemoveA1 <- intersect(which(apply(sparseA1, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA1, 2, function(Z){ all(Z == 0) })))
			idRemoveA2 <- intersect(which(apply(sparseA2, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA2, 2, function(Z){ all(Z == 0) })))
			diag(sparseP) <- 0 
			idRemoveP <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			idRemove <- intersect(intersect(idRemoveA1, idRemoveA2), idRemoveP)
			if (length(idRemove) > 0){ 
				sparseA1 <- sparseA1[-idRemove, -idRemove]
				sparseA2 <- sparseA2[-idRemove, -idRemove]
				sparseP <- sparseP[-idRemove, -idRemove] 
				nNames <- nNames[-idRemove]
			}
		}

		# store number of nodes		
		nNodes <- nrow(sparseA1)

		# default node size and font size if not provided
		if (vertex.label.cex <=0 ){ 
			vertex.label.cex <- max(6*(nrow(sparseA1))^(-0.8), 0.1) 
		}
		if (vertex.size <= 0){ 
			vertex.size <- max(75*(nrow(sparseA1))^(-0.7), 1) 
		}

		# default plot title if not provided
		if (is.null(main)){ 
			main <- "Cross-time relations of VAR(2) model" 
		}

		# specify plot layout
		if (nrow(sparseA1) %% 2 == 1){
			grid <- rbind(cbind(seq(0, -1, , length.out=nrow(sparseA1)), 
					seq(0.5, floor(nrow(sparseA1)/2) + 1/2, length.out=nrow(sparseA1))), 
					cbind(seq(-1, 0, , length.out=nrow(sparseA1)), 
					seq(ceiling(nrow(sparseA1)/2) + 1/2, nrow(sparseA1) + 1/2, 
					length.out=nrow(sparseA1))), cbind(1, 1:nrow(sparseA1)))
		}
		if (nrow(sparseA1) %% 2 == 0){
			grid <- rbind(cbind(seq(0, -1, , length.out=nrow(sparseA1)), 
					seq(0.5, floor(nrow(sparseA1)/2), length.out=nrow(sparseA1))), 
					cbind(seq(-1, 0, , length.out=nrow(sparseA1)), 
					seq(ceiling(nrow(sparseA1)/2 + 1), nrow(sparseA1) + 1/2, 
					length.out=nrow(sparseA1))), cbind(1, 1:nrow(sparseA1)))
		}
    
		# adjacency matrix 
		adjMat <- rbind(cbind(0*sparseA1, 0*sparseA2, t(sparseA1)), 
				cbind(0*sparseA1, 0*sparseA2, t(sparseA2)), 
				cbind(0*sparseA1, 0*sparseA2, 0*sparseP))
		adjMat[adjMat != 0] <- 1

		# convert adjacency matrix to graph-object
		gObj <- graph.adjacency(adjMat)

		# arrow heads of curved edges should be zero width
		negEdgesA1 <- which(sparseA1[which(sparseA1 != 0)] < 0)
		negEdgesA2 <- length(which(sparseA1 != 0)) + which(sparseA2[which(sparseA2 != 0)] < 0)
		negEdges <- c(negEdgesA1, negEdgesA2)
		posEdgesA1 <- which(sparseA1[which(sparseA1 != 0)] > 0)
		posEdgesA2 <- length(which(sparseA1 != 0)) + which(sparseA2[which(sparseA2 != 0)] > 0)
		posEdges <- c(posEdgesA1, posEdgesA2)

		# code sign of edges as solid/dashed
		igraph::E(gObj)[negEdges]$style <- 2
		igraph::E(gObj)[posEdges]$style <- 1

		# size of the edges       
		edge.widthA1 <- abs(sparseA1)[which(sparseA1 != 0, arr.ind=TRUE)]
		edge.widthA2 <- abs(sparseA2)[which(sparseA2 != 0, arr.ind=TRUE)] 
		edge.widthA1  <- edge.widthA1  / max(c(edge.widthA1, edge.widthA2))
		edge.widthA2  <- edge.widthA2  / max(c(edge.widthA2, edge.widthA1))
		if (edge.width <= 0){ 
			edge.width <- 40 * (nrow(sparseA1)^(-0.8)) * c(edge.widthA1, edge.widthA2)
		} else { 
			edge.width <- edge.width * c(edge.widthA1, edge.widthA2)  
		}
		if (edge.arrow.size <= 0){ 
			edge.arrow.size <- 10 * (nrow(sparseA1)^(-0.8)) * c(edge.widthA1, edge.widthA2) 
		} else { 
			edge.arrow.size <- edge.arrow.size * c(edge.widthA1, edge.widthA2)   
		}

		# make plot
		plot(gObj, 
			edge.lty=0, 
			edge.arrow.size=0, 
			layout=grid, 
			vertex.shape="none", 
			vertex.label=c(nNames, nNames), 
			vertex.label.family="sans", 
			main=main, 
			vertex.label.cex=vertex.label.cex, ...)
		for (e in seq_len(ecount(gObj))){
			graph2 <- delete.edges(gObj, igraph::E(gObj)[(1:ecount(gObj))[-e]])
			plot(graph2, 
				edge.arrow.size=edge.arrow.size[e], 
				layout=grid, 
				vertex.size=rep(vertex.size, 3*nNodes), 
				vertex.label.font=vertex.label.font, 
				vertex.color=c(rep(vertex.color.T1, nNodes), 
						rep(vertex.color.T0, nNodes), 
						rep(vertex.color.T2, nNodes)), 
				vertex.frame.color=vertex.frame.color,
				layout=grid, 
				edge.color="black", 
				vertex.label.cex=vertex.label.cex, 
				edge.width=edge.width[e], 
				vertex.label.family="sans", 
				edge.lty=igraph::E(graph2)$style,     
				vertex.label=c(nNames, nNames, nNames), 
				add=TRUE, ...)        
		}       
		text(0, -1.2, expression(paste(Y[italic(t)+1], "", sep="")), cex=1.2, font=3)
		text(0,  1.2, expression(paste(Y[italic(t)], "", sep="")), cex=1.2)
		text(1, -1.2, expression(paste(Y[italic(t)+2], "", sep="")), cex=1.2)
	}
    
	if (type=="globalPC"){     
		# prune covariate without connections (according to sparseA and sparseP)
		if (prune){ 
			idRemoveA1 <- intersect(which(apply(sparseA1, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA1, 2, function(Z){ all(Z == 0) })))
			idRemoveA2 <- intersect(which(apply(sparseA2, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA2, 2, function(Z){ all(Z == 0) })))
			diag(sparseP) <- 0 
			idRemoveP <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			idRemove <- intersect(intersect(idRemoveA1, idRemoveA2), idRemoveP)
			if (length(idRemove) > 0){ 
				sparseA1 <- sparseA1[-idRemove, -idRemove]
				sparseA2 <- sparseA2[-idRemove, -idRemove]
				sparseP <- sparseP[-idRemove, -idRemove] 
				nNames <- nNames[-idRemove]
			}
		}

		# store number of nodes
		nNodes <- nrow(sparseA1)

		# default node size and font size if not provided
		if (vertex.label.cex <= 0){ 
			vertex.label.cex <- max(8*(nrow(sparseA1))^(-0.6), 0.1) 
		}
		if (vertex.size <= 0){ 
			vertex.size <- max(75*(nrow(sparseA1))^(-0.6), 1) 
		}

		# default plot title if not provided
		if (is.null(main)){ 
			main <- "Partial correlation graph of VAR(2) model" 
		}

		# create graph object      
		CIG <- CIGofVAR2(sparseA1, sparseA2, sparseP, type="global")  
		diag(CIG) <- 0         
		gObj <- graph.adjacency(CIG, "undirected")

		# size of the edges       
		if (edge.width <= 0){ 
			edge.width <- 10 * (nrow(sparseA1)^(-0.8)) 
		} else {
			edge.width <- edge.width 
		}
        
		# actual plotting
		plot.igraph(gObj, 
			layout=layout.circle, 
			vertex.size=rep(vertex.size, nNodes), 
			vertex.label.font=vertex.label.font, 
			vertex.color=rep(vertex.color.T0, nNodes), 
			vertex.frame.color=vertex.frame.color, 
			edge.color="black",
			vertex.label=nNames, 
			vertex.label.cex=vertex.label.cex, 
			vertex.label.color=rep(vertex.label.color.T0, nNodes), 
			edge.width=edge.width, 
			vertex.label.family="sans", 
			main=main, ...)       
	}
    
	if (type=="contempPC"){
		# prune covariate without connections (according to sparseA1, sparseA2 and sparseP)
		if (prune){ 
			idRemoveA1 <- intersect(which(apply(sparseA1, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA1, 2, function(Z){ all(Z == 0) })))
			idRemoveA2 <- intersect(which(apply(sparseA2, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseA2, 2, function(Z){ all(Z == 0) })))
			diag(sparseP) <- 0 
			idRemoveP <- intersect(which(apply(sparseP, 1, function(Z){ all(Z == 0) })), 
						which(apply(sparseP, 2, function(Z){ all(Z == 0) })))
			idRemove <- intersect(intersect(idRemoveA1, idRemoveA2), idRemoveP)
			if (length(idRemove) > 0){ 
				sparseA1 <- sparseA1[-idRemove, -idRemove]
				sparseA2 <- sparseA2[-idRemove, -idRemove]
				sparseP <- sparseP[-idRemove, -idRemove] 
				nNames <- nNames[-idRemove]
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
			main <- "Contemp. cond. independence graph of VAR(2) model" 
		}
 
		# get global CIG
		CIG <- CIGofVAR2(sparseA1, sparseA2, sparseP, type="contemp")
		diag(CIG) <- 0
		gObj <- graph.adjacency(CIG, "undirected")

		# size of the edges       
		if (edge.width <= 0){ 
			edge.width <- 10 * (nrow(sparseP)^(-0.8)) 
		}
		if (edge.width > 0){ 
			edge.width <- edge.width 
		}
        
		# actual plotting
		plot.igraph(gObj, 
			layout=layout.circle, 
			vertex.size=rep(vertex.size, nNodes), 
			vertex.label.font=vertex.label.font, 
                	vertex.color=rep(vertex.color.T0, nNodes), 
			vertex.frame.color=vertex.frame.color, 
			edge.color="black",
                	vertex.label=nNames, 
			vertex.label.cex=vertex.label.cex, 
			vertex.label.color=rep(vertex.label.color.T0, nNodes), 
                	edge.width=edge.width, 
			vertex.label.family="sans", 
			main=main, ...)               
	}
}

