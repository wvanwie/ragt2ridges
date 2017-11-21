pruneMotifStats <- function(motifList, id){
	#######################################################################
	#
	# DESCRIPTION:
	# -> Function that subsets the motif list as returned by the
	#    motifStatsVAR1-function
	#
	# ARGUMENTS:
	# -> motiflist : A motif list as returned by the motifStatsVAR1-function.
	# -> id        : An integer representing the node of interest.
	# 
	# DEPENDENCIES:
	# ...          : ...
	#
	# REFERENCES:
	# -> Alon, U. (2007), "Network motifs: theory and experimental approaches", 
	#                     Nature Reviews Genetics, 8, 450-461.
	# 
	#######################################################################

	# input checks
	if (!is.list(motifList)){
		stop("Input (motifList) should be a list.")
	}

	# subsetting
	for (k in 1:length(motifList)){
		idInvolved <- which(unlist(lapply(motifList[[k]], 
                                                  function(Z, id){ any(as.numeric(Z[,-ncol(Z)]) == id) }, 
                                                  id=id)))
		motifList[[k]] <- motifList[[k]][idInvolved]
	}
	return(motifList)
}

