longitudinal2array <- function(Y){
	########################################################################
	#
	# DESCRIPTION:
	# Converts a 3-dim array (containing time-series data of multiple 
	# individuals) to an object of the 'longitudinal' class.
	#
	# ARGUMENTS:
	# -> Y			: Object of class 'longitudinal'. Essentially, 
	#                         this is a matrix with the 'timepoint'-slices 
	#                         of the array stacked on top of each other.
	#
	#
	# NOTES:
	# ...
	#
	########################################################################

	# reconstruct design matrix from rownames of longitudinal object
	design <- matrix(unlist(strsplit(rownames(Y), "-")), ncol=2, byrow=TRUE)
	design <- cbind(as.numeric(design[,1]), as.numeric(design[,2]))
	colnames(design) <- c("time", "sample")

	# recover time and sample dimension of experiment
	times <- unique(design[,1])
	T <- length(times)
	n <- length(unique(design[,2]))

	# reformat data as array
	Yarray <- array(NA, dim = c(ncol(Y), T, n))
	for (i in 1:n){
		timeSlh <- design[which(design[,2] == i),1]
		Ytemp <- t(Y[which(design[,2] == i),])
		Ytemp <- Ytemp[, order(timeSlh)]
		Yarray[ , match(timeSlh, times), i] <- Ytemp
	}
	
	# return array
	return(Yarray)
}
