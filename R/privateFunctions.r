.is.int <- function (x, tolerance = .Machine$double.eps){
    abs(x - round(x)) < tolerance
}


.LL <- function (S, P){
    LL <- -log(det(P)) + .trace(S %*% P)
    return(LL)
}

.trace <- function (M) {
    return(sum(diag(M)))
}

.isStationary <- function(A){
	#############################################################################
	# tests whether the VAR(1) process associated with matrix A is stationary
	#############################################################################
	
	# assess stationarity	
    	evs <- abs(eigen(A, only.values=TRUE)$values)
    	if (max(evs) > 1){ print("WARNING: the VAR(1) process associated with the supplied matrix A is not stable!") }
}

