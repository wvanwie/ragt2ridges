////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
SEASIDE VERSION OF THE DEFAULT TARGET FUNCTION (AS PROVIDED THROUGH THE RAGS2RIDGES-PACKAGE).
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline arma::mat armaP_defaultTarget(arma::mat S, std::string targetType, const double fraction, double const multiplier){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Rcpp version of its R-counterpart 'default.target'.
	// Key difference is the input checks.
	/////////////////////////////////////////////////////////////////////////////////////////

	// Diagonal matrix with average of inverse nonzero eigenvalues of S as entries
	if (targetType == "DAIE"){
		arma::vec eigvals;
		arma::mat eigvecs;
		arma::eig_sym(eigvals, eigvecs, S, "dc");
		double slh = mean(1/eigvals(find(eigvals >= (max(eigvals) * fraction))));
		S = arma::zeros(S.n_rows, S.n_rows);
		S.diag() += slh;
	}

	// Diagonal matrix with inverse of average of eigenvalues of S as entries
	if (targetType == "DIAES"){
		arma::vec eigvals;
		arma::mat eigvecs;
		arma::eig_sym(eigvals, eigvecs, S, "dc");
		const double slh = 1/mean(eigvals);
		S = arma::zeros(S.n_rows, S.n_rows);
		S.diag() += slh;		
	}

	// Diagonal matrix with unit partial variance as entries
	if (targetType == "DUPV"){
		S = arma::eye(S.n_rows, S.n_rows);
	}

	// Diagonal matrix with average empirical partial variances as entries
	if (targetType == "DAPV"){
		const double apv = mean(1/arma::diagvec(S));
		S = arma::zeros(S.n_rows, S.n_rows);
		S.diag() += apv;		
	}

	// Diagonal matrix with constant partial variance as entries
    	if (targetType == "DCPV"){
    		S = arma::zeros(S.n_rows, S.n_rows);
    		S.diag() += multiplier;
    }

    // Diagonal matrix with empirical partial variances as entries
    if (targetType == "DEPV"){
        S = arma::diagmat(1/arma::diagvec(S));
	}

    // Null matrix
    if (targetType == "Null"){
        S = arma::zeros(S.n_rows, S.n_rows);
    }

    // Return
    return(S);
}

