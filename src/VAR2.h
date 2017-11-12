////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
FUCTIONS FOR THE ESTIMATION OF THE VAR(2) MODEL: FOR USE ON THE SEASIDE, FOLLOWED BY WRAPPERS FOR DEPORTATION TO R.
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double armaVAR2_convergenceEvaluation(arma::mat& Ahat, arma::mat& Aprev, arma::mat& Phat, arma::mat& Pprev){
	/////////////////////////////////////////////////////////////////////////////////////////
	// assess maximum difference between the current and previous estimates
	/////////////////////////////////////////////////////////////////////////////////////////
	
	// extract element-wise maximum of both matrices
	arma::mat maxis = arma::max(arma::max(abs(Ahat - Aprev), 1), arma::max(abs(Phat - Pprev), 1));
	return maxis.max();

}

inline arma::mat armaVAR2_COVYhat(arma::cube& Y){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Covariance estimation of the multivariate time series data.
	/////////////////////////////////////////////////////////////////////////////////////////

	// variable declaration and initiation 
	int COVZdofLag1 = 0; int COVZdofLag2 = 0;
	int p = Y.n_rows; int T = Y.n_cols; int n = Y.n_slices;
	arma::vec missing = arma::ones(T);	
	arma::mat COVZlag1 = arma::mat(p, p, arma::fill::zeros);
	arma::mat COVZlag2 = arma::mat(p, p, arma::fill::zeros);
	arma::mat Yk; arma::uvec slh;

	// covariance estimation plus d.o.f. calculation 	
	for (int k = 0; k < n; ++k) {
		missing = arma::vec(T, arma::fill::ones);
		Yk = Y.slice(k);
		for (int l = 0; l < T; ++l){
			// ensure to exclude missings
			slh = arma::find_nonfinite(Yk.col(l));
			if (slh.n_elem != 0){ missing(l) = 0; } 
		}
        COVZdofLag1 = COVZdofLag1 + arma::sum(missing.subvec(0,T-2) % missing.subvec(1,T-1)); 		
        COVZdofLag2 = COVZdofLag2 + arma::sum(missing.subvec(0,T-3) % missing.subvec(2,T-1)); 		        
		Yk.elem( arma::find_nonfinite(Yk) ).zeros();
		COVZlag1 = COVZlag1 + Yk.submat(0, 1, p-1, T-1) * arma::trans(Yk.submat(0, 0, p-1, T-2));
		COVZlag2 = COVZlag2 + Yk.submat(0, 2, p-1, T-1) * arma::trans(Yk.submat(0, 0, p-1, T-3));		
	}
	COVZlag1 = COVZlag1 / COVZdofLag1;
	COVZlag2 = COVZlag2 / COVZdofLag2;
	return arma::join_rows(COVZlag1, COVZlag2);
}

inline arma::mat armaVAR2_VARYhat(arma::cube& Y, bool efficient){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Variance estimation of the multivariate time series data.
	/////////////////////////////////////////////////////////////////////////////////////////
 
	// variable declaration and initiation 
	int VARYdof = 0; int p = Y.n_rows; int T = Y.n_cols; int n = Y.n_slices;
	arma::vec missing; arma::mat VARY = arma::mat(2*p, 2*p, arma::fill::zeros);
	arma::mat Yk; arma::uvec slh;

	for (int k = 0; k < n; ++k) {
		Yk = Y.slice(k);		
		missing = arma::vec(T, arma::fill::ones);		
		if (efficient){
			// efficient variance calculation (uses data optimally)
			missing = arma::vec(T, arma::fill::ones);					
			for (int l = 0; l < T; ++l){
				// ensure to exclude missings
				slh = arma::find_nonfinite(Yk.col(l));
				if (slh.n_elem != 0){ missing(l) = 0; } 
			}
			Yk.elem( arma::find_nonfinite(Yk) ).zeros();
			VARYdof = VARYdof + arma::sum(missing.subvec(0,T-2) % missing.subvec(1,T-1)); 											
			VARY = VARY + arma::symmatl(arma::join_cols(Yk.submat(0,1,p-1,T-1), Yk.submat(0,0,p-1,T-2)) * arma::trans(arma::join_cols(Yk.submat(0,1,p-1,T-1), Yk.submat(0,0,p-1,T-2))));						
		} else {
			// "inefficient" variance calculation (in line with derivation of ML estimator of A)
			Yk.shed_col(0);			
			missing = arma::vec(T-1, arma::fill::ones);					
			for (int l = 0; l < (T-1); ++l){
				// ensure to exclude missings
				arma::uvec slh = arma::find_nonfinite(Yk.col(l));
				if (slh.n_elem > 0){ VARYdof = VARYdof + 1; }
			}
			Yk.elem( arma::find_nonfinite(Yk) ).zeros();
			VARYdof = VARYdof + arma::sum(missing.subvec(0,T-3) % missing.subvec(1,T-2)); 								
			VARY = VARY + arma::symmatl(arma::join_cols(Yk.submat(0,1,p-1,T-2), Yk.submat(0,0,p-1,T-3)) * arma::trans(arma::join_cols(Yk.submat(0,1,p-1,T-2), Yk.submat(0,0,p-1,T-3))));			
		}
	}
	return VARY / VARYdof;
}

inline arma::mat armaVAR2_Shat_ML(arma::cube& Y, arma::mat A1, arma::mat A2){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Least squares estimation of the regression parameter A of a VAR(1) model
	/////////////////////////////////////////////////////////////////////////////////////////

	// define relevant dimensions 
	int pY = A1.n_rows; int T = Y.n_cols; int n = Y.n_slices;

	// define dof and S
	arma::mat S = arma::symmatl(arma::zeros(pY, pY));
	int dof = 0;

	// evaluate dof and S
	arma::mat Yi;	
	for (int i = 0; i < n; ++i) {
		Yi = Y.slice(i);
		Yi = Yi.submat(0, 2, pY-1, T-1) - A1 * Yi.submat(0, 1, pY-1, T-2) - A2 * Yi.submat(0, 0, pY-1, T-3);
		arma::uvec ids = arma::find_finite(arma::sum(Yi));
		Yi = Yi.cols(ids);
		S = arma::symmatl(S) + arma::symmatl(Yi * arma::trans(Yi));
		dof = dof + ids.size();
	}
	return S / dof;
}

inline arma::mat armaVAR2_Ahat_ridgeSS(arma::mat& VARY, arma::mat& COVY, const double lambdaA1, const double lambdaA2, const arma::mat& targetA1, const arma::mat& targetA2){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Least squares estimation of the regression parameter A of a VAR(1) model
	/////////////////////////////////////////////////////////////////////////////////////////

	// estimate A by SS minimization
	arma::vec lambda = arma::vec(VARY.n_rows);
	lambda.fill(lambdaA2);
	lambda.subvec(0, COVY.n_rows-1) += lambdaA1 - lambdaA2;
	VARY.diag() += lambda;
	return (arma::join_rows(targetA1, targetA2) + COVY) * arma::inv_sympd(VARY);
}

inline arma::mat armaVAR2_Ahat_ridgeMLapprox(arma::mat& P, arma::mat& COVY, arma::mat& eigvecVARY, arma::vec eigvalVARY, const double lambdaA1, const double lambdaA2, const arma::mat& targetA1, const arma::mat& targetA2){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Ridge ML estimate the regression coefficient matrix A of a VAR model
	/////////////////////////////////////////////////////////////////////////////////////////

	// declare variables 
	int pY = P.n_rows;

	// eigendecomposition of error precision P
	arma::vec eigvalP; arma::mat eigvecP;
	arma::eig_sym(eigvalP, eigvecP, P, "dc");
	eigvalP = arma::real(eigvalP);

	// change-of-variable for penalty parameters
	double alpha = (lambdaA1 + lambdaA2) / 2;
	double beta = (lambdaA1 - lambdaA2) / 2;

	// estimate C with full support
	arma::mat AhatCorrection = arma::mat(pY, 2*pY);
	AhatCorrection.fill(alpha);  
	AhatCorrection += eigvalP * arma::trans(eigvalVARY);
	arma::mat AhatApprox = eigvecP * ((arma::trans(eigvecP) * (arma::join_rows(targetA1, targetA2) + P * COVY) * eigvecVARY) % 
				(1/(AhatCorrection) - beta / arma::square(AhatCorrection) + (beta * beta) / arma::pow(AhatCorrection, 3)) ) * arma::trans(eigvecVARY);
	AhatCorrection = 2 * eigvecP * ((arma::trans(eigvecP) * (targetA2 + P * COVY.cols(pY, 2*pY-1)) * eigvecVARY.rows(pY,2*pY-1)) % 
				(beta / arma::square(AhatCorrection))) * arma::trans(eigvecVARY.rows(pY, 2*pY-1));
	AhatApprox.cols(pY,2*pY-1) += AhatCorrection;
	return AhatApprox;
}

inline arma::mat armaVAR2_Ahat_zeros(arma::mat P, arma::mat& COVY, arma::mat& eigvecVARY, const arma::vec eigvalVARY, const double lambdaA1, const double lambdaA2, const arma::mat& targetA1, const arma::mat& targetA2, std::string fitA12, const arma::uvec& zerosRa, const arma::uvec& zerosCa, const arma::uvec& zerosRb, const arma::uvec& zerosCb, std::string zerosA1fit, std::string zerosA2fit){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Ridge ML estimate the regression coefficient matrix A of a VAR model with known support
	/////////////////////////////////////////////////////////////////////////////////////////

	// declare variables 
	int pY = P.n_rows;
	arma::mat Ahat = arma::mat(pY, 2*pY);

	// declare eigendecomposition of VARZ, blockwise	
	arma::mat VARY = eigvecVARY * arma::diagmat(eigvalVARY) * arma::trans(eigvecVARY);
	arma::mat eigvecVARYa1; arma::vec eigvalVARYa1;
	arma::eig_sym(eigvalVARYa1, eigvecVARYa1, VARY.submat(0, 0, pY-1, pY-1), "dc");
	arma::vec eigvalVARYa2; arma::mat eigvecVARYa2;
	arma::eig_sym(eigvalVARYa2, eigvecVARYa2, VARY.submat(pY, pY, 2*pY-1, 2*pY-1), "dc");

	// estimate C
	if (fitA12 == "ss"){
		Ahat = armaVAR2_Ahat_ridgeSS(VARY, COVY, lambdaA1, lambdaA2, targetA1, targetA2);
	}    
	if (fitA12 == "ml"){   
		Ahat = armaVAR2_Ahat_ridgeMLapprox(P, COVY, eigvecVARY, eigvalVARY, lambdaA1, lambdaA2, targetA1, targetA2);
	}

	// apply corrections (dependening on proportion of zeros).
	arma::vec eigvalP; arma::mat eigvecP;    
	if (zerosA1fit == "dense" || zerosA2fit == "dense"){      
	    arma::eig_sym(eigvalP, eigvecP, P, "dc");
	}
	if (zerosA1fit == "dense" && zerosRa.n_elem > 0){
		Ahat.submat(0, 0, pY-1, pY-1) = armaVARzerosCorrection_dense(Ahat.submat(0, 0, pY-1, pY-1), pY, pY, eigvecP, eigvalP, eigvecVARYa1, eigvalVARYa1, lambdaA1, zerosRa, zerosCa);
	}
	if (zerosA2fit == "dense" && zerosRb.n_elem > 0){
		Ahat.submat(0, pY, pY-1, 2*pY-1) = armaVARzerosCorrection_dense(Ahat.submat(0, pY, pY-1, 2*pY-1), pY, pY, eigvecP, eigvalP, eigvecVARYa2, eigvalVARYa2, lambdaA2, zerosRb, zerosCb);
	}
 
	if (zerosA1fit == "sparse" && zerosRa.n_elem > 0){
		Ahat.submat(0, 0, pY-1, pY-1) = armaVARzerosCorrection_sparse(Ahat.submat(0, 0, pY-1, pY-1), pY, pY, P, eigvecVARYa1, eigvalVARYa1, lambdaA1, zerosRa, zerosCa);
	}
	if (zerosA2fit == "sparse" && zerosRb.n_elem > 0){
		Ahat.submat(0, pY, pY-1, 2*pY-1) = armaVARzerosCorrection_sparse(Ahat.submat(0, pY, pY-1, 2*pY-1), pY, pY, P, eigvecVARYa2, eigvalVARYa2, lambdaA2, zerosRb, zerosCb);
	}

	return Ahat;
}

		
// [[Rcpp::export(".armaVAR2_ridgeML")]]
Rcpp::List armaVAR2_ridgeML(Rcpp::NumericVector& Yraw, const double lambdaA1, const double lambdaA2, const double lambdaP, const arma::mat& targetA1, const arma::mat& targetA2, arma::mat& targetP, std::string targetPtype, std::string fitA12, arma::mat unbalanced, bool diagP, bool efficient, const int nInit, const double minSuccDiff){

	// set profiles of missing (time, sample)-points to missing
	arma::cube Y;
	if (unbalanced.n_rows == 0){ Y = armaVAR_array2cube_withoutMissing(Yraw); }
	if (unbalanced.n_rows > 0){ Y = armaVAR_array2cube_withMissing(Yraw, arma::conv_to<arma::uvec >::from(unbalanced.col(0)), arma::conv_to<arma::uvec >::from(unbalanced.col(1))); }

	// dimension
	int pY = Y.n_rows;

	// estimate A by SS minimization
	arma::mat VARY = armaVAR2_VARYhat(Y, efficient);
	arma::mat COVY = armaVAR2_COVYhat(Y);
	arma::mat Ahat = armaVAR2_Ahat_ridgeSS(VARY, COVY, lambdaA1, lambdaA2, targetA1, targetA2);

	// calculate Se
	arma::mat Se = armaVAR2_Shat_ML(Y, Ahat.submat(0,0,pY-1,pY-1), Ahat.submat(0,pY,pY-1,2*pY-1));

	// ridge ML estimation of Se
	arma::mat Phat;
	if (targetPtype != "none"){ targetP = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); }
	if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
	if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }
	
	if (fitA12 == "ml"){
		///////////////////////////////////////////////////////////////////////////////
		// estimate parameters by ML, using the SS estimates as initials
		///////////////////////////////////////////////////////////////////////////////
 
		// eigen-decomposition of VARY
		arma::vec eigvalsVARY;
		arma::mat eigvecsVARY;
		arma::eig_sym(eigvalsVARY, eigvecsVARY, VARY, "dc");
    
		// declare storage variables
		arma::mat Aprev;
		arma::mat Pprev;
    
		for (int i = 0; i < nInit; ++i){
			// store latest estimates
			Aprev = Ahat;
			Pprev = Phat;

			// estimate A
			Ahat = armaVAR2_Ahat_ridgeMLapprox(Phat, COVY, eigvecsVARY, eigvalsVARY, lambdaA1, lambdaA2, targetA1, targetA2);
    
			// calculate Se
			arma::mat Se = armaVAR2_Shat_ML(Y, Ahat.submat(0, 0, pY-1, pY-1), Ahat.submat(0, pY, pY-1, 2*pY-1));

			// ridge ML estimation of precision
			if (targetPtype != "none"){ targetP = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); }
			if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
			if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }

			// assess convergence
			if (armaVAR2_convergenceEvaluation(Ahat, Aprev, Phat, Pprev) < minSuccDiff){ break; }
		}
	}
    
	// evaluate likelihood
	double val;
	double sign;
	arma::log_det(val, sign, Phat);
	double LL = (Y.n_cols - 1) * Y.n_slices * (val - arma::accu(Se % Phat)) / 2;

	return Rcpp::List::create(Rcpp::Named("A") = Ahat, Rcpp::Named("P") = Phat, Rcpp::Named("LL") = LL);
}


// [[Rcpp::export(".armaVAR2_ridgeML_zerosA")]]
Rcpp::List armaVAR2_ridgeML_zerosA(Rcpp::NumericVector& Yraw, const double lambdaA1, const double lambdaA2, const double lambdaP, const arma::mat& targetA1, const arma::mat& targetA2, arma::mat& targetP, std::string targetPtype, std::string fitA12, arma::mat unbalanced, bool diagP, bool efficient, const int nInit, const double minSuccDiff, const arma::uvec& zerosRa, const arma::uvec& zerosCa, const arma::uvec& zerosRb, const arma::uvec& zerosCb, std::string zerosA1fit, std::string zerosA2fit){

	// set profiles of missing (time, sample)-points to missing
	arma::cube Y;
	if (unbalanced.n_rows == 0){ Y = armaVAR_array2cube_withoutMissing(Yraw); } 
	if (unbalanced.n_rows > 0){ Y = armaVAR_array2cube_withMissing(Yraw, arma::conv_to<arma::uvec >::from(unbalanced.col(0)), arma::conv_to<arma::uvec >::from(unbalanced.col(1))); }

	// number of variates in Y and X
	const int pY = targetA1.n_rows; const int pX = targetA2.n_cols;

	// estimate A by SS minimization
	arma::mat VARY = armaVAR2_VARYhat(Y, efficient);
	arma::mat COVY = armaVAR2_COVYhat(Y);

	// declare eigendecomposition of VARZ
	arma::vec eigvalVARY; arma::mat eigvecVARY;
	arma::eig_sym(eigvalVARY, eigvecVARY, VARY, "dc");

	// estimate C by SS minimization	
	arma::mat Ahat = armaVAR2_Ahat_zeros(arma::eye<arma::mat>(pY,pY), COVY, eigvecVARY, eigvalVARY, lambdaA1, lambdaA2, targetA1, targetA2, fitA12, zerosRa, zerosCa, zerosRb, zerosCb, zerosA1fit, zerosA2fit);

	// calculate Se
	arma::mat Se = armaVAR2_Shat_ML(Y, Ahat.cols(0, pY-1), Ahat.cols(pY, pY+pX-1));

	// ridge ML estimation of Se
	arma::mat Phat;
	if (targetPtype != "none"){ targetP = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); }
	if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
	if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }

	if (fitA12 == "ml"){
		///////////////////////////////////////////////////////////////////////////////
		// estimate parameters by ML, using the SS estimates as initials
		///////////////////////////////////////////////////////////////////////////////

		// declare storage variables
		arma::mat Aprev;
		arma::mat Pprev;

		for (int i = 0; i < nInit; ++i){
			// store latest estimates
			Aprev = Ahat;
			Pprev = Phat;

			// estimate C
			Ahat = armaVAR2_Ahat_zeros(Phat, COVY, eigvecVARY, eigvalVARY, lambdaA1, lambdaA2, targetA1, targetA2, fitA12, zerosRa, zerosCa, zerosRb, zerosCb, zerosA1fit, zerosA2fit);

			// calculate Se
			arma::mat Se = armaVAR2_Shat_ML(Y, Ahat.cols(0, targetA1.n_cols-1), Ahat.cols(targetA1.n_cols, targetA1.n_cols + targetA2.n_cols-1));

			// construct target precision (if not provided)
			if (targetPtype != "none"){ targetP = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); }

			// ridge ML estimation of precision
			if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
			if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }

			// assess convergence
			if (armaVAR2_convergenceEvaluation(Ahat, Aprev, Phat, Pprev) < minSuccDiff){ break; }
		}
	}
    
	// evaluate likelihood
	double val;
	double sign;
	arma::log_det(val, sign, Phat);
	double LL = (Y.n_cols - 1) * Y.n_slices * (val - arma::accu(Se % Phat)) / 2;

	return Rcpp::List::create(Rcpp::Named("A") = Ahat, Rcpp::Named("P") = Phat, Rcpp::Named("LL")=LL);
}

// [[Rcpp::export(".armaVAR2_convergenceEvaluation")]]
double armaVAR2_convergenceEvaluation_forR(arma::mat& Ahat, arma::mat& Aprev, arma::mat& Phat, arma::mat& Pprev){
	// see armaVAR2_convergenceEvaluation
	return armaVAR2_convergenceEvaluation(Ahat, Aprev, Phat, Pprev);
}

// [[Rcpp::export(".armaVAR2_COVYhat")]]
arma::mat armaVARY_COVYhat_forR(Rcpp::NumericVector& Yraw){
	// Covariance estimation of the VAR(2) model from multivariate time series data.
	// First reformatting of data to an arma::cube format
	arma::cube Y = armaVAR_array2cube_withoutMissing(Yraw); 
	return armaVAR2_COVYhat(Y);
}

// [[Rcpp::export(".armaVAR2_VARYhat")]]
arma::mat armaVAR2_VARYhat_forR(Rcpp::NumericVector& Yraw, bool efficient){
	// Variance estimation of the VAR(2) model from multivariate time series data.
	// First reformatting of data to an arma::cube format
	arma::cube Y = armaVAR_array2cube_withoutMissing(Yraw); 
	return armaVAR2_VARYhat(Y, efficient);
}

// [[Rcpp::export(".armaVAR2_Ahat_ridgeSS")]]
arma::mat armaVAR2_Ahat_ridgeSS_forR(arma::mat& COVY, arma::mat& VARY, const double lambdaA1, const double lambdaA2, const arma::mat& targetA1, const arma::mat& targetA2){
	return armaVAR2_Ahat_ridgeSS(VARY, COVY, lambdaA1, lambdaA2, targetA1, targetA2);
}

// [[Rcpp::export(".armaVAR2_Ahat_ridgeML")]]
arma::mat armaVAR2_Ahat_ridgeML_forR(arma::mat& P, arma::mat& COVY, arma::mat& eigvecVARY, const arma::vec eigvalVARY, const double lambdaA1, const double lambdaA2, const arma::mat& targetA1, const arma::mat& targetA2){
	return armaVAR2_Ahat_ridgeMLapprox(P, COVY, eigvecVARY, eigvalVARY, lambdaA1, lambdaA2, targetA1, targetA2);
}

// [[Rcpp::export(".armaVAR2_Shat_ML")]]
arma::mat armaVAR2_Shat_ML_forR(Rcpp::NumericVector& Yraw, const arma::mat& A1, const arma::mat& A2){
	// Least squares estimation of the regression parameter A of a VAR(1) model
	// First reformatting of data to an arma::cube format
	arma::cube Y = armaVAR_array2cube_withoutMissing(Yraw); 
	return armaVAR2_Shat_ML(Y, A1, A2);
}

// [[Rcpp::export(".armaVAR2_Ahat_zeros")]]
arma::mat armaVAR2_Ahat_zeros_forR(arma::mat& P, arma::mat& COVY, arma::mat& eigvecVARY, arma::vec eigvalVARY, const double lambdaA1, const double lambdaA2, const arma::mat& targetA1, const arma::mat& targetA2, std::string fitA12, const arma::uvec& zerosRa, const arma::uvec& zerosCa, const arma::uvec& zerosRb, const arma::uvec& zerosCb, std::string zerosA1fit, std::string zerosA2fit){
	return armaVAR2_Ahat_zeros(P, COVY, eigvecVARY, eigvalVARY, lambdaA1, lambdaA2, targetA1, targetA2, fitA12, zerosRa, zerosCa, zerosRb, zerosCb, zerosA1fit, zerosA2fit);    
}

// [[Rcpp::export(".armaVAR2_loglik_LOOCVinternal")]]
double armaVAR2_loglikLOOCVinternal_forR(arma::vec Yt2, arma::vec Yt1, arma::vec Yt0, arma::mat& A1, arma::mat& A2, arma::mat& P){

	double LL = 0;
	int dofs = 0;
	double logDetP;
	double detPsign;
	arma::log_det(logDetP, detPsign, P);
		
	Yt2 = Yt2 - A1 * Yt1 - A2 * Yt0;
	return - arma::trace(Yt2.t() * P * Yt2) / 2 + logDetP / 2;
}


