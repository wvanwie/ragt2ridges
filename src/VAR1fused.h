// Rcpp::Rcout << "wessel is gek 1" << std::endl;					

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline double armaVAR1fused_convergenceEvaluation(arma::mat& Ahat, arma::mat& Aprev, arma::mat& Phat, arma::mat& Pprev){
	/////////////////////////////////////////////////////////////////////////////////////////
	// assess maximum difference between the current and previous estimates
	/////////////////////////////////////////////////////////////////////////////////////////

	// maximum of difference
	return (abs(arma::join_cols(Ahat, Phat) - arma::join_cols(Aprev, Pprev))).max();
}

inline double armaVAR1_convergenceEvaluationA(arma::mat& Ahat, arma::mat& Aprev){
	/////////////////////////////////////////////////////////////////////////////////////////
	// assess maximum difference between the current and previous estimates
	/////////////////////////////////////////////////////////////////////////////////////////

	// extract element-wise maximum of both matrices
	return abs(Ahat - Aprev).max();
}

inline arma::mat armaVAR1fused_Shat_ML(arma::cube& Y, arma::mat& As, arma::ivec grp){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Least squares estimation of the regression parameter A of a VAR(1) model
	/////////////////////////////////////////////////////////////////////////////////////////

	// covariance estimation 
	int p = Y.n_rows; int T = Y.n_cols; int n = Y.n_slices;

	// define dof and S
	arma::mat S = arma::symmatl(arma::zeros(p, p));
	int dof = 0; arma::uvec ids; arma::mat Yi;	
	for (int i = 0; i < n; ++i) {
	    arma::mat A = As.submat(p*grp(i), 0, p*(grp(i)+1)-1, p-1);
		Yi = Y.slice(i);
		Yi = Yi.submat(0, 1, p-1, T-1) - A * Yi.submat(0, 0, p-1, T-2);
		ids = arma::find_finite(arma::sum(Yi));
		Yi = Yi.cols(ids);
		S = arma::symmatl(S) + arma::symmatl(Yi * arma::trans(Yi));
		dof = dof + ids.size();
	}
	return S / dof;
}

// [[Rcpp::export(".armaVAR1fused_Shat_ML")]]
arma::mat armaVAR1fused_Shat_ML_forR(Rcpp::NumericVector& Yraw, arma::mat& As, arma::ivec grp){
	// Calculate the joint sample covariance matrix
	// First reformatting of data to an arma::cube format
	arma::cube Y = armaVAR_array2cube_withoutMissing(Yraw); 
	return armaVAR1fused_Shat_ML(Y, As, grp);
}

inline arma::mat armaVAR1f_Ahat_ridgeML(arma::mat& P, arma::mat COVY, const arma::mat eigvecVARY, const arma::vec eigvalVARY, const double lambdaA, arma::mat& targetA){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Ridge ML estimate the regression coefficient matrix A of a VAR model
	/////////////////////////////////////////////////////////////////////////////////////////

	// eigendecomposition of error precision P
	arma::vec eigvalP; arma::mat eigvecP;
	arma::eig_sym(eigvalP, eigvecP, P, "dc");
	eigvalP = arma::real(eigvalP);

	// evaluate in return the estimate of A
	return eigvecP * ((arma::trans(eigvecP) * (targetA + P * COVY) * eigvecVARY) / (eigvalP * arma::trans(eigvalVARY) + lambdaA)) * arma::trans(eigvecVARY);
}


inline arma::mat armaVAR1fused_Ahat_ridgeML_withoutZeros(arma::mat& As, arma::mat& P, arma::mat& COVYs, arma::mat& eigvecVARYs, arma::vec eigvalVARYs, const double lambdaA, const double lambdaF, const arma::mat& targetA, int nInitA, const double minSuccDiffA){
	// Fused ridge maximum likelihood estimation of the regression parameters A of several VAR(1) models with a common precision of the innovations

	// initiation of several variables
	int p = As.n_cols;
	int nGroups = As.n_rows / p;
	arma::mat targetAf;
	arma::mat Aprevs = As;

	// iterative estimation of the As
	for (int itA = 0; itA < nInitA; ++itA){
		for (int g = 0; g < nGroups; ++g){

			// modify targetA for estimation of A of group g
			targetAf = targetA;             
			for (int g1 = 0; g1 < nGroups; ++g1){
				if (g1 != g){
					targetAf = targetAf + lambdaF * As.submat(p*g1 ,0, p*(g1+1)-1, p-1);
				}
			}
            
			// update A estimate of group g    
			As.submat(p*g ,0, p*(g+1)-1, p-1) = armaVAR1f_Ahat_ridgeML(P, COVYs.submat(p*g ,0, p*(g+1)-1, p-1), eigvecVARYs.submat(p*g ,0, p*(g+1)-1, p-1), 
                                                        eigvalVARYs.subvec(p*g, p*(g+1)-1), lambdaA+(nGroups-1)*lambdaF, targetAf);
		}            
		// convergence criterion
		if (armaVAR1_convergenceEvaluationA(As, Aprevs) < minSuccDiffA){ break; }
	}
	return As;
}

inline arma::mat armaVAR1f_Ahat_zeros(const arma::mat& P, arma::mat COVY, const arma::mat eigvecVARY, const arma::vec eigvalVARY, const double lambdaA, const arma::mat& targetA, std::string fitA, const arma::uvec zerosR, const arma::uvec zerosC, std::string zerosAfit){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Ridge ML estimate the regression coefficient matrix A of a VAR model with known support
	/////////////////////////////////////////////////////////////////////////////////////////

	// declare variables 
	int pY = COVY.n_rows; arma::vec eigvalP; arma::mat eigvecP; arma::mat Ahat;

	if (fitA == "ss"){
		// eigendecomposition of error precision P
		eigvalP = arma::ones<arma::vec>(pY);
		eigvecP = P;
    	 
		// estimate A with full support
		Ahat = (((targetA + COVY) * eigvecVARY) / (eigvalP * arma::trans(eigvalVARY) + lambdaA)) * arma::trans(eigvecVARY);
	}
	if (fitA == "ml"){
		// eigendecomposition of error precision P
		arma::eig_sym(eigvalP, eigvecP, P, "dc");
		eigvalP = arma::real(eigvalP);
 
		// estimate A with full support
		Ahat = eigvecP * ((arma::trans(eigvecP) * (targetA + P * COVY) * eigvecVARY) / (eigvalP * arma::trans(eigvalVARY) + lambdaA)) * arma::trans(eigvecVARY);
	}

	if (zerosR.n_elem > 0){    
		if (zerosAfit == "dense"){
			Ahat = armaVARzerosCorrection_dense(Ahat, pY, pY, eigvecP, eigvalP, eigvecVARY, eigvalVARY, lambdaA, zerosR, zerosC);
		} 
		if (zerosAfit == "sparse"){
			Ahat = armaVARzerosCorrection_sparse(Ahat, pY, pY, P, eigvecVARY, eigvalVARY, lambdaA, zerosR, zerosC);
		}
	}

	return Ahat;
}

inline arma::mat armaVAR1fused_Ahat(arma::mat& As, arma::mat P, arma::mat& COVYs, arma::mat eigvecVARYs, arma::vec eigvalVARYs, const double lambdaA, const double lambdaF, const arma::mat& targetA, std::string fitA, const arma::uvec& zerosR, const arma::uvec& zerosC, std::string zerosAfit, int nInitA, const double minSuccDiffA){

	// Fused ridge maximum likelihood estimation of the regression parameters A of several VAR(1) models with a common precision of the innovations

	// initiation of several variables
	int p = As.n_cols;
	int nGroups = As.n_rows / p;
	arma::mat targetAf;

	// iterative estimation of the As
	for (int itA = 0; itA < nInitA; ++itA){
		arma::mat Aprevs = As;	
		for (int g = 0; g < nGroups; ++g){
			// modify targetA for estimation of A of group g
			targetAf = targetA;             
			for (int g1 = 0; g1 < nGroups; ++g1){
				if (g1 != g){
					targetAf = targetAf + lambdaF * As.submat(p*g1 ,0, p*(g1+1)-1, p-1);
				}
			}
            
			// update A estimate of group g    
			As.submat(p*g ,0, p*(g+1)-1, p-1) = armaVAR1f_Ahat_zeros(P, COVYs.submat(p*g ,0, p*(g+1)-1, p-1), eigvecVARYs.submat(p*g ,0, p*(g+1)-1, p-1), 
                                                        eigvalVARYs.subvec(p*g, p*(g+1)-1), lambdaA+(nGroups-1)*lambdaF, targetAf, fitA, zerosR, zerosC, zerosAfit);          
		}            

		// convergence criterion
		if (armaVAR1_convergenceEvaluationA(As, Aprevs) < minSuccDiffA){ break; }
	}
	return As;
}

// [[Rcpp::export(".armaVAR1fused_Ahat")]]
arma::mat armaVAR1fused_Ahat_forR(arma::mat& As, arma::mat P, arma::mat& COVYs, arma::mat eigvecVARYs, arma::vec eigvalVARYs, const double lambdaA, const double lambdaF, const arma::mat& targetA, std::string fitA, const arma::uvec& zerosR, const arma::uvec& zerosC, std::string zerosAfit, int nInitA, const double minSuccDiffA){
	// Fused ridge maximum likelihood estimation of the regression parameters A of several VAR(1) models with a common precision of the innovations
	return armaVAR1fused_Ahat(As, P, COVYs, eigvecVARYs, eigvalVARYs, lambdaA, lambdaF, targetA, fitA, zerosR, zerosC, zerosAfit, nInitA, minSuccDiffA);
}

inline Rcpp::List armaEigenDecomp_stackedCovariances(arma::mat& mats){
	/////////////////////////////////////////////////////////////////////////////////////////
	// eigendecomposition of multiple covariance matrices: stored in stacked form
	/////////////////////////////////////////////////////////////////////////////////////////
	arma::vec eigvals = arma::zeros(mats.n_rows);
	arma::mat eigvecs = arma::zeros(mats.n_rows, mats.n_cols);
	arma::vec eigvalsSLH;
	arma::mat eigvecsSLH;
	
	// initiation of several variables
	int p = mats.n_cols;
	int nGroups = mats.n_rows / p;

	for (unsigned int k = 0; k < nGroups; ++k) {
		arma::eig_sym(eigvalsSLH, eigvecsSLH, mats.submat(k*p, 0, (k+1)*p-1, p-1), "dc");
		eigvecs.submat(k*p, 0, (k+1)*p-1, p-1) = eigvecsSLH;
		eigvals.subvec(k*p, (k+1)*p-1) = eigvalsSLH;
	}
	return Rcpp::List::create(Rcpp::Named("values") = eigvals, Rcpp::Named("vectors") = eigvecs);
}

// [[Rcpp::export(".armaEigenDecomp_stackedCovariances")]]
Rcpp::List armaEigenDecomp_stackedCovariances_forR(arma::mat& mats){
	// eigendecomposition of multiple covariance matrices: stored in stacked form
	return armaEigenDecomp_stackedCovariances(mats);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
WRAPPERS FUNCTIONs FOR THE VAR(1) MODEL WITH FUSED RIDGE ESTIMATION
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline arma::mat armaVAR1f_Ahat_ridgeSS(arma::mat VARY, arma::mat COVY, const double lambdaA, const arma::mat& targetA){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Least squares estimation of the regression parameter A of a VAR(1) model
	/////////////////////////////////////////////////////////////////////////////////////////

	// evaluate the estimate of A
	VARY.diag() += lambdaA;
	return (targetA + COVY) * arma::inv_sympd(VARY);
}


inline Rcpp::List armaVAR1fused_ridgeML(Rcpp::NumericVector& Yraw, arma::ivec id, const double lambdaA, const double lambdaF, const double lambdaP, const arma::mat& targetA, const arma::mat& targetP, std::string targetPtype, std::string fitA, arma::mat unbalanced, bool diagP, bool efficient, const arma::uvec& zerosR, const arma::uvec& zerosC, std::string zerosAfit, const int nInit, const int nInitA, const double minSuccDiff, const double minSuccDiffA){

	// set profiles of missing (time, sample)-points to missing
	arma::cube Y;
	if (unbalanced.n_rows == 0){ Y = armaVAR_array2cube_withoutMissing(Yraw); } 
	if (unbalanced.n_rows > 0){ Y = armaVAR_array2cube_withMissing(Yraw, arma::conv_to<arma::uvec >::from(unbalanced.col(0)), arma::conv_to<arma::uvec >::from(unbalanced.col(1))); }

	// extract dimensions
	int p = Y.n_rows;

	// estimate As by SS minimization and store as a long matrix
	arma::ivec idUniq = arma::unique(id);			    
	arma::mat VARYs = arma::mat(p*idUniq.n_elem,p); arma::mat COVYs = arma::mat(p*idUniq.n_elem,p); arma::mat Ahats = arma::mat(p*idUniq.n_elem,p); arma::uvec sliceID;
	for (int g = 0; g < idUniq.n_elem; ++g){	
		sliceID = arma::find(id == g);
		VARYs.submat(p*g,0,p*(g+1)-1,p-1) = armaVAR1_VARYhat(Y.slices(sliceID.min(), sliceID.max()), efficient, unbalanced);
		COVYs.submat(p*g,0,p*(g+1)-1,p-1) = armaVAR1_COVYhat(Y.slices(sliceID.min(), sliceID.max()));
		Ahats.submat(p*g,0,p*(g+1)-1,p-1) = armaVAR1f_Ahat_ridgeSS(VARYs.submat(p*g,0,p*(g+1)-1,p-1), COVYs.submat(p*g,0,p*(g+1)-1,p-1), lambdaA, targetA);
	}

	// eigenvalue decompositions of groups-wise process variance estimates    
	Rcpp::List eigDecomps = armaEigenDecomp_stackedCovariances(VARYs);

	// estimate As in fused sum-of-squares fashion
	if (fitA == "ss"){ 
		Ahats = armaVAR1fused_Ahat(Ahats, arma::eye(p,p), COVYs, eigDecomps[1], eigDecomps[0], lambdaA, lambdaF, targetA, fitA, zerosR, zerosC, zerosAfit, nInitA, minSuccDiffA);
	}

	// calculate Se
	arma::mat Se = armaVAR1fused_Shat_ML(Y, Ahats, id);
    
	// ridge ML estimation of precision
	arma::mat Phat;
	if (targetPtype == "none"){ 
		if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
		if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }
	}
	if (targetPtype != "none"){ 
		arma::mat targetPnew = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); 
		if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetPnew), lambdaP); }
		if (!diagP){ Phat = armaRidgeP(Se, targetPnew, lambdaP); }
	}
 
	if (fitA == "ml"){ 
		///////////////////////////////////////////////////////////////////////////////
		// estimate parameters by ML, using the SS estimates are initials
		///////////////////////////////////////////////////////////////////////////////

		// declare storage variables
		arma::mat Aprevs; arma::mat Pprev;
    
		for (int i = 0; i < nInit; ++i){
			// store latest estimates
			Aprevs = Ahats; Pprev = Phat;

			// estimate As in fused fashion
			Ahats = armaVAR1fused_Ahat(Ahats, Phat, COVYs, eigDecomps[1], eigDecomps[0], lambdaA, lambdaF, targetA, fitA, zerosR, zerosC, zerosAfit, nInitA, minSuccDiffA);
                
			// calculate Se
			Se = armaVAR1fused_Shat_ML(Y, Ahats, id);

			// ridge ML estimation of precision
			arma::mat Phat;
			if (targetPtype == "none"){ 
				if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
				if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP, 2); }
			}
			if (targetPtype != "none"){ 
				arma::mat targetPnew = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); 
				if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetPnew), lambdaP); }
				if (!diagP){ Phat = armaRidgeP(Se, targetPnew, lambdaP); }
			}
        
			// assess convergence
			if (armaVAR1fused_convergenceEvaluation(Ahats, Aprevs, Phat, Pprev) < minSuccDiff){ break; }
		}
	}
   
	/* 
	// evaluate likelihood
	double logDetP; double detPsign;
	arma::log_det(logDetP, detPsign, Phat);
	double LL = (Y.n_cols - 1) * Y.n_slices * (logDetP - arma::accu(Se % Phat)) / 2;
	*/

	return Rcpp::List::create(Rcpp::Named("As") = Ahats, Rcpp::Named("P") = Phat, Rcpp::Named("LL")=0);
}

// [[Rcpp::export(".armaVAR1fused_ridgeML")]]
Rcpp::List armaVAR1fused_ridgeML_forR(Rcpp::NumericVector& Yraw, arma::ivec id, const double lambdaA, const double lambdaF, const double lambdaP, const arma::mat& targetA, const arma::mat& targetP, std::string targetPtype, std::string fitA, arma::mat unbalanced, bool diagP, bool efficient, const arma::uvec& zerosR, const arma::uvec& zerosC, std::string zerosAfit, const int nInit, const int nInitA, const double minSuccDiff, const double minSuccDiffA){
	// eigendecomposition of multiple covariance matrices: stored in stacked form
	return armaVAR1fused_ridgeML(Yraw, id, lambdaA, lambdaF, lambdaP, targetA, targetP, targetPtype, fitA, unbalanced, diagP, efficient, zerosR, zerosC, zerosAfit, nInit, nInitA, minSuccDiff, 				minSuccDiffA);	
}

// [[Rcpp::export(".armaVAR1fused_convergenceEvaluation")]]
double armaVAR1fused_convergenceEvaluation_forR(arma::mat& Ahat, arma::mat& Aprev, arma::mat& Phat, arma::mat& Pprev){
	/////////////////////////////////////////////////////////////////////////////////////////
	// assess maximum difference between the current and previous estimates
	/////////////////////////////////////////////////////////////////////////////////////////

	// maximum of difference
	return armaVAR1fused_convergenceEvaluation(Ahat, Aprev, Phat, Pprev);
}

