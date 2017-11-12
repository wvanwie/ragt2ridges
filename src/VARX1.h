////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
FUCTIONS FOR THE ESTIMATION OF THE VARX(1) MODEL: FOR USE ON THE SEASIDE, FOLLOWED BY WRAPPERS FOR DEPORTATION TO R.
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double armaVARX1_convergenceEvaluation(arma::mat& Chat, arma::mat& Cprev, arma::mat& Phat, arma::mat& Pprev){
	/////////////////////////////////////////////////////////////////////////////////////////
	// assess maximum difference between the current and previous estimates
	/////////////////////////////////////////////////////////////////////////////////////////

	// extract element-wise maximum of both matrices
	arma::mat maxis = arma::max(arma::max(abs(Chat - Cprev), 1), arma::max(abs(Phat - Pprev), 1));
	return maxis.max();
}

inline arma::mat armaVARX1_COVZhat(arma::cube& Z, const int pY){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Covariance estimation of the multivariate time series data.
	/////////////////////////////////////////////////////////////////////////////////////////

	// variable declaration and initiation 
	int COVZdof = 0; int p = Z.n_rows; int T = Z.n_cols; int n = Z.n_slices;
	arma::vec missing = arma::ones(T);
	arma::mat COVZ = arma::mat(pY, p, arma::fill::zeros);
	arma::mat Zk; arma::uvec slh;
	
	// covariance estimation plus d.o.f. calculation 	
	for (int k = 0; k < n; ++k) {
		missing = arma::vec(T, arma::fill::ones);
		Zk = Z.slice(k);
		for (int l = 0; l < T; ++l){
			// ensure to exclude missings
			slh = arma::find_nonfinite(Zk.col(l));
			if (slh.n_elem != 0){ missing(l) = 0; } 
		}
        	COVZdof = COVZdof + arma::sum(missing.subvec(0,T-2) % missing.subvec(1,T-1)); 		
		Zk.elem( arma::find_nonfinite(Zk) ).zeros();
		COVZ = COVZ + Zk.submat(0, 1, pY-1, T-1) * arma::trans(Zk.submat(0, 0, p-1, T-2));
	}
	return COVZ / COVZdof;
}

inline arma::mat armaVARX1_Shat_ML(arma::cube& Z, arma::mat A, arma::mat B){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Least squares estimation of the regression parameter A of a VAR(1) model
	/////////////////////////////////////////////////////////////////////////////////////////

	// define relevant dimensions 
	int pY = A.n_rows; int pX = B.n_cols; int T = Z.n_cols; int n = Z.n_slices;

	// define dof and S
	arma::mat S = arma::symmatl(arma::zeros(pY, pY));
	int dof = 0;

	// evaluate dof and S
	arma::mat Zi;	
	for (int i = 0; i < n; ++i) {
		Zi = Z.slice(i);
		Zi = Zi.submat(0, 1, pY-1, T-1) - A * Zi.submat(0, 0, pY-1, T-2) - B * Zi.submat(pY, 0, pY + pX-1, T-2);
		arma::uvec ids = arma::find_finite(arma::sum(Zi));
		Zi = Zi.cols(ids);
		S = arma::symmatl(S) + arma::symmatl(Zi * arma::trans(Zi));
		dof = dof + ids.size();
	}
	return S / dof;
}

inline arma::mat armaVARX1_Chat_ridgeSS(arma::mat VARZ, arma::mat COVZ, const double lambdaA, const double lambdaB, const arma::mat targetA, const arma::mat targetB){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Least squares estimation of the regression parameter A of a VAR(1) model
	/////////////////////////////////////////////////////////////////////////////////////////

	// evaluate the estimate of A
	arma::vec lambda = arma::vec(VARZ.n_rows);
	lambda.fill(lambdaB); 
	lambda.subvec(0, COVZ.n_rows-1) += lambdaA - lambdaB; 
	VARZ.diag() += lambda;
 	return (arma::join_rows(targetA, targetB) + COVZ) * arma::inv_sympd(VARZ);

}

inline arma::mat armaVARX1_Chat_ridgeMLapprox(arma::mat P, arma::mat COVZ, arma::mat eigvecVARZ, const arma::vec eigvalVARZ, const double lambdaA, const double lambdaB, const arma::mat  targetA, const arma::mat targetB){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Ridge ML estimate the regression coefficient matrix A of a VAR model
	/////////////////////////////////////////////////////////////////////////////////////////

	// declare variables 
	int pX = targetB.n_cols; int pY = P.n_rows;

	// eigendecomposition of error precision P
	arma::vec eigvalP; arma::mat eigvecP;
	arma::eig_sym(eigvalP, eigvecP, P, "dc");
	eigvalP = arma::real(eigvalP);

	// change-of-variable for penalty parameters
	double alpha = (lambdaA + lambdaB) / 2;
	double beta = (lambdaA - lambdaB) / 2;

	// estimate C with full support
	arma::mat ChatCorrection = arma::mat(pY, pX+pY);
	ChatCorrection.fill(alpha);  
	ChatCorrection += eigvalP * arma::trans(eigvalVARZ);
	arma::mat ChatApprox = eigvecP * ((arma::trans(eigvecP) * (arma::join_rows(targetA, targetB) + P * COVZ) * eigvecVARZ) % 
				(1/(ChatCorrection) - beta / arma::square(ChatCorrection) + (beta * beta) / arma::pow(ChatCorrection, 3)) ) * arma::trans(eigvecVARZ);
	ChatCorrection = 2 * eigvecP * ((arma::trans(eigvecP) * (targetB + P * COVZ.cols(pY,pX+pY-1)) * eigvecVARZ.rows(pY,pX+pY-1)) % 
				(beta / arma::square(ChatCorrection))) * arma::trans(eigvecVARZ.rows(pY,pX+pY-1));
	ChatApprox.cols(pY,pX+pY-1) += ChatCorrection;
	return ChatApprox;
}


inline arma::mat armaVARX1_Chat_zeros(arma::mat& P, arma::mat& COVZ, arma::mat& eigvecVARZ, const arma::vec eigvalVARZ, const double lambdaA, const double lambdaB, const arma::mat& targetA, const arma::mat& targetB, std::string fitAB, const arma::uvec& zerosRa, const arma::uvec& zerosCa, const arma::uvec& zerosRb, const arma::uvec& zerosCb, std::string zerosAfit, std::string zerosBfit){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Ridge ML estimate the regression coefficient matrix A of a VAR model with known support
	/////////////////////////////////////////////////////////////////////////////////////////

	// declare variables 
	int pX = targetB.n_cols; int pY = P.n_rows;
	arma::mat Chat = arma::mat(pY, pY+pX);

	// declare eigendecomposition of VARZ, blockwise	
	arma::mat VARZ = eigvecVARZ * arma::diagmat(eigvalVARZ) * arma::trans(eigvecVARZ);
	arma::mat eigvecVARZa; arma::vec eigvalVARZa;
	arma::eig_sym(eigvalVARZa, eigvecVARZa, VARZ.submat(0, 0, pY-1, pY-1), "dc");
	arma::vec eigvalVARZb; arma::mat eigvecVARZb;
	arma::eig_sym(eigvalVARZb, eigvecVARZb, VARZ.submat(pY, pY, pY+pX-1, pY+pX-1), "dc");

	// estimate C
	if (fitAB == "ss"){
		Chat = armaVARX1_Chat_ridgeSS(VARZ, COVZ, lambdaA, lambdaB, targetA, targetB);
	}    
	if (fitAB == "ml"){   
		Chat = armaVARX1_Chat_ridgeMLapprox(P, COVZ, eigvecVARZ, eigvalVARZ, lambdaA, lambdaB, targetA, targetB);
	}

	// apply corrections (dependening on proportion of zeros).
	arma::vec eigvalP; arma::mat eigvecP;    
	if (zerosAfit == "dense" || zerosBfit == "dense"){      
	    arma::eig_sym(eigvalP, eigvecP, P, "dc");
	}
    if (zerosAfit == "dense" && zerosRa.n_elem > 0){
			Chat.submat(0, 0, pY-1, pY-1) = armaVARzerosCorrection_dense(Chat.submat(0, 0, pY-1, pY-1), pY, pY, eigvecP, eigvalP, eigvecVARZa, eigvalVARZa, lambdaA, zerosRa, zerosCa);
	}
	if (zerosBfit == "dense" && zerosRb.n_elem > 0){
		Chat.submat(0, pY, pY-1, pX+pY-1) = armaVARzerosCorrection_dense(Chat.submat(0, pY, pY-1, pY+pX-1), pY, pX, eigvecP, eigvalP, eigvecVARZb, eigvalVARZb, lambdaB, zerosRb, zerosCb);
	}
 
	if (zerosAfit == "sparse" && zerosRa.n_elem > 0){
		Chat.submat(0, 0, pY-1, pY-1) = armaVARzerosCorrection_sparse(Chat.submat(0, 0, pY-1, pY-1), pY, pY, P, eigvecVARZa, eigvalVARZa, lambdaA, zerosRa, zerosCa);
	}
	if (zerosBfit == "sparse" && zerosRb.n_elem > 0){
		Chat.submat(0, pY, pY-1, pX+pY-1) = armaVARzerosCorrection_sparse(Chat.submat(0, pY, pY-1, pY+pX-1), pY, pX, P, eigvecVARZb, eigvalVARZb, lambdaB, zerosRb, zerosCb);
	}

	return Chat;
}

// [[Rcpp::export(".armaVARX1_ridgeML")]]
Rcpp::List armaVARX1_ridgeML(Rcpp::NumericVector& Zraw, const double lambdaA, const double lambdaB, const double lambdaP, const arma::mat& targetA, const arma::mat& targetB, const arma::mat& targetP, std::string targetPtype, std::string fitAB, arma::mat unbalanced, bool diagP, bool efficient, const int nInit, const double minSuccDiff){

	// number of variates in Y and X
	const int pY = targetA.n_rows; 	const int pX = targetB.n_cols;

	// set profiles of missing (time, sample)-points to missing
	arma::cube Z;
	if (unbalanced.n_rows == 0){ Z = armaVAR_array2cube_withoutMissing(Zraw); } 
	if (unbalanced.n_rows > 0){ Z = armaVAR_array2cube_withMissing(Zraw, arma::conv_to<arma::uvec >::from(unbalanced.col(0)), arma::conv_to<arma::uvec >::from(unbalanced.col(1))); }

	// estimate A by SS minimization
	arma::mat VARZ = armaVAR1_VARYhat(Z, efficient, unbalanced);
	arma::mat COVZ = armaVARX1_COVZhat(Z, pY);
	arma::mat Chat = armaVARX1_Chat_ridgeSS(VARZ, COVZ, lambdaA, lambdaB, targetA, targetB);

	// calculate Se
	arma::mat Se = armaVARX1_Shat_ML(Z, Chat.cols(0, pY-1), Chat.cols(pY, pY+pX-1));

	// ridge ML estimation of Se
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

	if (fitAB == "ml"){
		///////////////////////////////////////////////////////////////////////////////
		// estimate parameters by ML, using the SS estimates as initials
		///////////////////////////////////////////////////////////////////////////////
    
		// declare eigendecomposition of VARZ
		arma::vec eigvalsVARZ;
		arma::mat eigvecsVARZ;
		arma::eig_sym(eigvalsVARZ, eigvecsVARZ, VARZ, "dc");

		// declare storage variables
		arma::mat Cprev; arma::mat Pprev;

		for (int i = 0; i < nInit; ++i){
			// store latest estimates
			Cprev = Chat; Pprev = Phat;

			// estimate A
			Chat = armaVARX1_Chat_ridgeMLapprox(Phat, COVZ, eigvecsVARZ, eigvalsVARZ, lambdaA, lambdaB, targetA, targetB);

			// calculate Se
			arma::mat Se = armaVARX1_Shat_ML(Z, Chat.cols(0, pY-1), Chat.cols(pY, pY+pX-1));

			// ridge ML estimation of Se
			if (targetPtype == "none"){ 
				if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
				if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }
			}
			if (targetPtype != "none"){ 
				arma::mat targetPnew = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); 
				if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetPnew), lambdaP); }
				if (!diagP){ Phat = armaRidgeP(Se, targetPnew, lambdaP); }
			}

			// assess convergence
			if (armaVARX1_convergenceEvaluation(Chat, Cprev, Phat, Pprev) < minSuccDiff){ break; }
		}
	}

	// evaluate likelihood
	double val;
	double sign;	
	arma::log_det(val, sign, Phat);
	double LL = (Z.n_cols - 1) * Z.n_slices * (val - arma::accu(Se % Phat)) / 2;

	return Rcpp::List::create(Rcpp::Named("A") = Chat.cols(0, pY-1), Rcpp::Named("B") = Chat.cols(pY, pY+pX-1), Rcpp::Named("P") = Phat, Rcpp::Named("LL")=LL);
}

// [[Rcpp::export(".armaVARX1_ridgeML_zerosC")]]
Rcpp::List armaVARX1_ridgeML_zerosC(Rcpp::NumericVector& Zraw, const double lambdaA, const double lambdaB, const double lambdaP, const arma::mat& targetA, const arma::mat& targetB, const arma::mat& targetP, std::string targetPtype, std::string fitAB, arma::mat unbalanced, bool diagP, bool efficient, const int nInit, const double minSuccDiff, const arma::uvec& zerosRa, const arma::uvec& zerosCa, const arma::uvec& zerosRb, const arma::uvec& zerosCb, std::string zerosAfit, std::string zerosBfit){

	// set profiles of missing (time, sample)-points to missing
	arma::cube Z;
	if (unbalanced.n_rows == 0){ Z = armaVAR_array2cube_withoutMissing(Zraw); } 
	if (unbalanced.n_rows > 0){ Z = armaVAR_array2cube_withMissing(Zraw, arma::conv_to<arma::uvec >::from(unbalanced.col(0)), arma::conv_to<arma::uvec >::from(unbalanced.col(1))); }

	// number of variates in Y and X
	const int pY = targetA.n_rows; 	const int pX = targetB.n_cols;

	// estimate A by SS minimization
	arma::mat VARZ = armaVAR1_VARYhat(Z, efficient, unbalanced);
	arma::mat COVZ = armaVARX1_COVZhat(Z, pY);

	// declare eigendecomposition of VARZ
	arma::vec eigvalVARZ; arma::mat eigvecVARZ;
	arma::eig_sym(eigvalVARZ, eigvecVARZ, VARZ, "dc");

	// estimate C by SS minimization
	arma::mat Ipp =	arma::eye<arma::mat>(pY,pY);
	arma::mat Chat = armaVARX1_Chat_zeros(Ipp, COVZ, eigvecVARZ, eigvalVARZ, lambdaA, lambdaB, targetA, targetB, fitAB, zerosRa, zerosCa, zerosRb, zerosCb, zerosAfit, zerosBfit);

	// calculate Se
	arma::mat Se = armaVARX1_Shat_ML(Z, Chat.cols(0, pY-1), Chat.cols(pY, pY+pX-1));

	// ridge ML estimation of Se
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

	if (fitAB == "ml"){
		///////////////////////////////////////////////////////////////////////////////
		// estimate parameters by ML, using the SS estimates as initials
		///////////////////////////////////////////////////////////////////////////////

		// declare storage variables
		arma::mat Cprev;
		arma::mat Pprev;

		for (int i = 0; i < nInit; ++i){
			// store latest estimates
			Cprev = Chat;
			Pprev = Phat;

			// estimate C
			Chat = armaVARX1_Chat_zeros(Phat, COVZ, eigvecVARZ, eigvalVARZ, lambdaA, lambdaB, targetA, targetB, fitAB, zerosRa, zerosCa, zerosRb, zerosCb, zerosAfit, zerosBfit);

			// calculate Se
			arma::mat Se = armaVARX1_Shat_ML(Z, Chat.cols(0, targetA.n_cols-1), Chat.cols(targetA.n_cols, targetA.n_cols + targetB.n_cols-1));

			// ridge ML estimation of Se
			if (targetPtype == "none"){ 
				if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
				if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }
			}
			if (targetPtype != "none"){ 
				arma::mat targetPnew = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); 
				if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetPnew), lambdaP); }
				if (!diagP){ Phat = armaRidgeP(Se, targetPnew, lambdaP); }
			}

			// assess convergence
			if (armaVARX1_convergenceEvaluation(Chat, Cprev, Phat, Pprev) < minSuccDiff){ break; }
		}
	}
    
	// evaluate likelihood
	double val;
	double sign;
	arma::log_det(val, sign, Phat);
	double LL = (Z.n_cols - 1) * Z.n_slices * (val - arma::accu(Se % Phat)) / 2;

	return Rcpp::List::create(Rcpp::Named("A") = Chat.cols(0, pY-1), Rcpp::Named("B") = Chat.cols(pY, pY+pX-1), Rcpp::Named("P") = Phat, Rcpp::Named("LL")=LL);
}

// [[Rcpp::export(".armaVARX1_Chat_zeros")]]
arma::mat armaVARX1_Chat_zeros_forR(arma::mat& P, arma::mat& COVZ, arma::mat& eigvecVARZ, const arma::vec eigvalVARZ, const double lambdaA, const double lambdaB, const arma::mat& targetA, const arma::mat& targetB, std::string fitAB, const arma::uvec& zerosRa, const arma::uvec& zerosCa, const arma::uvec& zerosRb, const arma::uvec& zerosCb, std::string zerosAfit, std::string zerosBfit){
	return armaVARX1_Chat_zeros(P, COVZ, eigvecVARZ, eigvalVARZ, lambdaA, lambdaB, targetA, targetB, fitAB, zerosRa, zerosCa, zerosRb, zerosCb, zerosAfit, zerosBfit);    
}

// [[Rcpp::export(".armaVARX1_Chat_ridgeSS")]]
arma::mat armaVARX1_Chat_ridgeSS_forR(arma::mat& COVZ, arma::mat& VARZ, const double lambdaA, const double lambdaB, const arma::mat& targetA, const arma::mat& targetB){
	return armaVARX1_Chat_ridgeSS(VARZ, COVZ, lambdaA, lambdaB, targetA, targetB);
}

// [[Rcpp::export(".armaVARX1_Chat_ridgeML")]]
arma::mat armaVARX1_Chat_ridgeML_forR(arma::mat& P, arma::mat& COVZ, arma::mat& eigvecVARZ, const arma::vec eigvalVARZ, const double lambdaA, const double lambdaB, const arma::mat& targetA, const arma::mat& targetB){
	return armaVARX1_Chat_ridgeMLapprox(P, COVZ, eigvecVARZ, eigvalVARZ, lambdaA, lambdaB, targetA, targetB);
}

// [[Rcpp::export(".armaVARX1_convergenceEvaluation")]]
double armaVARX1_convergenceEvaluation_forR(arma::mat& Chat, arma::mat& Cprev, arma::mat& Phat, arma::mat& Pprev){
	// see armaVARX1_convergenceEvaluation
	return armaVARX1_convergenceEvaluation(Chat, Cprev, Phat, Pprev);
}

// [[Rcpp::export(".armaVARX1_COVZhat")]]
arma::mat armaVARX1_COVZhat_forR(Rcpp::NumericVector& Zraw, const int pY){
	// Covariance estimation of the VARX(1) model from multivariate time series data.
	// First reformatting of data to an arma::cube format
	arma::cube Z = armaVAR_array2cube_withoutMissing(Zraw); 
	return armaVARX1_COVZhat(Z, pY);
}

// [[Rcpp::export(".armaVARX1_Shat_ML")]]
arma::mat armaVARX1_Shat_ML_forR(Rcpp::NumericVector& Zraw, arma::mat& A, arma::mat& B){
	// Least squares estimation of the regression parameter A of a VAR(1) model
	// First reformatting of data to an arma::cube format
	arma::cube Z = armaVAR_array2cube_withoutMissing(Zraw); 
	return armaVARX1_Shat_ML(Z, A, B);
}

// [[Rcpp::export(".armaVARX1_loglik_LOOCVinternal")]]
double armaVARX1_loglikLOOCVinternal_forR(arma::vec Yt1, arma::vec Yt0, arma::vec Xt0, arma::mat& A, arma::mat& B, arma::mat& P){

	double LL = 0;
	int dofs = 0;
	Yt1 = Yt1 - A * Yt0 - B * Xt0;
	double logDetP;
	double detPsign;
	arma::log_det(logDetP, detPsign, P);
		
	return - arma::trace(Yt1.t() * P * Yt1) / 2 + logDetP / 2;
}


