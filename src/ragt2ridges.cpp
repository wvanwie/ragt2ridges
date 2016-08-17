// Rcpp::Rcout << "wessel is gek 1" << std::endl;					

// We only include RcppArmadillo.h which pulls Rcpp.h in for us
// # include <Rcpp.h>
# include <RcppArmadillo.h>

// These are not needed:
// using namespace std;
// using namespace Rcpp;
// using namespace RcppArmadillo;
// [[Rcpp::depends("Rcpp")]]
// [[Rcpp::depends("RcppArmadillo")]]


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
FUNCTIONS FOR THE RIDGE ML ESTIMATION OF THE PRECISION MATRIX.
TAKEN DIRECTLY FROM THE RAGS2RIDGES-PACKAGE AND TO BE REMOVED WHEN THEY CAN BE DIRECTLY IMPORTED FROM THAT PACKAGE.
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline arma::mat rev_eig(const arma::vec eigval, const arma::mat eigvec) {
  /* ---------------------------------------------------------------------------
   "Reverse" the eigen decomposition, i.e. perform the multiplcation
   - eigval > a vector of eigenvalues
   - eigvec > a matrix of corresponding eigenvectors
  --------------------------------------------------------------------------- */

  return eigvec*diagmat(eigval)*eigvec.t();  // Could be more efficient
}

// [[Rcpp::export(.armaRidgePAnyTarget)]]
arma::mat armaRidgePAnyTarget(const arma::mat & S,
                              const arma::mat & target,
                              const double lambda,
                              int invert = 2) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for general/arbitrary targets.
   Depending on the value of "invert"" using matrix inversion (via
   diagonalization) or avoiding it.
   - S      > A sample covariance matrix. Should not contain NAs, Infs, or NaNs!
   - target > The target matrix with the same size as S
   - lambda > The the ridge penalty
   - invert > integer. Should the estimate be compute using inversion?
              0 = "no", 1 = "yes", 2 = "automatic" (default).
  --------------------------------------------------------------------------- */

  arma::vec eigvals;
  arma::mat eigvecs = S - lambda*target;
  if (!eigvecs.is_finite()) {
    return target;
  }
  eig_sym(eigvals, eigvecs, eigvecs, "dc");
  eigvals = 0.5*eigvals;
  arma::vec sqroot = sqrt(lambda + pow(eigvals, 2.0));

  // Return target if shrunken evals are infinite and lambda is "large"
  // Usually happens for lambda >= 1e154
  if (lambda > 1e6 && (!eigvals.is_finite() || !sqroot.is_finite())) {
    return target;
  }

  arma::vec D_inv = 1.0/(sqroot + eigvals); // inversion diagonal

  if (invert == 2) {   // Determine to invert or not
    if (lambda > 1) {  // Generally, don't use inversion for "large" lambda
      invert = 0;
    } else {
      if (!D_inv.is_finite()) {
        invert = 0;
      } else {
        invert = 1;
      }
    }
  }

  // Determine to invert or not
  if (invert == 1) {
    return rev_eig(D_inv, eigvecs);  // Proper inversion
  } else {
    arma::vec D_noinv = (sqroot - eigvals)/lambda; // inversion-less diagonal
    return rev_eig(D_noinv, eigvecs);  // Inversion by proposion
  }

}

// [[Rcpp::export(.armaRidgePScalarTarget)]]
arma::mat armaRidgePScalarTarget(const arma::mat & S,
                                 const double alpha,
                                 const double lambda,
                                 int invert = 2) {
  /* ---------------------------------------------------------------------------
   Compute the ridge estimate for rotational equivariant targets.
   Depending on the value of "invert"" using matrix inversion (via
   diagonalization) or avoiding it.
   - S      > A sample covariance matrix.
   - alpha  > The scaling of the identity matrix. Shoud not contain NaNs, Infs,
              or NA.s
   - lambda > The ridge penalty. Can be set to Inf (on the R side)
   - invert > Should the estimate be compute using inversion?
              0 = "no", 1 = "yes", 2 = "automatic", (default).
  --------------------------------------------------------------------------- */

  arma::vec eigvals;
  arma::mat eigvecs;
  arma::eig_sym(eigvals, eigvecs, S, "dc");

  eigvals = 0.5*(eigvals - lambda*alpha);
  arma::vec sqroot = sqrt(lambda + pow(eigvals, 2.0));

  // Return target if shrunken evals are infinite and lambda is "large"
  // Usually happens for lambda >= 1e154
  if (lambda > 1e6 && (!eigvals.is_finite() || !sqroot.is_finite())) {
    const int p = S.n_rows;
    return alpha*arma::eye<arma::mat>(p, p);
  }

  arma::vec D_inv = 1.0/(sqroot + eigvals); // inversion diagonal

  if (invert == 2) {   // Determine to invert or not
    if (lambda > 1) {  // Generally, don't use inversion for "large" lambda
      invert = 0;
    } else {
      if (!D_inv.is_finite()) {
        invert = 0;
      } else {
        invert = 1;
      }
    }
  }

  // Determine to invert or not
  if (invert == 1) {
    return rev_eig(D_inv, eigvecs);  // Proper inversion
  } else {
    arma::vec D_noinv = (sqroot - eigvals)/lambda; // inversion-less diagonal
    return rev_eig(D_noinv, eigvecs);  // Inversion by proposion
  }

}

// [[Rcpp::export(.armaRidgeP)]]
arma::mat armaRidgeP(const arma::mat & S,
                     const arma::mat & target,
                     const double lambda,
                     int invert = 2) {
  /* ---------------------------------------------------------------------------
   The ridge estimator in C++. Wrapper for the subroutines
   - S      > The sample covariance matrix (a numeric matrix on the R side)
   - target > Target matrix (a numeric matrix on the R side, same size as S)
   - lambda > The penalty (a numeric of length one on the R side)
   - invert > Should the estimate be compute using inversion?
              0 = "no", 1 = "yes", 2 = "automatic", (default).
  --------------------------------------------------------------------------- */

 if (lambda <= 0) {
    Rcpp::stop("The penalty (lambda) must be strictly postive");
  }

  if (lambda == arma::datum::inf) {
    return target;
  }

  const int p = S.n_rows;
  const double alpha = target(0, 0);
  arma::mat alphaI = arma::zeros<arma::mat>(p, p);
  alphaI.diag() += alpha;

  if (arma::all(arma::all(target == alphaI))) {
    return armaRidgePScalarTarget(S, alpha, lambda, invert);
  } else {
    return armaRidgePAnyTarget(S, target, lambda, invert);
  }

}


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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
SUPPORT FUNCTIONS FOR RIDGE PRECISION ESTIMATION WITH CHORDAL SUPPORT
INTENDED TO BE INCLUDED IN RAGS2RIDGES-PACKAGE WHEN DIRECT IMPORTATION FROM THAT SOURCE IS POSSIBLE
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export(".armaRidgePchordalInit")]]
arma::mat armaRidgePchordalInitWorkhorse(arma::mat S, const double lambda, arma::mat target, std::string type, Rcpp::List Cliques, Rcpp::List Separators){
	////////////////////////////////////////////////////////////////////////////
	// construct (guess of) ridge MLE of precision
	// various types (Alt, ArchI, ArchII) of the ridge estimor are implemented
	////////////////////////////////////////////////////////////////////////////

	// start with null matrix
	arma::mat P = arma::zeros<arma::mat>(S.n_rows, S.n_rows);

	// first archetypal ridge estimator
	if (type == "ArchI"){ 
		arma::mat targetInv = arma::inv_sympd(target); 
		// account for clique contributions
		for (int k = 0; k < Cliques.size(); ++k){
			arma::ivec slh = Cliques[k];
			slh = slh - 1;
			arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
			P.submat(cliqNodes, cliqNodes) =  P.submat(cliqNodes, cliqNodes) + 
				arma::inv_sympd((1-lambda) * S.submat(cliqNodes, cliqNodes) + lambda * targetInv.submat(cliqNodes, cliqNodes));
		}
		// account for separator contributions
	 	for (int k = 0; k < Separators.size(); ++k){
			arma::ivec slh = Separators[k];
			if (slh.n_elem > 0){
				slh = slh - 1;
				arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
        	 		P.submat(sepNodes, sepNodes) =  arma::symmatl(P.submat(sepNodes, sepNodes)) - 
					arma::inv_sympd((1-lambda) * S.submat(sepNodes, sepNodes) + lambda * targetInv.submat(sepNodes, sepNodes));
			}
		}
	}

	// second archetypal ridge estimator
	if (type == "ArchII"){ 
		// account for clique contributions
		for (int k = 0; k < Cliques.size(); ++k){
			arma::ivec slh = Cliques[k];
			slh = slh - 1;
			arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
			P.submat(cliqNodes, cliqNodes) = arma::symmatl(P.submat(cliqNodes, cliqNodes)) + 
				arma::inv_sympd(S.submat(cliqNodes, cliqNodes) + lambda * arma::eye(cliqNodes.n_elem, cliqNodes.n_elem));
		}
		// account for separator contributions
	 	for (int k = 0; k < Separators.size(); ++k){
			arma::ivec slh = Separators[k];
			if (slh.n_elem > 0){
				slh = slh - 1;
				arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
				P.submat(sepNodes, sepNodes) = arma::symmatl(P.submat(sepNodes, sepNodes)) - 
					arma::inv_sympd(S.submat(sepNodes, sepNodes) + lambda * arma::eye(sepNodes.n_elem, sepNodes.n_elem));
			}
		}
	}

	// alternative ridge estimator
	if (type == "Alt"){ 
		// account for clique contributions
		for (int k = 0; k < Cliques.size(); ++k){
			arma::ivec slh = Cliques[k];
			slh = slh - 1;
			arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
			if (cliqNodes.n_elem > 1){ P.submat(cliqNodes, cliqNodes) = P.submat(cliqNodes, cliqNodes) + 
				armaRidgeP(S.submat(cliqNodes, cliqNodes), target.submat(cliqNodes, cliqNodes), lambda); 
			}
			if (cliqNodes.n_elem == 1){  
				arma::mat slh2 = S.submat(cliqNodes, cliqNodes) - lambda * target.submat(cliqNodes, cliqNodes);
				double slh3 = slh2(0,0);
				P.submat(cliqNodes, cliqNodes) = P.submat(cliqNodes, cliqNodes) + 1/(sqrt(lambda + 0.25f * slh3 * slh3) + slh3 / 2);
			}
		}
		// account for separator contributions
	 	for (int k = 0; k < Separators.size(); ++k){
			arma::ivec slh = Separators[k];
			if (slh.n_elem > 0){
				slh = slh - 1;
				arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
				if (sepNodes.n_elem > 1){ P.submat(sepNodes, sepNodes) = P.submat(sepNodes, sepNodes) - 
					armaRidgeP(S.submat(sepNodes, sepNodes), target.submat(sepNodes, sepNodes), lambda); 
				}
				if (sepNodes.n_elem == 1){  
					arma::mat slh2 = S.submat(sepNodes, sepNodes) - lambda * target.submat(sepNodes, sepNodes);
					double slh3 = slh2(0,0);
				 	P.submat(sepNodes, sepNodes) = P.submat(sepNodes, sepNodes) - 1/(sqrt(lambda + 0.25f * slh3 * slh3) + slh3 / 2);
				}
			}
		}
	}
	
	// return initial ridge guess for chordal precision matrix
	return(P);
}

// [[Rcpp::export(".armaPenLLreparPforNLM")]]
Rcpp::NumericVector armaPenLLreparPforNLM(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// ridge penalized log-likelihood (and its gradient) for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005).
	// passed on to the 'nlm' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

	// return (minus) penalized log-likelihood
	Rcpp::NumericVector penLL(1);
	penLL = -log(det(P)) + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));

	// return gradient of (minus) penalized log-likelihood: nonzero elements only
	int p = S.n_rows;
	const arma::mat gradMat = 2 * (S - arma::inv_sympd(P) + lambda * (P - target));
	arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
	arma::vec grad = gradMat.elem(elemID);
	penLL.attr("gradient") = grad;

	return penLL;
}

// [[Rcpp::export(".armaPenLLreparP")]]
const double armaPenLLreparP(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// ridge penalized log-likelihood for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005)
	// passed on to the 'optim' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

	// return (minus) penalized log-likelihood
    // double logDetP; double detPsign;
    // arma::log_det(logDetP, detPsign, P);	
	// const double penLL = -logDetP + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));
	const double penLL = -log(det(P)) + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));	
	return penLL;
}

// [[Rcpp::export(".armaPenLLreparPgrad")]]
arma::vec armaPenLLreparPgrad(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// gradient of ridge penalized log-likelihood for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005)
	// passed on to the 'optim' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

	// return gradient of (minus) penalized log-likelihood: nonzero elements only
	int p = S.n_rows;
	const arma::mat gradMat = 2 * (arma::symmatl(S) - arma::symmatl(arma::inv_sympd(P)) + lambda * (arma::symmatl(P) - arma::symmatl(target)));
	arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
	return gradMat.elem(elemID);
}

// [[Rcpp::export(".armaPenLLreparGradArchI")]]
arma::vec armaPenLLreparGradArchI(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// gradient of 'ArchI'-ridge penalized log-likelihood for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005)
	// passed on to the 'optim' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);
 
	// return gradient of (minus) penalized log-likelihood: nonzero elements only
	int p = S.n_rows;
	const arma::mat gradMat = 2 * ((1-lambda) * arma::symmatl(S) + lambda * arma::symmatl(arma::inv_sympd(target)) - arma::symmatl(arma::inv_sympd(P)));
	arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
	return gradMat.elem(elemID);
}

// [[Rcpp::export(".armaPenLLreparGradArchII")]]
arma::vec armaPenLLreparGradArchII(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// gradient of 'ArchII-ridge penalized log-likelihood for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005)
	// passed on to the 'optim' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

	// return gradient of (minus) penalized log-likelihood: nonzero elements only
	int p = S.n_rows;
	const arma::mat gradMat = 2 * (arma::symmatl(S) + lambda * arma::eye(p,p) - arma::inv_sympd(P));
	arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
	return gradMat.elem(elemID);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
AUXILIARY FUCTIONS FOR THE ESTIMATION OF THE VAR-TYPE MODEL: FOR USE ON THE SEASIDE ONLY. WRAPPERS FOR EXPORTATION TO R PROVIDED AT END OF FILE, IF CONSIDERED OF USE.
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline arma::mat armaP_ridge_diag(const arma::vec & S, const arma::vec & target, const double lambda){
	/* ---------------------------------------------------------------------------
	The ridge estimator in C++. Wrapper for the subroutines
	- S      > The sample covariance matrix (a numeric matrix on the R side)
	- target > Target matrix (a numeric matrix on the R side, same size as S)
	- lambda > The penalty (a numeric of length one on the R side)
	- invert > Should the estimate be compute using inversion?
              0 = "no", 1 = "yes", 2 = "automatic", (default).
	--------------------------------------------------------------------------- */

	if (lambda <= 0) {
		Rcpp::stop("The penalty (lambda) must be strictly postive");
	}

	if (lambda == arma::datum::inf) {
		return target;
	}
	arma::vec U = (S - lambda * target) / 2;
	return arma::diagmat((sqrt(U % U + lambda) - U) / lambda);
}

inline Rcpp::List armaEigenDecomp_blockDiagOnly(const arma::mat symMat, const arma::ivec blockDims){
	/////////////////////////////////////////////////////////////////////////////////////////
	// eigendecomposition of symmetric matrix: block diagonals only
	/////////////////////////////////////////////////////////////////////////////////////////
	arma::vec eigvals = arma::zeros(arma::sum(blockDims));
	arma::mat eigvecs = arma::zeros(arma::sum(blockDims), arma::sum(blockDims));
	arma::vec eigvalsSLH;
	arma::mat eigvecsSLH;
	arma::uvec blockDimsU = arma::conv_to<arma::uvec>::from(blockDims);	
	unsigned int blockStart = 0;
	unsigned int blockEnd = blockDimsU[0]-1;
	for (unsigned int k = 0; k < blockDimsU.n_elem; ++k) {
		arma::eig_sym(eigvalsSLH, eigvecsSLH, symMat.submat(blockStart, blockStart, blockEnd, blockEnd), "dc");
		eigvecs.submat(blockStart, blockStart, blockEnd, blockEnd) = eigvecsSLH;
		eigvals.subvec(blockStart, blockEnd) = eigvalsSLH;
		blockStart = blockStart + blockDimsU[k];
		if (k+1 < blockDims.n_elem){
			blockEnd = blockEnd + blockDimsU[k+1];
		}
	}
	return Rcpp::List::create(Rcpp::Named("values") = eigvals, Rcpp::Named("vectors") = eigvecs);
}

inline arma::cube armaArray2cube(Rcpp::NumericVector Yraw){
	/////////////////////////////////////////////////////////////////////////////////////////
	// reformatting of data to an arma::cube format
	/////////////////////////////////////////////////////////////////////////////////////////

	Rcpp::NumericVector vecArray(Yraw);
	Rcpp::IntegerVector arrayDims = vecArray.attr("dim");
	arma::cube Y(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
	return Y;
}

inline arma::mat armaVARzerosCorrection_dense(arma::mat Ahat, int pRows, int pCols, const arma::mat eigvecP, const arma::vec eigvalP, const arma::mat eigvecVAR, const arma::vec eigvalVAR, const double lambda, const arma::uvec zerosR, const arma::uvec zerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// evaluates and applies the imposed zero constraint correction, assuming few zeros
	/////////////////////////////////////////////////////////////////////////////////////////

	// efficient calculation of the correction under the assumption of few zeros
	arma::uvec zerosID = (zerosC - 1) * pRows + zerosR - 1;
	arma::mat lambdaMat = arma::mat(pRows, pCols); 
	lambdaMat.fill(lambda);
	arma::vec eigvalShrunkSLH = arma::vectorise(eigvalP * arma::trans(eigvalVAR) + lambdaMat);
	arma::mat Vu = arma::mat(pRows * pCols, zerosR.n_elem);
	arma::mat eigvalShrunk = arma::mat(pRows * pCols, zerosR.n_elem);
	for (unsigned int z = 0; z < zerosR.n_elem; ++z) {
		Vu.col(z) = arma::vectorise(arma::trans(eigvecP.row(zerosR(z)-1)) * eigvecVAR.row(zerosC(z)-1));
		eigvalShrunk.col(z) = eigvalShrunkSLH;
	}
	arma::mat Q = arma::inv(arma::trans(Vu / eigvalShrunk) * Vu);
	arma::mat slh;
	for (unsigned int z = 0; z < zerosR.n_elem; ++z) {
	    slh = arma::reshape(Vu.col(z) / eigvalShrunk.col(z), pRows, pCols);
	    Vu.col(z) = arma::reshape(eigvecP * slh * arma::trans(eigvecVAR), pRows * pCols, 1);
	}
	return Ahat - reshape(Vu * Q * Ahat.elem(zerosID), pRows, pCols);
}

inline arma::mat armaVARzerosCorrection_sparse(arma::mat Ahat, int pRows, int pCols, const arma::mat P, const arma::mat eigvecVAR, const arma::vec eigvalVAR, const double lambda, const arma::uvec zerosR, const arma::uvec zerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// evaluates and applies the imposed zero constraint correction, assuming many zeros
	/////////////////////////////////////////////////////////////////////////////////////////

	// efficient calculate of the correction under the assumption of few nonzeros
	arma::mat nonzerosAmat = arma::mat(pRows, pCols, arma::fill::ones);
	arma::uvec zerosID = (zerosC - 1) * pRows + zerosR - 1;
	nonzerosAmat.elem(zerosID) = arma::zeros(zerosID.n_elem);
	arma::uvec nonzerosAvec = find(nonzerosAmat != 0);
	arma::uvec nonzerosC = arma::floor(nonzerosAvec / pRows);
	arma::uvec nonzerosR = nonzerosAvec - pRows * nonzerosC;    
	arma::mat VAR = eigvecVAR * arma::diagmat(eigvalVAR) * arma::trans(eigvecVAR);
	arma::vec lambdaNZ = arma::vec(nonzerosC.n_elem);
	lambdaNZ.fill(lambda);
	arma::mat slh1 = arma::inv(arma::diagmat(lambdaNZ) + VAR(nonzerosC, nonzerosC) % P(nonzerosR, nonzerosR));
	arma::vec slh2 = arma::vec(nonzerosAvec.n_elem, arma::fill::zeros);
	arma::mat slh3 = arma::mat(nonzerosAvec.n_elem, nonzerosAvec.n_elem);
	for (unsigned int z = 0; z < nonzerosAvec.n_elem; ++z) {        
		slh3 = arma::trans(P.row(nonzerosR(z))) * VAR.row(nonzerosC(z));
		slh2(z) += sum(slh3.elem(zerosID) % Ahat.elem(zerosID));
	}
	Ahat.elem(zerosID) = arma::zeros(zerosID.n_elem);
	Ahat.elem(nonzerosAvec) += slh1 * slh2;
	return Ahat;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
AUXILIARY FUCTIONS FOR THE ESTIMATION OF THE VAR(1) MODEL: FOR USE ON THE SEASIDE ONLY. WRAPPERS FOR EXPORTATION TO R PROVIDED AT END OF FILE.
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline arma::cube armaVAR_array2cube_withMissing(Rcpp::NumericVector Yraw, const arma::uvec unbalancedR, const arma::uvec unbalancedC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// set profiles of missing (time, sample)-points to missing
	/////////////////////////////////////////////////////////////////////////////////////////

	// reformatting of data to an arma::cube format
	const Rcpp::NumericVector vecArray(Yraw);
	Rcpp::NumericVector vecArraySlh = Rcpp::clone(vecArray);
	Rcpp::IntegerVector arrayDims = vecArraySlh.attr("dim");
  	arma::cube Y(vecArraySlh.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

	// insert missings per slice 
	int nMissing = unbalancedR.n_elem;
	for (int k = 0; k < nMissing; ++k){
		arma::mat Yk = Y.slice(unbalancedC(k)-1);
		arma::vec w = Yk.col(0);
		std::fill(w.begin(), w.end(), NA_REAL);
		Yk.col(unbalancedR(k)-1) = w;
		Y.slice(unbalancedC(k)-1) = Yk;
	}
	return Y;
}

inline arma::cube armaVAR_array2cube_withoutMissing(Rcpp::NumericVector Yraw){
	/////////////////////////////////////////////////////////////////////////////////////////
	// set profiles of missing (time, sample)-points to missing
	/////////////////////////////////////////////////////////////////////////////////////////

	// reformatting of data to an arma::cube format
	const Rcpp::NumericVector vecArray(Yraw);
	Rcpp::NumericVector vecArraySlh = Rcpp::clone(vecArray);
	Rcpp::IntegerVector arrayDims = vecArraySlh.attr("dim");
	arma::cube Y(vecArraySlh.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
	return Y;
}

inline arma::mat armaVAR1_Shat_ML(arma::cube Y, const arma::mat A){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Least squares estimation of the regression parameter A of a VAR(1) model
	/////////////////////////////////////////////////////////////////////////////////////////

	// covariance estimation 
	int p = Y.n_rows; int T = Y.n_cols; int n = Y.n_slices;

	// define dof and S
	arma::mat S = arma::symmatl(arma::zeros(p, p));
	int dof = 0; arma::uvec ids; arma::mat Yi;	
	for (int i = 0; i < n; ++i) {
		Yi = Y.slice(i);
		Yi = Yi.submat(0, 1, p-1, T-1) - A * Yi.submat(0, 0, p-1, T-2);
		ids = arma::find_finite(arma::sum(Yi));
		Yi = Yi.cols(ids);
		S = arma::symmatl(S) + arma::symmatl(Yi * arma::trans(Yi));
		dof = dof + ids.size();
	}
	return S / dof;
}

inline const double armaVAR1_loglik(arma::cube Y, arma::mat A, arma::mat P){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Evaluate the log-likelihood of the VAR1 model
	/////////////////////////////////////////////////////////////////////////////////////////

	// covariance estimation plus d.o.f. calculation 
	int p = Y.n_rows; int T = Y.n_cols; int n = Y.n_slices; arma::mat Yi;

	// declare variables
	double LL = 0; int dofs = 0;
	for (int i = 0; i < n; ++i){
		// extract data of sample i
		Yi = Y.slice(i);

		// obtain residuals
		Yi = Yi.cols(1, T-1) - A * Yi.cols(0, T-2);

		// limit residuals to those time points without missing
		arma::mat colsums = arma::sum(Yi);
		arma::uvec idT = find_finite(colsums);
		Yi = Yi.cols(idT);

		// update log-likelihood and d.o.f.
		LL = LL - arma::trace(Yi.t() * P * Yi) / 2;
		dofs = dofs + idT.n_elem;
	}
	
	// return likelihood of the VAR(1) model
	double logDetP; double detPsign;
	arma::log_det(logDetP, detPsign, P);		
	return - dofs * p * log(2 * arma::datum::pi) / 2 + dofs * logDetP /2 + LL;
}

inline arma::mat armaVAR1_COVYhat(arma::cube Y){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Covariance estimation of the multivariate time series data.
	/////////////////////////////////////////////////////////////////////////////////////////

	// covariance estimation plus d.o.f. calculation 
	int COVYdof = 0; int p = Y.n_rows;
	int T = Y.n_cols; int n = Y.n_slices;
	arma::vec missing = arma::ones(T);
	arma::mat COVY = arma::mat(p, p, arma::fill::zeros);
	arma::mat Yk; arma::uvec slh;
	for (int k = 0; k < n; ++k) {
		missing = arma::vec(T, arma::fill::ones);
		Yk = Y.slice(k);
		for (int l = 0; l < T; ++l){
			// ensure to exclude missings
			slh = arma::find_nonfinite(Yk.col(l));
			if (slh.n_elem != 0){ missing(l) = 0; } 
		}
        COVYdof = COVYdof + arma::sum(missing.subvec(0,T-2) % missing.subvec(1,T-1)); 		
		Yk.elem( arma::find_nonfinite(Yk) ).zeros();
		COVY = COVY + Yk.submat(0, 1, p-1, T-1) * arma::trans(Yk.submat(0, 0, p-1, T-2));
	}
	return COVY / COVYdof;
}

inline arma::mat armaVAR1_VARYhat(arma::cube Y, bool efficient, arma::mat unbalanced ){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Variance estimation of the multivariate time series data.
	/////////////////////////////////////////////////////////////////////////////////////////
 
	// variance estimation plus d.o.f. calculation
	int VARYdof = 0; int p = Y.n_rows; int T = Y.n_cols; int n = Y.n_slices; arma::uvec idKeep;	
	if (efficient){		
		idKeep = arma::find(unbalanced.col(0) != T);	
		if (idKeep.n_rows > 0){            	
			unbalanced = unbalanced.rows(idKeep);
		}
	}
	int nMissing = unbalanced.n_rows;    
	arma::mat VARY = arma::mat(p, p, arma::fill::zeros);
	arma::mat Yk;
	for (int k = 0; k < n; ++k) {
		Yk = Y.slice(k);
		Yk.elem( arma::find_nonfinite(Yk) ).zeros();
		if (!efficient){ Yk.shed_col(T-1); } 
		VARY = VARY + arma::symmatl(Yk * arma::trans(Yk));
		if (efficient){
		    VARYdof = n * T - nMissing;
		} else {
		    VARYdof = n * (T-1) - nMissing;			
		}
	}
	return VARY / VARYdof;
}

inline arma::mat armaVAR1_Ahat_ridgeML(arma::mat P, arma::mat COVY, const arma::mat eigvecVARY, const arma::vec eigvalVARY, const double lambdaA, arma::mat targetA){
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

inline arma::mat armaVAR1_Ahat_ridgeML_speed(arma::mat P, arma::mat COVY, const arma::mat eigvecVARY, const arma::vec eigvalVARY, const double lambdaA, arma::mat targetA){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Ridge ML estimate the regression coefficient matrix A of a VAR model
	/////////////////////////////////////////////////////////////////////////////////////////

	// eigendecomposition of error precision P
	arma::vec eigvalP; arma::mat eigvecP;
	arma::eig_sym(eigvalP, eigvecP, P, "dc");
	eigvalP = arma::real(eigvalP);

	// evaluate in return the estimate of A
	return eigvecP * ((arma::trans(eigvecP) * (targetA + P * COVY)) / (eigvalP * arma::trans(eigvalVARY) + lambdaA)) * arma::trans(eigvecVARY);
}

inline arma::mat armaVAR1_Ahat_ridgeSS(arma::mat VARY, const arma::mat & COVY, const double & lambdaA, arma::mat & targetA){
	/////////////////////////////////////////////////////////////////////////////////////////
	// Least squares estimation of the regression parameter A of a VAR(1) model
	/////////////////////////////////////////////////////////////////////////////////////////

	// evaluate the estimate of A
	VARY.diag() += lambdaA;
	return (targetA + COVY) * arma::inv_sympd(VARY);
}

inline const double armaVAR1_convergenceEvaluation(arma::mat Ahat, arma::mat Aprev, arma::mat Phat, arma::mat Pprev){
	/////////////////////////////////////////////////////////////////////////////////////////
	// assess maximum difference between the current and previous estimates
	/////////////////////////////////////////////////////////////////////////////////////////

	// extract element-wise maximum of both matrices
	arma::mat maxis = arma::max(abs(Ahat - Aprev), abs(Phat - Pprev));
	return maxis.max();
}

inline arma::mat armaVAR1_Ahat_zeros(const arma::mat & P, arma::mat & COVY, const arma::mat & eigvecVARY, const arma::vec & eigvalVARY, const double & lambdaA, const arma::mat & targetA, std::string fitA, const arma::uvec & zerosR, const arma::uvec & zerosC, std::string & zerosAfit){
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
    
	if (zerosAfit == "dense"){
		Ahat = armaVARzerosCorrection_dense(Ahat, pY, pY, eigvecP, eigvalP, eigvecVARY, eigvalVARY, lambdaA, zerosR, zerosC);
	} 
	if (zerosAfit == "sparse"){
		Ahat = armaVARzerosCorrection_sparse(Ahat, pY, pY, P, eigvecVARY, eigvalVARY, lambdaA, zerosR, zerosC);
	}

	return Ahat;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
WRAPPERS FUNCTIONs FOR THE VAR(1) MODEL ESTIMATION
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export(".armaVAR1_ridgeML")]]
Rcpp::List armaVAR1_ridgeML(Rcpp::NumericVector & Yraw, const double & lambdaA, const double lambdaP, arma::mat & targetA, arma::mat & targetP, std::string targetPtype, 
std::string fitA, arma::mat & unbalanced, bool & diagP, bool & efficient, const int & nInit, const double & minSuccDiff){

	// set profiles of missing (time, sample)-points to missing
	arma::cube Y;
	if (unbalanced.n_rows == 0){ Y = armaVAR_array2cube_withoutMissing(Yraw); } 
	if (unbalanced.n_rows > 0){ Y = armaVAR_array2cube_withMissing(Yraw, arma::conv_to<arma::uvec >::from(unbalanced.col(0)), arma::conv_to<arma::uvec >::from(unbalanced.col(1))); }

	// estimate A by SS minimization
	arma::mat VARY = armaVAR1_VARYhat(Y, efficient, unbalanced);
	arma::mat COVY = armaVAR1_COVYhat(Y);
	arma::mat Ahat = armaVAR1_Ahat_ridgeSS(VARY, COVY, lambdaA, targetA);

	// calculate Se
	arma::mat Se = armaVAR1_Shat_ML(Y, Ahat);

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

	if (fitA == "ml"){ 
		///////////////////////////////////////////////////////////////////////////////
		// estimate parameters by ML, using the SS estimates are initials
		///////////////////////////////////////////////////////////////////////////////

		// eigen-decomposition of VARY
		arma::vec eigvalsVARY; arma::mat eigvecsVARY;
		arma::eig_sym(eigvalsVARY, eigvecsVARY, VARY, "dc");

		// declare storage variables
		arma::mat Aprev; arma::mat Pprev;
        
        // for speed a multiplication is taken outside the for-loop 
		COVY = COVY * eigvecsVARY;
		targetA = targetA * eigvecsVARY;
    
		for (int i = 0; i < nInit; ++i){
			// store latest estimates
			Aprev = Ahat; Pprev = Phat;

			// estimate A
			Ahat = armaVAR1_Ahat_ridgeML_speed(Phat, COVY, eigvecsVARY, eigvalsVARY, lambdaA, targetA);
    
			// calculate Se
			Se = armaVAR1_Shat_ML(Y, Ahat);

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
			if (armaVAR1_convergenceEvaluation(Ahat, Aprev, Phat, Pprev) < minSuccDiff){ break; }
		}
	}
    
	// evaluate likelihood
    double logDetP; double detPsign;
    arma::log_det(logDetP, detPsign, Phat);
	double LL = (Y.n_cols - 1) * Y.n_slices * (logDetP - arma::accu(Se % Phat)) / 2;

	return Rcpp::List::create(Rcpp::Named("A") = Ahat, Rcpp::Named("P") = Phat, Rcpp::Named("LL")=LL);
}

// [[Rcpp::export(".armaVAR1_ridgeML_zerosA")]]
Rcpp::List armaVAR1_ridgeML_zerosA(Rcpp::NumericVector Yraw, const double lambdaA, const double lambdaP, arma::mat targetA, arma::mat targetP, std::string targetPtype, 
std::string fitA, arma::mat unbalanced, bool diagP, bool efficient, const int nInit, const double minSuccDiff, const arma::uvec zerosR, const arma::uvec zerosC, std::string zerosAfit){

	// set profiles of missing (time, sample)-points to missing
	arma::cube Y;
	if (unbalanced.n_rows == 0){ Y = armaVAR_array2cube_withoutMissing(Yraw); }
	if (unbalanced.n_rows > 0){ Y = armaVAR_array2cube_withMissing(Yraw, arma::conv_to<arma::uvec >::from(unbalanced.col(0)), arma::conv_to<arma::uvec >::from(unbalanced.col(1))); }

	// for computational efficiency redefine targetA
	// targetA = lambdaA * targetA;

	// estimate A by SS minimization
	arma::mat VARY = armaVAR1_VARYhat(Y, efficient, unbalanced);
	arma::mat COVY = armaVAR1_COVYhat(Y);
	
   	// eigen-decomposition of VARY
	arma::vec eigvalsVARY; arma::mat eigvecsVARY;
	arma::eig_sym(eigvalsVARY, eigvecsVARY, VARY, "dc");
	
	// estimate A by SS minimization	
	arma::mat Ahat = armaVAR1_Ahat_zeros(arma::eye<arma::mat>(targetA.n_rows,targetA.n_rows), COVY, eigvecsVARY, eigvalsVARY, lambdaA, targetA, fitA, zerosR, zerosC, zerosAfit);

	// calculate Se
	arma::mat Se = armaVAR1_Shat_ML(Y, Ahat);

	// ridge ML estimation of Se
	arma::mat Phat;
	if (targetPtype != "none"){ targetP = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); }
	if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
	if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }

	if (fitA == "ml"){
		///////////////////////////////////////////////////////////////////////////////
		// estimate parameters by ML, using the SS estimates as initials
		///////////////////////////////////////////////////////////////////////////////
    
		// declare storage variables
		arma::mat Aprev; arma::mat Pprev;
    
		for (int i = 0; i < nInit; ++i){
			// store latest estimates
			Aprev = Ahat; Pprev = Phat;

			// estimate A
			Ahat = armaVAR1_Ahat_zeros(Phat, COVY, eigvecsVARY, eigvalsVARY, lambdaA, targetA, fitA, zerosR, zerosC, zerosAfit);
    
			// calculate Se
			Se = armaVAR1_Shat_ML(Y, Ahat);

			// construct target precision (if not provided)
			if (targetPtype != "none"){ targetP = armaP_defaultTarget(Se, targetPtype, 0.0001, 0); }

			// ridge ML estimation of precision
			if (diagP){ Phat = armaP_ridge_diag(arma::diagvec(Se), arma::diagvec(targetP), lambdaP); }
			if (!diagP){ Phat = armaRidgeP(Se, targetP, lambdaP); }

			// assess convergence
			if (armaVAR1_convergenceEvaluation(Ahat, Aprev, Phat, Pprev) < minSuccDiff){ break; }
		}
	}
    
	// evaluate likelihood
	double val;
	double sign;
	log_det(val, sign, Phat);
	double LL = (Y.n_cols - 1) * Y.n_slices * (val - arma::accu(Se % Phat)) / 2;

	return Rcpp::List::create(Rcpp::Named("A") = Ahat, Rcpp::Named("P") = Phat, Rcpp::Named("LL")=LL);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
WRAPPERS FOR R
-> FUNCTIONS DEFINED FOR USE ON THE SEASIDE EMPLOYING 'INLINE' FOR SPEED
-> CODE BELOW MAKES THEM AVAILABLE IN R
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export(".armaEigenDecomp_blockDiagOnly")]]
Rcpp::List armaEigenDecomp_blockDiagOnly_forR(const arma::mat symMat, const arma::ivec blockDims){
	// see armaEigenDecomp_blockDiagOnly
	return armaEigenDecomp_blockDiagOnly(symMat, blockDims);
}

// [[Rcpp::export(".armaVAR1_Shat_ML")]]
arma::mat armaVAR1_Shat_ML_forR(Rcpp::NumericVector Yraw, const arma::mat A){
	// Least squares estimation of the regression parameter A of a VAR(1) model
	// First reformatting of data to an arma::cube format
	arma::cube Y = armaArray2cube(Yraw); 
	return armaVAR1_Shat_ML(Y, A);
}

// [[Rcpp::export(".armaVAR1_COVYhat")]]
arma::mat armaVAR1_COVYhat_forR(Rcpp::NumericVector Yraw){
	// Covariance estimation of the multivariate time series data.
	// First reformatting of data to an arma::cube format
	arma::cube Y = armaArray2cube(Yraw); 
	return armaVAR1_COVYhat(Y);
}

// [[Rcpp::export(".armaVAR1_VARYhat")]]
arma::mat armaVAR1_VARYhat_forR(Rcpp::NumericVector Yraw, bool efficient, arma::mat unbalanced){
	// Variance estimation of the multivariate time series data.
	// First reformatting of data to an arma::cube format
	arma::cube Y = armaArray2cube(Yraw); 
	return armaVAR1_VARYhat(Y, efficient, unbalanced);
}

// [[Rcpp::export(".armaVAR1_Ahat_zeros")]]
arma::mat armaVAR1_Ahat_zeros_forR(const arma::mat & P, arma::mat & COVY, const arma::mat & eigvecVARY, const arma::vec & eigvalVARY, const double & lambdaA, const arma::mat & targetA, std::string fitA, const arma::uvec & zerosR, const arma::uvec & zerosC, std::string & zerosAfit){
	// see armaVAR1_Ahat_zeros
	return armaVAR1_Ahat_zeros(P, COVY, eigvecVARY, eigvalVARY, lambdaA, targetA, fitA, zerosR, zerosC, zerosAfit);
}

// [[Rcpp::export(".armaP_defaultTarget")]]
arma::mat armaP_defaultTarget_forR(arma::mat S, std::string targetType, const double fraction, double const multiplier){
	// Rcpp version of 'default.target': key difference is the input checks.
	return armaP_defaultTarget(S, targetType, fraction, multiplier);
}

// [[Rcpp::export(".armaVAR_array2cube_withMissing")]]
arma::cube armaVAR_array2cube_withMissing_forR(Rcpp::NumericVector & Yraw, const arma::uvec & unbalancedR, const arma::uvec & unbalancedC){
	// set profiles of missing (time, sample)-points to missing
	// First reformatting of data to an arma::cube format
	return armaVAR_array2cube_withMissing(Yraw, unbalancedR, unbalancedC);
}

// [[Rcpp::export(".armaVAR_array2cube_withoutMissing")]]
arma::cube armaVAR_array2cube_withoutMissing_forR(Rcpp::NumericVector & Yraw){
	// reformatting of data to an arma::cube format
	return armaVAR_array2cube_withoutMissing_forR(Yraw);
}

// [[Rcpp::export(".armaEigenDecomp")]]
Rcpp::List armaEigenDecomp_forR(const arma::mat symMat){
	// eigendecomposition of a symmetric matrix
	arma::vec eigvals; arma::mat eigvecs;
	arma::eig_sym(eigvals, eigvecs, symMat, "dc");
	return Rcpp::List::create(Rcpp::Named("values") = eigvals, Rcpp::Named("vectors") = eigvecs);
}

// [[Rcpp::export(".armaVAR1_convergenceEvaluation")]]
const double armaVAR1_convergenceEvaluation_forR(arma::mat Ahat, arma::mat Aprev, arma::mat Phat, arma::mat Pprev){
	// assess maximum difference between the current and previous estimates
	return armaVAR1_convergenceEvaluation(Ahat, Aprev, Phat, Pprev);
}

// [[Rcpp::export(".armaVAR1_Ahat_ridgeML")]]
arma::mat armaVAR1_Ahat_ridgeML_forR(arma::mat P, arma::mat COVY, const arma::mat eigvecVARY, const arma::colvec eigvalVARY, const double lambdaA, arma::mat targetA){
	return armaVAR1_Ahat_ridgeML(P, COVY, eigvecVARY, eigvalVARY, lambdaA, targetA);
}

// [[Rcpp::export(".armaVAR1_Ahat_ridgeSS")]]
arma::mat armaVAR1_Ahat_ridgeSS_forR(arma::mat COVY, arma::mat VARY, const double & lambdaA, arma::mat & targetA){
	// see armaVARX1_Ahat_ridgeSS
	return armaVAR1_Ahat_ridgeSS(VARY, COVY, lambdaA, targetA);
}

// [[Rcpp::export(".armaVAR1_loglik_LOOCVinternal")]]
const double armaVAR1_loglikLOOCVinternal_forR(arma::vec Yt1, arma::vec Yt0, arma::mat A, arma::mat P){

	double logDetP; double detPsign;
	arma::log_det(logDetP, detPsign, P);
	
	Yt1 = Yt1 - A * Yt0;
	return - arma::trace(Yt1.t() * P * Yt1) / 2 + logDetP / 2;
}

// [[Rcpp::export(".armaVAR1_loglik")]]
const double armaVAR1_loglik_forR(Rcpp::NumericVector Yraw, arma::mat A, arma::mat P){
	// Evaluate the log-likelihood of the VAR1 model
	// First reformatting of data to an arma::cube format
	arma::cube Y = armaArray2cube(Yraw);
	return armaVAR1_loglik(Y, A, P);
}
