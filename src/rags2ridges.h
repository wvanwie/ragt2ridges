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
  arma::vec sqroot = arma::sqrt(lambda + arma::pow(eigvals, 2.0));

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
  arma::vec sqroot = arma::sqrt(lambda + arma::pow(eigvals, 2.0));

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

