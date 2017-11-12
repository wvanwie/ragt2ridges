\name{loglikLOOCVVARX1}
\alias{loglikLOOCVVARX1}
\title{
Leave-one-out (minus) cross-validated log-likelihood of VARX(1) model
}
\description{
Evaluation of the (minus) leave-one-out cross-validated log-likelihood of the VARX(1) model for given choices of the ridge penalty parameters (\eqn{\lambda_a}, \eqn{\lambda_b} and \eqn{\lambda_{\omega}} for the autoregression coefficient matrix \eqn{\mathbf{A}}, regression coefficient matrix \eqn{\mathbf{B}} of time-varying covariates, and the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}, respectively). The functions also works with a (possibly) unbalanced experimental set-up. The VARX(1)-process is assumed to have mean zero.
}
\usage{
loglikLOOCVVARX1(lambdas, Y, X, unbalanced=matrix(nrow=0, ncol=2), lagX=0, ...)
}
\arguments{
\item{lambdas}{ A \code{numeric} of length three, comprising positive values only. It contains the ridge penalty parameters to be used in the estimation of \eqn{\mathbf{A}}, \eqn{\mathbf{B}} and the precision matrix of the errors, respectively. }
\item{Y}{ Three-dimensional \code{array} containing the response data. The first, second and third dimensions correspond to variates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{X}{ Three-dimensional \code{array} containing the time-varying covariate data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{lagX}{ \code{Integer}, either \code{0} or {1}, specifying whether \eqn{\mathbf{X}_t} or \eqn{\mathbf{X}_{t-1}} affects \eqn{\mathbf{Y}_t}, respectively. }
\item{...}{  Other arguments to be passed to \code{ridgeVARX1}. }
}
\value{
A \code{numeric} of length one: the minus (!) LOOCV log-likelihood.
}
\note{
The minus LOOCV log-likelihood is returned as standard optimization procedures in R like \code{nlminb} and \code{constrOptim} minimize (rather then maximize). Hence, by providing the minus LOOCV log-likelihood the function \code{\link{loglikLOOCVVARX1}} can directly used by these optimization procedures.
}

\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link[rags2ridges:ridgeP]{ridgeP}} and \code{\link{ridgeVARX1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
Ax     <- createA(3, "chain")

# generate time-varying covariate data
X <- dataVAR1(n, T, Ax, SigmaE)

# regression parameter matrices of VARX(1) model
A <- createA(p, topology="clique", nonzeroA=0.1, nClique=1)
B <- createA(p, topology="hub", nonzeroA=0.1, nHubs=1)

# generate data
Y <- dataVARX1(X, A, B, SigmaE, lagX=0)

## determine optimal values of the penalty parameters
\dontrun{optLambdas <- constrOptim(c(1,1, 1), loglikLOOCVVARX1, gr=NULL, }
\dontrun{              ui=diag(3), ci=c(0,0,0), Y=Y, X=X, lagX=0,}
\dontrun{              control=list(reltol=0.01))$par }

# ridge ML estimation of the VAR(1) parameter estimates with 
# optimal penalty parameters
optLambdas <- c(0.1, 0.1, 0.1)
ridgeVARX1(Y, X, optLambdas[1], optLambdas[2], optLambdas[3], lagX=0)$A
}

