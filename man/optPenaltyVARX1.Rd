\name{optPenaltyVARX1}
\alias{optPenaltyVARX1}
\title{
Automatic penalty parameter selection for the VARX(1) model.
}
\description{
Automatic penalty parameter selection for the VARX(1) model through maximization of the leave-one-out cross-validated (LOOCV) log-likelihood. 
}
\usage{
optPenaltyVARX1(Y, X, lambdaMin, lambdaMax, 
                lambdaInit=(lambdaMin+lambdaMax)/2, 
                optimizer="nlm", ...)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the response data. The first, second and third dimensions correspond to variates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{X}{ Three-dimensional \code{array} containing the time-varying covariate data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{lambdaMin }{ A \code{numeric} of length three, containing the minimum values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}}, the matrix with autoregression coefficients, the second to matrix \eqn{\mathbf{B}} containing the regression coefficients of the time-varying covariates, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{lambdaMax}{ A \code{numeric} of length three, containing the maximum values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}}, the matrix with autoregression coefficients, the second to matrix \eqn{\mathbf{B}} containing the regression coefficients of the time-varying covariates, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{lambdaInit }{ A \code{numeric} of length three, containing the initial values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}}, the matrix with autoregression coefficients, the second to matrix \eqn{\mathbf{B}} containing the regression coefficients of the time-varying covariates, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{optimizer}{ A \code{character} (either \code{nlm} (default) or \code{optim}) specifying which optimization function should be used: \code{\link[stats:nlminb]{nlminb}} (default) or \code{\link[stats:constrOptim]{constrOptim}}? }
\item{...}{ Additional arguments passed on to the \code{\link{loglikLOOCVVARX1}}-function. }
}
\value{
A \code{numeric} with the LOOCV optimal choice for the ridge penalty parameter.
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link{loglikLOOCVVARX1}}, \code{\link{ridgeVARX1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
Ax     <- createA(3, "chain")

# generate time-varying covariate data
X <- dataVAR1(n, T, Ax, SigmaE)

# (auto)regression parameter matrices of VARX(1) model
A <- createA(p, topology="clique", nonzeroA=0.1, nClique=1)
B <- createA(p, topology="hub", nonzeroA=0.1, nHubs=1)

# generate data
Y <- dataVARX1(X, A, B, SigmaE, lagX=0)

# determine the optimal penalty parameter
optLambda <- optPenaltyVARX1(Y, X, rep(10^(-10), 3), rep(1000, 3), 
                             optimizer="nlm", lagX=0)

# fit VAR(1) model
ridgeVARX1(Y, X, optLambda[1], optLambda[2], optLambda[3], lagX=0)$A
}
