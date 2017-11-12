\name{optPenaltyVAR1}
\alias{optPenaltyVAR1}
\title{
Automatic penalty parameter selection for the VAR(1) model.
}
\description{
Automatic penalty parameter selection for the VAR(1) model through maximization of the leave-one-out cross-validated (LOOCV) log-likelihood. 
}
\usage{
optPenaltyVAR1(Y, lambdaMin, lambdaMax, 
               lambdaInit=(lambdaMin+lambdaMax)/2, 
               optimizer="nlm", ...)
}
\arguments{
\item{Y}{ Three-dimensional array containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to centered covariate-wise. }
\item{ lambdaMin }{ A \code{numeric} of length two, containing the minimum values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}}, the matrix with autoregression coefficients, while the second parameter relates to the penalty on (\eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{lambdaMax}{ A \code{numeric} of length two, containing the maximum values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}}, the matrix with autoregression coefficients, while the second parameter relates to the penalty on (\eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{ lambdaInit }{ A \code{numeric} of length two, containing the initial values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}}, the matrix with autoregression coefficients, while the second parameter relates to the penalty on (\eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{optimizer}{ A \code{character} (either \code{nlm} (default) or \code{optim}) specifying which optimization function should be used: \code{\link[stats:nlminb]{nlminb}} (default) or \code{\link[stats:constrOptim]{constrOptim}}? }
\item{...}{ Additional arguments passed on to loglikLOOCVVAR1. }
}
\value{
A \code{numeric} with the LOOCV optimal choice for the ridge penalty parameter.
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link{loglikLOOCVVAR1}}, \code{\link{ridgeVAR1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A      <- createA(3, "chain")

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

# determine the optimal penalty parameter
optLambda <- optPenaltyVAR1(Y, rep(10^(-10), 2), rep(1000, 2))

# fit VAR(1) model
ridgeVAR1(Y, optLambda[1], optLambda[2])$A
}



