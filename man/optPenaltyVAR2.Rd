\name{optPenaltyVAR2}
\alias{optPenaltyVAR2}
\title{
Automatic penalty parameter selection for the VAR(2) model.
}
\description{
Automatic penalty parameter selection for the VAR(2) model through maximization of the leave-one-out cross-validated (LOOCV) log-likelihood. 
}
\usage{
optPenaltyVAR2(Y, lambdaMin, lambdaMax, 
               lambdaInit=(lambdaMin+lambdaMax)/2, 
               optimizer="nlm", ...)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the response data. The first, second and third dimensions correspond to variates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{lambdaMin }{ A \code{numeric} of length three, containing the minimum values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}_1}, the matrix with lag one autoregression coefficients, the second to matrix \eqn{\mathbf{A}_2} containing the lag two autoregression coefficients, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{lambdaMax}{ A \code{numeric} of length three, containing the maximum values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}_1}, the matrix with lag one autoregression coefficients, the second to matrix \eqn{\mathbf{A}_2} containing the lag two autoregression coefficients, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{lambdaInit }{ A \code{numeric} of length three, containing the initial values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}_1}, the matrix with lag one autoregression coefficients, the second to matrix \eqn{\mathbf{A}_2} containing the lag two autoregression coefficients, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{optimizer}{ A \code{character} (either \code{nlm} (default) or \code{optim}) specifying which optimization function should be used: \code{\link[stats:nlminb]{nlminb}} (default) or \code{\link[stats:constrOptim]{constrOptim}}? }

\item{...}{ Additional arguments passed on to loglikLOOCVVAR2. }
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
\code{\link{loglikLOOCVVAR2}}, \code{\link{ridgeVAR2}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A1 <- createA(p, topology="clique", nonzeroA=0.1, nClique=1)
A2 <- createA(p, topology="hub", nonzeroA=0.1, nHubs=1)

# generate data
Y <- dataVAR2(n, T, A1, A2, SigmaE)

# determine the optimal penalty parameter
optLambda <- optPenaltyVAR2(Y, rep(10^(-10), 3), rep(1000, 3), 
                            optimizer="nlm")

# fit VAR(2) model
ridgeVAR2(Y, optLambda[1], optLambda[2], optLambda[3])
}
