\name{optPenaltyVAR1fused}
\alias{optPenaltyVAR1fused}
\title{
Automatic penalty parameter selection for multiple VAR(1) models.
}
\description{
Automatic penalty parameter selection for multiple VAR(1) models through maximization of the leave-one-out cross-validated (LOOCV) log-likelihood. 
}
\usage{
optPenaltyVAR1fused(Y, id, lambdaMin, lambdaMax, 
                    lambdaInit=(lambdaMin+lambdaMax)/2, 
                    optimizer="nlm", ...)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the response data. The first, second and third dimensions correspond to variates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{id}{ A vector with group indices comprising of integers only. First group is represented by '0', the next by '1', and so on until the last. }
\item{lambdaMin }{ A \code{numeric} of length three, containing the minimum values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}_g}'s, the matrices with lag one autoregression coefficients, the second to the fused ridge parameter for these \eqn{\mathbf{A}_g}'s, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{lambdaMax}{ A \code{numeric} of length three, containing the maximum values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}_g}'s, the matrices with lag one autoregression coefficients, the second to the fused ridge parameter for these \eqn{\mathbf{A}_g}'s, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{lambdaInit }{ A \code{numeric} of length three, containing the initial values of ridge penalty parameters to be considered. The first element is the ridge parameter corresponding to the penalty on \eqn{\mathbf{A}_g}'s, the matrices with lag one autoregression coefficients, the second to the fused ridge parameter for these \eqn{\mathbf{A}_g}'s, while the third parameter relates to the penalty on \eqn{\mathbf{\Omega}_{\varepsilon}}, the precision matrix of the errors. }
\item{optimizer}{ A \code{character} (either \code{nlm} (default) or \code{optim}) specifying which optimization function should be used: \code{\link[stats:nlminb]{nlminb}} (default) or \code{\link[stats:constrOptim]{constrOptim}}? }

\item{...}{ Additional arguments passed on to loglikLOOCVVAR1fused. }
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
\code{\link{loglikLOOCVVAR1fused}}, \code{\link{ridgeVAR1fused}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points, G=groups)
p <- 3; n <- 12; T <- 10; G <- 3

# set model parameters
SigmaE       <- matrix(1/2, p, p)
diag(SigmaE) <- 1
A1           <- -createA(p, "clique", nCliques=1, nonzeroA=0.1)
A2           <- t(createA(p, "chain", nBands=1, nonzeroA=0.1))
A3           <- (A1 + A2) / 2

# generate data
Y1 <- dataVAR1(n/G, T, A1, SigmaE)
Y2 <- dataVAR1(n/G, T, A2, SigmaE)
Y3 <- dataVAR1(n/G, T, A3, SigmaE)
Y  <- abind::abind(Y1, Y2, Y3, along=3)
id <- c(rep(1, n/G), rep(2, n/G), rep(3, n/G))-1

# determine the optimal penalty parameter
\dontrun{optLambdas <- optPenaltyVAR1fused(Y, rep(10^(-10), 3), }
\dontrun{               rep(1000, 3), optimizer="nlm") }

# ridge ML estimation of the VAR(1) parameter estimates with 
# optimal penalty parameters
optLambdas <- c(0.1, 0.1, 0.1)
VAR1hats <- ridgeVAR1fused(Y, id, optLambdas[1], optLambdas[2], optLambdas[3])
}
