\name{mutualInfoVAR2}
\alias{mutualInfoVAR2}
\title{
Mutual information analysis of the VAR(2) model
}
\description{
Evaluate, within the VAR(2) model, the mutual information between a variate at one time point and the variates at a future time point.
}
\usage{
mutualInfoVAR2(A1, A2, SigmaE, T)
}
\arguments{
\item{A1}{ A \code{matrix} \eqn{\mathbf{A}_1} of lag one autoregression parameters. }
\item{A2}{ A \code{matrix} \eqn{\mathbf{A}_2} of lag two autoregression parameters. }
\item{SigmaE}{ Covariance \code{matrix} of the errors (innovations). }
\item{T}{ Positive \code{integer} of length one specifying the lag between time points for which the mutual informations are to be evaluated. }
}
\value{
Object of class \code{numeric} with elements corresponding to the mutual informations. 
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link{mutualInfoVAR1}}, \code{\link{ridgeVAR2}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A1     <- -createA(p, "clique", nCliques=1, nonzeroA=0.1)
A2     <- t(createA(p, "chain", nBands=1, nonzeroA=0.1))

# generate data
Y <- dataVAR2(n, T, A1, A2, SigmaE)

# fit VAR(1) model
VAR2hat <- ridgeVAR2(Y, 1, 1, 1)

# impulse response analysis
mutualInfoVAR2(VAR2hat$A1, VAR2hat$A2, solve(symm(VAR2hat$P)), 10)
}

