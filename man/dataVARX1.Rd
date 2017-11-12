\name{dataVARX1}
\alias{dataVARX1}
\title{
Sample data from a VARX(1) model
}
\description{
Sample data from a VARX(1) model.
}
\usage{
dataVARX1(X, A, B, SigmaE, lagX)
}
\arguments{
\item{X}{ Three-dimensional \code{array} containing the time-varying covariates. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{A}{ A \code{matrix} \eqn{\mathbf{A}} of autoregression parameters. }
\item{B}{ A \code{matrix} \eqn{\mathbf{B}} of regression parameters of the time-varying covariates stored in the array \code{X}. }
\item{SigmaE}{ Covariance \code{matrix} of the errors (innovations). }
\item{lagX}{ An \code{integer}, either \code{0} or {1}, specifying whether \eqn{\mathbf{X}_t} or \eqn{\mathbf{X}_{t-1}} affects \eqn{\mathbf{Y}_t}, respectively. }
}
\value{
A three dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link{createA}}, \code{\link[rags2ridges:createS]{createS}}, \code{\link{dataVAR1}}, \code{\link{dataVAR2}}, \code{\link{ridgeVARX1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
Ax     <- createA(3, "chain")

# generate time-varying covariate data
X <- dataVAR1(n, T, Ax, SigmaE)

# autoregression parameter matrices of VARX(1) model
A <- createA(p, topology="clique", nonzeroA=0.1, nClique=1)
B <- createA(p, topology="hub",    nonzeroA=0.1, nHubs=1)

# generate data
Y <- dataVARX1(X, A, B, SigmaE, lagX=0)
}

