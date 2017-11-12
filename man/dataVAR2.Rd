\name{dataVAR2}
\alias{dataVAR2}
\title{
Sample data from a VAR(2) model
}
\description{
Sample data from a VAR(2) model.
}
\usage{
dataVAR2(n, T, A1, A2, SigmaE, TburnIn=1000)
}
\arguments{
\item{n}{ Positive \code{numeric} of length one: number of individuals to be sampled. }
\item{T}{ Positive \code{numeric} of length one: number of time points (per individual) to be sampled. }
\item{A1}{ A \code{matrix} \eqn{\mathbf{A}_1} of lag one autoregression parameters. }
\item{A2}{ A \code{matrix} \eqn{\mathbf{A}_2} of lag two autoregression parameters. }
\item{SigmaE}{ Covariance \code{matrix} of the errors (innovations). }
\item{TburnIn}{ Positive \code{numeric} of length one: number of time points used to burn in the process. }
}
\value{
A three dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link{dataVAR1}}, \code{\link{dataVAR2}}, \code{\link{ridgeVARX1}}, .
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
}

