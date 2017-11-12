\name{centerVAR1data}
\alias{centerVAR1data}
\title{
Zero-centering of time-course data
}
\description{
Per individual, covariate-wise zero centering of the time-series data.
}
\usage{
centerVAR1data(Y)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. }
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\value{
An \code{array} with dimensions as the input.
}
\seealso{
\code{\link{dataVAR1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A      <- createA(p, "chain")

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

# plot data sampled from the VAR(1) model.
Ycentered <- centerVAR1data(Y)
}

