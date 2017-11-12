\name{plotVAR1data}
\alias{plotVAR1data}
\title{
Time series plot
}
\description{
Plot of time series data. Per variate and individual a line connecting the observations at each time point is plotted.
}
\usage{
plotVAR1data(Y, lwd=1)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be zero centered covariate-wise, per sample. }
\item{lwd}{ A \code{numeric} of length one, specifying the line width. }
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
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
plotVAR1data(Y, lwd=2)
}

