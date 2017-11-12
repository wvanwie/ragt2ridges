\name{longitudinal2array}
\alias{longitudinal2array}
\title{
Convert a longitudinal object into an array.
}
\description{
Converts an object of the \code{longitudinal}-class into a 3-dim array (containing time-series data of multiple individuals).
}
\usage{
longitudinal2array(Y)
}
\arguments{
\item{Y}{ Object of the \code{\link[longitudinal:longitudinal]{longitudinal}} class. Essentially, this is a matrix with the timepoint slices of the array stacked on top of each other. }
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\value{
A \code{array}-object containing time-series data.
}
\seealso{
\code{\link{array2longitudinal}}, \code{\link{dataVAR1}}, \code{\link[longitudinal:longitudinal]{longitudinal}}, \code{\link[longitudinal:as.longitudinal]{as.longitudinal}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A      <- createA(p, "chain")

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

# convert data to a longitudinal-object
Ylongitudinal <- array2longitudinal(Y)

# convert back to array-format
Yback <- longitudinal2array(Ylongitudinal)
}

