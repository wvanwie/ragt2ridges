\name{array2longitudinal}
\alias{array2longitudinal}
\title{
Convert a time-series array to a longitudinal-object.
}
\description{
Converts a 3-dim array (containing time-series data of multiple individuals) to an object of that can directly be converted to the \code{\link[longitudinal:longitudinal]{longitudinal}} class through the \code{\link[longitudinal:as.longitudinal]{as.longitudinal}} function.
}
\usage{
array2longitudinal(Y, keepMissings=TRUE)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. }
\item{keepMissings}{ The \code{array}-format assumes a balanced layout of the time-course experiment. The experiment may have failed for some design points and no data is available. In the \code{longitudinal}-format these design points may be left out. This \code{logical} indicates whether they should be kept in (or not). }
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\value{
A \code{longitudinal}-object containing time-series data.
}
\seealso{
\code{\link{dataVAR1}}, \code{\link[longitudinal:longitudinal]{longitudinal}}, \code{\link{longitudinal2array}}, \code{\link[longitudinal:as.longitudinal]{as.longitudinal}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A      <- createA(p, "chain")

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

# convert to longitudinal
Ylongitudinal <- array2longitudinal(Y)
}

