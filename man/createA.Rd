\name{createA}
\alias{createA}
\title{
Generation of the VAR(1) autoregression coefficient matrix.
}
\description{
Generates autoregression coefficient matrices of the VAR(1) with various type of topologies
}
\usage{
createA(p, topology, nonzeroA=0, nCliques=1, nHubs=1, 
        nBands=1, percZeros=0.9, stationary=TRUE)
}
\arguments{
\item{p}{ A positive \code{integer} specifying the dimension of the square matrix \eqn{\mathbf{A}}. }
\item{topology}{ Topology to impose on \eqn{\mathbf{A}}: a \code{character} equalling either \code{clique},  \code{hub},  \code{chain}, or \code{random}. }
\item{nonzeroA}{ Numeric, value that nonzero elements of \eqn{\mathbf{A}} will assume. If equal to zero, a random value from the interval [-1,1] is sampled. }
\item{nCliques}{ When \code{topology="clique"}, this positive integer specifies number of cliques. }
\item{nHubs}{ When \code{topology="hub"}, this positive integer specifies number of hubs. }
\item{nBands}{ When \code{topology="chain"}, this positive integer specifies number of bands. }
\item{percZeros}{ When \code{topology="random"}, the probability with which zero elements of \eqn{\mathbf{A}} are to be sampled. }
\item{stationary}{ A \code{logical}: should the generated \eqn{\mathbf{A}} be stationary? }
}
\value{
A \code{matrix} with autoregression coefficient matrix \eqn{\mathbf{A}} of the VAR(1) model.
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2016), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
}
\author{
Viktorian Miok, Wessel N. van Wieringen <w.vanwieringen@vumc.nl>.
}
\seealso{
\code{\link[rags2ridges:createS]{createS}}, \code{\link{dataVAR1}}.
}
\examples{
# create the VAR(1) parameters
A  <- createA(10, topology="clique", nonzeroA=0.1, nClique=4)
Se <- createS(1000, 10, "star")

# sample data from the VAR(1) model with above parameters
Y <- dataVAR1(4, 8, A, Se)
}




