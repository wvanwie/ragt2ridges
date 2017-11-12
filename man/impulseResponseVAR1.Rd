\name{impulseResponseVAR1}
\alias{impulseResponseVAR1}
\title{
Impulse response analysis of the VAR(1) model
}
\description{
Evaluate the impulse responses of the VAR(1) model. It assesses the effect of an innovation (error) at one time point on the 
variates at future time points. In the VAR(1) model this amounts to studying powers of \eqn{\mathbf{A}}, the matrix of autoregression coefficients. 
}
\usage{
impulseResponseVAR1(A, T)
}
\arguments{
\item{A}{ A \code{matrix} \eqn{\mathbf{A}} of autoregression parameters. }
\item{T}{ Positive \code{numeric} of length one specifying the time points for which the impulse response is to be evaluated. }
}
\value{
Object of class \code{matrix}. Rows and columns correspond to covariates, elements to the impulse response of 'row variate' on the 'columns variate' on \code{T} time points from the current. 
}
\references{
Hamilton, J. D. (1994), Time series analysis. Princeton: Princeton university press.

Lutkepohl, H. (2005), New Introduction to Multiple Time Series Analysis. Springer, Berlin.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link{ridgeVAR1}}.
}
\examples{
# set dimensions
p <- 3
n <- 4
T <- 10

# set model parameters
SigmaE <- diag(p)/4
A      <- matrix(c(-0.1, -0.3, 0.6, 0.5, -0.4, 0, 0.3, -0.5, -0.2), 
                 byrow=TRUE, ncol=3)

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

# fit VAR(1) model
VAR1hat <- ridgeVAR1(Y, 1, 1)

# impulse response analysis
impulseResponseVAR1(VAR1hat$A, 10)
}

