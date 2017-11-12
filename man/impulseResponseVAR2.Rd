\name{impulseResponseVAR2}
\alias{impulseResponseVAR2}
\title{
Impulse response analysis of the VAR(2) model
}
\description{
Evaluate the impulse responses of the VAR(2) model. It assesses the effect of an innovation (error) at one time point on the variates at future time points. In the VAR(2) model this amounts evaluating a recursive relationship in \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2}, the matrices of lag 1 and lag 2 autoregression coefficients.
}
\usage{
impulseResponseVAR2(A1, A2, T)
}
\arguments{
\item{A1}{ A \code{matrix} \eqn{\mathbf{A}_1} with lag 1 autoregression parameters. }
\item{A2}{ A \code{matrix} \eqn{\mathbf{A}_2} with lag 2 autoregression parameters. }
\item{T}{ Non-negative \code{integer} of length one specifying the time point at which the impulse responses is to be evaluated. }
}
\value{
A \code{matrix} with the impulse response of the innovation vector.
}
\references{
Hamilton, J. D. (1994). \emph{Time series analysis}. Princeton: Princeton university press.

Lutkepohl, H. (2005). \emph{New Introduction to Multiple Time Series Analysis}. Springer, Berlin.

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link{impulseResponseVAR1}}, \code{\link{impulseResponseVARX1}}, \code{\link{ridgeVAR2}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A1     <- createA(p, "clique", nCliques=1)
A2     <- createA(p, "hub", nHubs=1)

# generate time-varying covariates in accordance with VAR(2) process
Y <- dataVAR2(n, T, A1, A2, SigmaE)

# fit VAR(2) model
VAR2hat <- ridgeVAR2(Y, 1, 1, 1)

# impulse response analysis
impulseResponseVAR2(VAR2hat$A1, VAR2hat$A2, 10)
}

