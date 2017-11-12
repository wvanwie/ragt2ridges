\name{impulseResponseVARX1}
\alias{impulseResponseVARX1}
\title{
Impulse response analysis of the VARX(1) model
}
\description{
Evaluate the impulse responses of the VARX(1) model. It assesses the effect of an innovation (error) at one time point on the variates at future time points. In the VARX(1) model this amounts to studying powers of \eqn{\mathbf{A}}, the matrix of autoregression coefficients, for the innovation. For the impulse response of the time-varying covariates these powers are post-multiplied by \eqn{\mathbf{B}}. 
}
\usage{
impulseResponseVARX1(A, B, T)
}
\arguments{
\item{A}{ A \code{matrix} \eqn{\mathbf{A}} with autoregression parameters. }
\item{B}{ A \code{matrix} \eqn{\mathbf{B}} with regression parameters time-varying covariates. }
\item{T}{ Non-negative \code{integer} of length one specifying the time point at which the impulse responses is to be evaluated. }
}
\value{
A \code{list}-object with slots:
\item{\code{impResponse}}{ A \code{matrix} with the impulse response of the innovation vector.}
\item{\code{impResponseX}}{ A \code{matrix} with the impulse response of the time-varying covariate vector. }
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
\code{\link{ridgeVAR1}}, \code{\link{ridgeVAR2}}, \code{\link{ridgeVARX1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
Ax     <- createA(p, "chain", nBands=1)

# generate time-varying covariates in accordance with VAR(1) process
X <- dataVAR1(n, T, Ax, SigmaE)

# set model parameters
B <- createA(p, "clique", nCliques=1)
A <- createA(p, "hub", nHubs=1)

# generate time-varying covariates in accordance with VAR(1) process
Y <- dataVARX1(X, A, B, SigmaE, lagX=0)

# fit VARX(1) model
VARX1hat <- ridgeVARX1(Y, X, 1, 1, 1, lagX=0)

# impulse response analysis
impulseResponseVARX1(VARX1hat$A, VARX1hat$B, 10)
}

