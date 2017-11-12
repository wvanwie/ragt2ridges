\name{mutualInfoVAR1}
\alias{mutualInfoVAR1}
\title{
Mutual information analysis of the VAR(1) model
}
\description{
Evaluate, within the VAR(1) model, the mutual information between each variate at the current time point and those at a future time point.
}
\usage{
mutualInfoVAR1(A, SigmaE, T)
}
\arguments{
\item{A}{ \code{Matrix} \eqn{\mathbf{A}} of autoregression parameters. }
\item{SigmaE}{ Covariance \code{matrix} of the errors (innovations). }
\item{T}{ Positive \code{numeric} of length one specifying the future time point with which the mutual informations are to be evaluated. }
}
\value{
Object of class \code{numeric} with elements corresponding to the mutual informations. The \eqn{j}-th element represents the mutual information of the \eqn{j}-th variate at the current time point with all variates at the \code{T}-th time point from now. 
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
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
A <- matrix(c(-0.1, -0.3,  0.6, 
               0.5, -0.4,    0, 
               0.3, -0.5, -0.2), byrow=TRUE, ncol=3)

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

# fit VAR(1) model
VAR1hat <- ridgeVAR1(Y, 1, 1)

# impulse response analysis
mutualInfoVAR1(VAR1hat$A, solve(symm(VAR1hat$P)), T=5)
}

