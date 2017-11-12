\name{ridgePathVAR1}
\alias{ridgePathVAR1}
\title{
Visualize the ridge regularization paths of the parameters of the VAR(1) model
}
\description{
Function that visualizes the regularization paths of the parameters of the VAR(1) model. The elements of the ridge ML estimate of either \eqn{\mathbf{A}} or (possibly standardized, inverse of) 
\eqn{\mathbf{\Sigma}_{\varepsilon}} are plotted against a specified range of their penalty parameter (keeping the other penalty parameter fixed). 
}
\usage{
ridgePathVAR1(Y, lambdaAgrid, lambdaPgrid, pathType="A", 
              plotTypeSigmaE="pcor", diag=FALSE, verbose=TRUE, ...)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{lambdaAgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to the \eqn{\lambda_a} (the penalty parameter for the autoregression coefficient matrix \eqn{\mathbf{A}}). }
\item{lambdaPgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to the \eqn{\lambda_{\omega}} (the penalty parameters for the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}). }
\item{pathType}{ A \code{character} indicating of which parameter to plot its ridge regularization paths. Either "\code{A}" or "\code{SigmaE}". } 
\item{plotTypeSigmaE}{ A \code{character} indicating the type of element for which a visualization of the regularization paths (of \eqn{\mathbf{\Sigma}_{\varepsilon}}) is desired. Must
be one of: "\code{pcor}", "\code{cor}", "\code{cov}", "\code{prec}". }
\item{diag}{ A \code{logical} indicating if the diagonal elements should be retained for visualization of the regularization path of \eqn{\mathbf{\Sigma}_{\varepsilon}}. }
\item{verbose}{ A \code{logical} indicator: should intermediate output be printed on the screen? }
\item{...}{  Other arguments to be passed to \code{ridgeVAR1}. }
}
\details{
If \code{pathType="A"}, the regularization paths of \eqn{\mathbf{A}} will be evaluated for \eqn{\lambda_{\omega}} equal to the minimum value of \code{lambdaPgrid}. If \code{pathType="SigmaE"}, the regularization paths of (inverse of / possibly standardarized) \eqn{\mathbf{\Sigma}_{\varepsilon}} will be evaluated for \eqn{\lambda_a} equal to the minimum value of \code{lambdaAgrid}.

Regularization paths may be visualized for (partial) correlations, covariances and precision elements. The type of element for which a visualization of the regularization paths is desired can be indicated by the argument \code{plotType}.

The arguments \code{diag} and \code{plotTypeSigmaE} are ignored when \code{pathType="A"}.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>.
}
\seealso{
\code{\link{ridgePathS}}, \code{\link[rags2ridges:ridgeP]{ridgeP}}, \code{\link{ridgeVAR1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10
     
# set model parameters
SigmaE <- diag(p)/4
A      <- createA(p, "chain")
     
# generate data
Y <- dataVAR1(n, T, A, SigmaE)

## Visualize regularization paths
lambdaAgrid <- seq(0.01, 1, length.out=20) 
lambdaPgrid <- seq(0.01, 1, length.out=20) 
ridgePathVAR1(Y, lambdaAgrid, lambdaPgrid, pathType="A")
}
