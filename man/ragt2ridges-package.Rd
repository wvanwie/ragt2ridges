\name{ragt2ridges-package}
\alias{ragt2ridges-package}
\alias{ragt2ridges}
\docType{package}
\title{
Ridge Estimation of Vector Auto-Regressive (VAR) Processes
}
\description{
Ridge maximum likelihood estimation of vector auto-regressive processes and supporting functions for their exploitation. Currently, it includes: 
\itemize{
\item Ridge estimation of the parameters of the first-order Vector Auto-Regressive model, commonly referred to as the VAR(1) model, through the function \code{\link{ridgeVAR1}}. This function is complemented by \code{\link{optPenaltyVAR1}}, a function for penalty parameters selection through leave-one-out cross-validation (with supporting functions \code{\link{loglikVAR1}} and \code{\link{loglikLOOCVVAR1}}).

\item Functions for simulating VAR(1) data (\code{\link{createA}} and \code{\link{dataVAR1}}), data visualization (\code{\link{plotVAR1data}}), and some simple data manipulations (\code{\link{centerVAR1data}}, \code{\link{array2longitudinal}}, and \code{\link{longitudinal2array}}).

\item Ridge estimation of the precision matrix with known support (\code{\link{ridgePchordal}}, \code{\link{optPenaltyPchordal}}, and \code{\link{support4ridgeP}}).

\item Some diagnostics provided through \code{\link{evaluateVAR1fit}}, \code{\link{loglikLOOCVcontourVAR1}}, and \code{\link{ridgePathVAR1}}.

\item Several post-estimation analyses to exploit the fitted model. Among others: support determination of the VAR(1) model parameters (\code{\link{sparsifyVAR1}}), visualization of the (aspects of the) time-series chain graph (\code{\link{graphVAR1}} and \code{\link{CIGofVAR1}}), and summary statistics per variate in terms of the VAR(1) model and its associated time-series chain graph (\code{\link{nodeStatsVAR1}}, \code{\link{impulseResponseVAR1}}, and \code{\link{mutualInfoVAR1}}). 

\item Time-series omics data (\code{\link{hpvP53}}).
}
Future versions aim to include more functionality for time-series models.

The \code{ragt2ridges}-package is a sister-package to the \code{\link[rags2ridges:rags2ridges]{rags2ridges}}-package, augmenting the latter 'base' package with functionality for time-course studies. Being its sibling \code{ragt2ridges} mimicks \code{\link[rags2ridges:rags2ridges]{rags2ridges}} in the function names (compare e.g. \code{\link[rags2ridges:ridgeP]{ridgeP}} to \code{\link{ridgeVAR1}}). Moreover, some parts of the \code{ragt2ridges}-package not specific to time-series will probably be moved to the \code{\link[rags2ridges:rags2ridges]{rags2ridges}}-package in the near future, e.g. \code{\link{momentS}}, \code{\link{ridgePchordal}}, \code{\link{optPenaltyPchordal}}, and \code{\link{support4ridgeP}}.
}
\note{
The \code{(R)cpp}-code of \code{ragt2ridges} includes parts of the \code{Rcpp}-module of \code{\link[rags2ridges:rags2ridges]{rags2ridges}}, for it is currently impossible to import this directly.
}
\details{
\tabular{ll}{
Package: \tab ragt2ridges\cr
Type: \tab Package\cr
Version: \tab 0.1.15\cr
Date: \tab 2016-01-24\cr
License: \tab GPL (>= 2) \cr
}
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.ml>
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191.
}
\keyword{ package }
\seealso{
The \code{rags2ridges}-package.
}

