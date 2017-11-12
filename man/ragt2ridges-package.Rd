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
\item Ridge estimation of the parameters of Vector Auto-Regressive models, commonly referred to as VAR models, through the functions \code{\link{ridgeVAR1}}, \code{\link{ridgeVAR2}}, \code{\link{ridgeVARX1}} and \code{\link{ridgeVAR1fused}}. These functions are complemented by \code{\link{optPenaltyVAR1}}, \code{\link{optPenaltyVAR2}}, \code{\link{optPenaltyVARX1}} and \code{\link{optPenaltyVAR1fused}}, functions for penalty parameters selection through (leave-one-out) cross-validation (with supporting functions \code{\link{loglikVAR1}}, \code{\link{loglikLOOCVVAR1}},
\code{\link{loglikLOOCVVAR2}}, \code{\link{loglikLOOCVVARX1}} and \code{\link{loglikLOOCVVAR1fused}}).

\item Functions for simulating VAR-type data (\code{\link{createA}}, \code{\link{dataVAR1}}, \code{\link{dataVAR2}} and \code{\link{dataVARX1}}), data visualization (\code{\link{plotVAR1data}}), and some simple data manipulations (\code{\link{centerVAR1data}}, \code{\link{array2longitudinal}}, and \code{\link{longitudinal2array}}).

\item Some diagnostics provided through \code{\link{evaluateVAR1fit}}, \code{\link{loglikLOOCVcontourVAR1}}, and \code{\link{ridgePathVAR1}}.

\item Several post-estimation analyses to exploit the fitted model. Among others: support determination of the various VAR model parameters (\code{\link{sparsifyVAR1}}, \code{\link{sparsifyVAR2}}, \code{\link{sparsifyVARX1}}), visualization of the (aspects of the) time-series chain graph (\code{\link{graphVAR1}}, \code{\link{graphVAR2}}, \code{\link{graphVARX1}}, \code{\link{CIGofVAR1}} and \code{\link{CIGofVAR2}}), and summary statistics per variate in terms of the VAR(1) model and its associated time-series chain graph (\code{\link{nodeStatsVAR1}}, \code{\link{impulseResponseVAR1}}, and \code{\link{mutualInfoVAR1}}). The latter are also available for the VAR(2) and VARX(1) models: \code{\link{impulseResponseVAR2}}, \code{\link{impulseResponseVARX1}} and \code{\link{mutualInfoVAR2}}).

\item Time-series omics data (\code{\link{hpvP53}}).
}
Future versions aim to include more functionality for time-series models.

The \code{ragt2ridges}-package is a sister-package to the \code{\link[rags2ridges:rags2ridges]{rags2ridges}}-package, augmenting the latter 'base' package with functionality for time-course studies. Being its sibling \code{ragt2ridges} mimicks \code{\link[rags2ridges:rags2ridges]{rags2ridges}} in the function names (compare e.g. \code{\link[rags2ridges:ridgeP]{ridgeP}} to \code{\link{ridgeVAR1}}).
}
\note{
The \code{(R)cpp}-code of \code{ragt2ridges} includes parts of the \code{Rcpp}-module of \code{\link[rags2ridges:rags2ridges]{rags2ridges}}, for it is currently impossible to import this directly.
}
\details{
\tabular{ll}{
Package: \tab ragt2ridges\cr
Type: \tab Package\cr
Version: \tab 0.3.0\cr
Date: \tab 2017-11-02\cr
License: \tab GPL (>= 2) \cr
}
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.ml>
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of network models from time-course omics data", \emph{submitted}. 
}
\keyword{ package }
\seealso{
The \code{\link[rags2ridges:rags2ridges]{rags2ridges}}-package.
}

