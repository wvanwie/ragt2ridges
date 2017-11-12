\name{CIGofVAR1}
\alias{CIGofVAR1}
\title{
Conditional independence graphs of the VAR(1) model
}
\description{
Constructs the global or contemporaneous conditional independence graph (CIG) of the VAR(1) model, as implied by the partial correlations.
}
\usage{
CIGofVAR1(sparseA, sparseP, type="global")
}
\arguments{
\item{sparseA}{ A \code{matrix} \eqn{\mathbf{A}} of autoregression parameters, which is assumed to be sparse. }
\item{sparseP}{ Precision \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon}} the error, which is assumed to be sparse. }
\item{type}{ A \code{character} indicating whether the \code{global} or \code{contemp} (contemporaneous) CIG should be plotted. }
}
\references{
Dahlhaus (2000), "Graphical interaction models for multivariate time series", \emph{Metrika}, 51, 157-172.

Dahlhaus, Eichler (2003), "Causality and graphical models in time series analysis", \emph{Oxford Statistical Science Series}, 115-137.

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link{graphVAR1}}, \code{\link[rags2ridges:sparsify]{sparsify}}, \code{\link{sparsifyVAR1}}.
}
\examples{
# specify VAR(1) model parameters
A <- matrix(c(-0.1, -0.3,    0, 
               0.5,    0,    0, 
                 0,    0, -0.4), byrow=TRUE, ncol=3)
P <- matrix(c(   1,  0.5,    0, 
               0.5,    1,    0,  
                 0,    0,    1), byrow=TRUE, ncol=3)

# adjacency matrix of (global) conditional independencies.
CIGofVAR1(A, P, type="global")
}

