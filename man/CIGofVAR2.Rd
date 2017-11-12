\name{CIGofVAR2}
\alias{CIGofVAR2}
\title{
Conditional independence graphs of the VAR(2) model
}
\description{
Constructs the global or contemporaneous conditional independence graph (CIG) of the VAR(2) model, as implied by the partial correlations.
}
\usage{
CIGofVAR2(sparseA1, sparseA2, sparseP, type="global")
}
\arguments{
\item{sparseA1}{ A \code{matrix} \eqn{\mathbf{A}_1} of lag one autoregression parameters, which is assumed to be sparse. }
\item{sparseA2}{ A \code{matrix} \eqn{\mathbf{A}_2} of lag-two autoregression parameters, which is assumed to be sparse. }
\item{sparseP}{ Precision \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon}} the error, which is assumed to be sparse. }
\item{type}{ A \code{character} indicating whether the \code{global} or \code{contemp} (contemporaneous) CIG should be plotted. }
}
\references{
Dahlhaus (2000), "Graphical interaction models for multivariate time series", \emph{Metrika}, 51, 157-172.

Dahlhaus, Eichler (2003), "Causality and graphical models in time series analysis", \emph{Oxford Statistical Science Series}, 115-137.

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of network models from time-course omics data", \emph{submitted}.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link{CIGofVAR1}},  \code{\link{graphVAR2}}, \code{\link[rags2ridges:sparsify]{sparsify}}, \code{\link{sparsifyVAR2}}.
}
\examples{
# specify VAR(2) model parameters
A1 <- matrix(c(-0.1, -0.3,    0, 
                0.5,    0,    0, 
                  0,    0, -0.4), byrow=TRUE, ncol=3)
A2 <- matrix(c(   0,    0,    0, 
                  0,    0,  0.2, 
                  0, -0.4,    0), byrow=TRUE, ncol=3)
P  <- matrix(c(   1,  0.5,    0, 
                0.5,    1,    0,  
                  0,    0,    1), byrow=TRUE, ncol=3)

# adjacency matrix of (global) conditional independencies.
CIGofVAR2(A1, A2, P, type="global")
}

