\name{graphVAR1}
\alias{graphVAR1}
\title{
Graphs of the temporal (or contemporaneous) relations implied by the VAR(1) model
}
\description{
Graphs of the (conditional independence, temporal, or contemporaneous) relations among the variates as implied by the parameters of the VAR(1) model.
}
\usage{
graphVAR1(sparseA, sparseP, type="TSCG", side="left", prune=TRUE, 
          nNames=NULL, main=NULL, vertex.color.T0="lightcyan2", 
          vertex.color.T1="lightcyan2", vertex.frame.color="steelblue", 
          vertex.label.cex=-1, vertex.label.color.T0="black", 
          vertex.label.color.T1="black", vertex.label.font=1.5, 
          vertex.size=-1, edge.arrow.size=-1, edge.width=-1, ...)
}
\arguments{
\item{sparseA}{ A \code{matrix} \eqn{\mathbf{A}} of autoregression parameters, which is assumed to be sparse. }
\item{sparseP}{ Precision \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon}} the error, which is assumed to be sparse. }
\item{type}{ A \code{character} indicating what should be plotted. If \code{type="TSCG"}. the time series chain graph is plotted, while if \code{type="Aonly"} limits this graph
             to the temporal relations. If \code{type="globalPC"} or \code{type="contempPC"}, the global or contemporaneous (respectively) partial correlation graph is plotted. }
\item{side}{ A \code{character} indicating whether the contemporaneous dependencies should be plotted on the left-hand (time t) or the right-hand (time t+1) side. Only active when \code{type="TSCG"}. }
\item{prune}{  A \code{logical} indicating whether to remove covariates without any temporal (or contemporaneous) relations (as implied by \code{sparseA} and \code{sparseP}). }
\item{nNames}{ A \code{character} containing the covariate names to be written inside the nodes. }
\item{main}{ The \code{character} to be plotted as title above the graph. }
\item{vertex.color.T0}{ Color of nodes at time point t. }
\item{vertex.color.T1}{ Color of nodes at time point t+1. This is ignored when \code{type="globalPC"} or \code{type="contempPC"}. }
\item{vertex.frame.color}{ Refer to \code{plot.igraph}. }
\item{vertex.label.cex}{ Refer to \code{plot.igraph}. }
\item{vertex.label.color.T0}{ Color of the node label at time point t. }
\item{vertex.label.color.T1}{ Color of the node label at time point t+1. Ignored when \code{type="globalPC"} or \code{type="contempPC"}. }
\item{vertex.label.font}{ Refer to \code{plot.igraph}. }
\item{vertex.size}{ Refer to \code{plot.igraph}. }
\item{edge.arrow.size}{ Refer to \code{plot.igraph}. }
\item{edge.width}{ Refer to \code{plot.igraph}. }
\item{...}{ Other arguments to be passed on to \code{plot.igraph}. }
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link{CIGofVAR1}}, \code{\link[rags2ridges:sparsify]{sparsify}}, \code{\link{sparsifyVAR1}}, \code{\link[igraph:plot.igraph]{plot.igraph}}.
}
\examples{
# specify VAR(1) model parameters
A <- matrix(c(-0.1, -0.3,    0,  0.5,    0,    0, 
                 0,    0, -0.4, -0.1, -0.3,    0,  
               0.5,    0,    0,    0,    0, -0.4, 
              -0.1, -0.3,    0,  0.5,    0,    0, 
                 0,    0, -0.4, -0.1, -0.3,    0, 
               0.5,    0,    0,    0,    0, -0.4), byrow=TRUE, ncol=6)
P <- matrix(c(   2,    0, -0.5,  0.5,    0,  0.5, 
                 0,    1,  0.5,  0.5,  0.5,    0, 
              -0.5,  0.5,    1,    0,    0,  0.5, 
               0.5,  0.5,    0,    1,    0,    0, 
                 0,  0.5,    0,    0,    1,    0, 
               0.5,    0,  0.5,    0,    0,    1), byrow=TRUE, ncol=6)

# adjacency matrix of (global) conditional independencies.
graphVAR1(A, P, type="TSCG")
}

