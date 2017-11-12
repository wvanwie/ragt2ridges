\name{graphVARX1}
\alias{graphVARX1}
\title{
Graphs of the temporal (or contemporaneous) relations implied by the VARX(1) model
}
\description{
Graphs of the (temporal and/or contemporaneous) relations among the variates as implied by the parameters of the VARX(1) model.
}
\usage{
graphVARX1(sparseA, sparseB, sparseP, type="TSCG", side="right", 
           prune=TRUE, nNamesY=NULL, nNamesX=NULL, main=NULL, 
           vertex.color.X="lightcyan2", vertex.color.T0="lightcyan2", 
           vertex.color.T1="lightcyan2", vertex.frame.color="steelblue", 
           vertex.label.cex=-1, vertex.label.color.X="black", 
           vertex.label.color.T0="black", vertex.label.color.T1="black", 
           vertex.label.font=1.5, vertex.size=-1, edge.arrow.size=-1, 
           edge.width=-1, ...)
}
\arguments{
\item{sparseA}{ A \code{matrix} \eqn{\mathbf{A}} of lag one autoregression parameters, which is assumed to be sparse. }
\item{sparseB}{ A \code{matrix} \eqn{\mathbf{B}} of regression parameters of the time-varying covariates, which is assumed to be sparse. }
\item{sparseP}{ Precision \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon}} the error, which is assumed to be sparse. }
\item{type}{ A \code{character} indicating what should be plotted. If \code{type="TSCG"}. the time series chain graph is plotted, while if \code{type="ABonly"} limits this graph
             to the temporal relations. If \code{type="contempPC"}, the global or contemporaneous (respectively) partial correlation graph is plotted. }
\item{side}{ A \code{character} indicating whether the contemporaneous dependencies should be plotted on the \code{bottomleft} (time t) or the \code{right}-hand (time t+1) side. Only active when \code{type="TSCG"}. }
\item{prune}{  A \code{logical} indicating whether to remove covariates without any temporal (or contemporaneous) relations (as implied by \code{sparseA}, \code{sparseB} and \code{sparseP}). }
\item{nNamesY}{ A \code{character} containing the variate names to be written inside the nodes. }
\item{nNamesX}{ A \code{character} containing the covariate names to be written inside the nodes. }
\item{main}{ The \code{character} to be plotted as title above the graph. }
\item{vertex.color.X}{ Color of covariate nodes. }
\item{vertex.color.T0}{ Color of nodes at time point t. This is ignored when \code{type="contempPC"}. }
\item{vertex.color.T1}{ Color of nodes at time point t+1. This is ignored when \code{type="contempPC"}. }
\item{vertex.frame.color}{ Refer to \code{plot.igraph}. }
\item{vertex.label.cex}{ Refer to \code{plot.igraph}. }
\item{vertex.label.color.X}{ Color of the covariate node label. }
\item{vertex.label.color.T0}{ Color of the node label at time point t. Ignored when \code{type="contempPC"}. }
\item{vertex.label.color.T1}{ Color of the node label at time point t+1. Ignored when \code{type="contempPC"}. }
\item{vertex.label.font}{ Refer to \code{plot.igraph}. }
\item{vertex.size}{ Refer to \code{plot.igraph}. }
\item{edge.arrow.size}{ Refer to \code{plot.igraph}. }
\item{edge.width}{ Refer to \code{plot.igraph}. }
\item{...}{ Other arguments to be passed on to \code{plot.igraph}. }
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of network models from time-course omics data", \emph{submitted}. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link{ridgeVARX1}}, \code{\link[ragt2ridges:graphVAR1]{graphVAR1}}, \code{\link[ragt2ridges:graphVAR2]{graphVAR2}}, \code{\link{sparsifyVARX1}}, \code{\link[igraph:plot.igraph]{plot.igraph}}.
}
\examples{
# specify VARX(1) model parameters
A <- matrix(c(-0.1, -0.3,    0,  0.5,    0,    0, 
                 0,    0, -0.4, -0.1, -0.3,    0,  
               0.5,    0,    0,    0,    0, -0.4, 
              -0.1, -0.3,    0,  0.5,    0,    0, 
                 0,    0, -0.4, -0.1, -0.3,    0, 
               0.5,    0,    0,    0,    0, -0.4), byrow=TRUE, ncol=6)
B <- matrix(c(   0,    0,  0.2,    0,  0.3, -0.7, 
                 0,    0,    0,  0.4,    0,    0, 
                 0, -0.3,    0,    0,    0,    0, 
                 0,    0,  0.5,    0,  0.1,    0, 
                 0,    0,    0,    0,  0.4,    0, 
                 0,    0,    0,  0.4,    0,    0), byrow=TRUE, ncol=6)
P <- matrix(c(   2,    0, -0.5,  0.5,    0,  0.5, 
                 0,    1,  0.5,  0.5,  0.5,    0, 
              -0.5,  0.5,    1,    0,    0,  0.5, 
               0.5,  0.5,    0,    1,    0,    0, 
                 0,  0.5,    0,    0,    1,    0, 
               0.5,    0,  0.5,    0,    0,    1), byrow=TRUE, ncol=6)

# time-series chain graph  of the VARX(1) model
graphVARX1(A, B, P, type="TSCG")
}

