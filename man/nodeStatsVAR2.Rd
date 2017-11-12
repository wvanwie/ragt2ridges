\name{nodeStatsVAR2}
\alias{nodeStatsVAR2}
\title{
VAR(2) model node statistics
}
\description{
Function that calculates for each variate various statistics from a sparse VAR(2) model
}
\usage{
nodeStatsVAR2(sparseA1, sparseA2, sparseP, as.table = FALSE)
}
\arguments{
\item{sparseA1}{ A \code{matrix} \eqn{\mathbf{A}_1} of autoregression parameters, which is assumed to be sparse. }
\item{sparseA2}{ A \code{matrix} \eqn{\mathbf{A}_2} of autoregression parameters, which is assumed to be sparse. }
\item{sparseP}{ Precision \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon}} the error, which is assumed to be sparse. } 
\item{as.table}{ A \code{logical} indicating if the output should be in tabular format. }
}
\details{
The function calculates various node statistics from a sparse VAR(2) model. The input matrices \eqn{\mathbf{A}_1}, \eqn{\mathbf{A}_2} and \eqn{\mathbf{\Omega}_{\varepsilon}} 
are assumed to be sparse autoregression coefficient and error precision matrices. From these matrices the global and contemporaneous conditional independence graphs are obtained (Dahlhaus, 2000; Dahlhaus, Eichler, 2003). 

For both graph types the function calculates various measures of centrality: node degree, betweenness centrality, closeness centrality, and eigenvalue
centrality. It also calculates the number of positive and the number of negative edges for each node. For more information on network measures, consult, e.g., Newman (2010).

In addition, for each variate the mutual information (with all other variates) (Miok et al., 2017), mean impulse response (Hamilton, 1994; Lutkepohl 2005), the (error) variance, and the partial error variance are presented. 
}
\value{
An object of class \code{list} (when \code{as.table = FALSE}) with slots:
\item{degreeA1in}{ A \code{numeric} vector with the number of one-lag (temporal) edges pointing to each node ('in'-degree). }
\item{degreeA1out}{ A \code{numeric} vector with the number of one-lag (temporal) edges leaving each node ('out'-degree). }
\item{nNegA1in}{ A \code{numeric} vector with the number of negative one-lag (temporal) edges pointing to each node. }
\item{nPosA1in}{ A \code{numeric} vector with the number of positive one-lag (temporal) edges pointing to each node ('in'-degree) }
\item{nNegA1out}{ A \code{numeric} vector with the number of negative one-lag (temporal) edges leaving each node ('out'-degree) }
\item{nPosA1out}{ A \code{numeric} vector with the number of positive one-lag (temporal) edges leaving each node ('out'-degree) }
\item{degreeA2in}{ A \code{numeric} vector with the number of two-lag (temporal) edges pointing to each node ('in'-degree). }
\item{degreeA2out}{ A \code{numeric} vector with the number of two-lag (temporal) edges leaving each node ('out'-degree). }
\item{nNegA2in}{ A \code{numeric} vector with the number of negative two-lag (temporal) edges pointing to each node. }
\item{nPosA2in}{ A \code{numeric} vector with the number of positive two-lag (temporal) edges pointing to each node ('in'-degree) }
\item{nNegA2out}{ A \code{numeric} vector with the number of negative two-lag (temporal) edges leaving each node ('out'-degree) }
\item{nPosA2out}{ A \code{numeric} vector with the number of positive two-lag (temporal) edges leaving each node ('out'-degree) }
\item{degreePe}{ A \code{numeric} vector with the number of contemporaneous edges of each node (as implied by the error precision matrix) }
\item{betweennessPe}{ A \code{numeric} vector representing the contemporaneous betweenness centrality for each node. }
\item{closenessPe}{ A \code{numeric} vector representing the contemporaneous closeness centrality for each node. }
\item{eigenCentralityPe}{ A \code{numeric} vector representing the contemporaneous eigen centrality for each node. }
\item{nNegPe}{ A \code{numeric} vector representing the number of negative contemporaneous edges for each node. }
\item{nPosPe}{ A \code{numeric} vector representing the number of positive contemporaneous edges for each node. }
\item{variancePe}{ A \code{numeric} vector representing the error variance of each node. }
\item{partialVarPe}{ A \code{numeric} vector representing the partial error variance of each node. }
\item{varianceY}{ A \code{numeric} vector representing the variance of each node. }
\item{degreePy}{ A \code{numeric} number of edges of each node in the global Markov graph. }
\item{betweennessPy}{ A \code{numeric} vector representing the betweenness centrality for each node in the global Markov graph. }
\item{closenessPy}{ A \code{numeric} vector representing the closeness centrality for each node in the global Markov graph. }
\item{eigenCentralityPy}{ A \code{numeric} vector representing the eigen centrality for each node in the global Markov graph. }
\item{mutualInfo_Tplus1}{ A \code{numeric} vector with for each node its mutual information with all other nodes at the next (t+1) time point. }
\item{mutualInfo_Tplus2}{ A \code{numeric} vector with for each node its mutual information with all other nodes at the (t+2)-th time point. }
\item{itemResponse_Tplus1}{ A \code{numeric} vector with for each node its  mean absolute impulse response on all other nodes at the next (t+1) time point. }
\item{itemResponse_Tplus2}{ A \code{numeric} vector with for each node its  mean absolute impulse response on all other nodes at the (t+2)-th time point. }
When \code{as.table = TRUE} the list items above are represented in tabular form as an object of 
class \code{matrix}.

Future versions of this function may include additional statistics
}
\references{
Dahlhaus (2000), "Graphical interaction models for multivariate time series", Metrika, 51, 157-172.

Dahlhaus, Eichler (2003), "Causality and graphical models in time series analysis", Oxford Statistical Science Series, 115-137.

Hamilton, J. D. (1994), Time series analysis. Princeton: Princeton university press.

Lutkepohl, H. (2005), New Introduction to Multiple Time Series Analysis. Springer, Berlin.

Newman, M.E.J. (2010). "Networks: an introduction", Oxford University Press.

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of network models from time-course omics data", \emph{submitted}. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>, Carel F.W. Peeters.
}
\seealso{
\code{\link{ridgeVAR2}}, \code{\link{sparsifyVAR2}}, \code{\link{graphVAR2}}
}
\examples{
# specify dimension
p <- 3

# specify VAR(1) model parameters
A1 <- -createA(p, "clique", nCliques=1, nonzeroA=0.1)
A2 <- t(createA(p, "chain", nBands=1, nonzeroA=0.1))
P <- matrix(c(  1, 0.5,  0, 
              0.5,   1,  0, 
                0,   0,  1), byrow=TRUE, ncol=3)

# adjacency matrix of (global) conditional independencies.
# nodeStatsVAR2(A1, A2, P)
}
