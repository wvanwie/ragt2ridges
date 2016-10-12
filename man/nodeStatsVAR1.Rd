\name{nodeStatsVAR1}
\alias{nodeStatsVAR1}
\title{
VAR(1) model node statistics
}
\description{
Function that calculates for each variate various statistics from a sparse VAR(1) model
}
\usage{
nodeStatsVAR1(sparseA, sparseP, as.table = FALSE)
}
\arguments{
\item{sparseA}{ A \code{matrix} \eqn{\mathbf{A}} of auto-regression parameters, which is assumed to be sparse. }
\item{sparseP}{ Precision \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon}} the error, which is assumed to be sparse. } 
\item{as.table}{ A \code{logical} indicating if the output should be in tabular format. }
}
\details{
The function calculates various node statistics from a sparse VAR(1) model. The input matrices \eqn{\mathbf{A}} and \eqn{\mathbf{\Omega}_{\varepsilon}} 
are assumed to be sparse regression coefficient and error precision matrices. From these matrices the global and contemporaneous conditional independence graphs are obtained (Dahlhaus, 2000; Dahlhaus, Eichler, 2003). 

For both graph types the function calculates various measures of centrality: node degree, betweenness centrality, closeness centrality, and eigenvalue
centrality. It also calculates the number of positive and the number of negative edges for each node. For more information on network measures, consult, e.g., Newman (2010).

In addition, for each variate the mutual information (with all other variates) (Miok et al., 2016), mean impulse response (Hamilton, 1994; Lutkepohl 2005), the (error) variance, and the partial error variance are presented. 
}
\value{
An object of class \code{list} (when \code{as.table = FALSE}) with slots:
\item{degreeAin}{ A \code{numeric} vector with the number of (temporal) edges pointing to each node ('in'-degree). }
\item{degreeAout}{ A \code{numeric} vector with the number of (temporal) edges leaving each node ('out'-degree). }
\item{nNegAin}{ A \code{numeric} vector with the number of negative (temporal) edges pointing to each node. }
\item{nPosAin}{ A \code{numeric} vector with the number of positive (temporal) edges pointing to each node ('in'-degree) }
\item{nNegAout}{ A \code{numeric} vector with the number of negative (temporal) edges leaving each node ('out'-degree) }
\item{nPosAout}{ A \code{numeric} vector with the number of positive (temporal) edges leaving each node ('out'-degree) }
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

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2016), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, accepted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>, Carel F.W. Peeters.
}
\seealso{
\code{\link{ridgeVAR1}}, \code{\link{sparsifyVAR1}}, \code{\link{graphVAR1}}
}
\examples{
# specify VAR(1) model parameters
A <- matrix(c(-0.1, -0.3, 0, 0.5, 0, 0, 0, 0, -0.4), byrow=TRUE, ncol=3)
P <- matrix(c(1, 0.5, 0, 0.5, 1, 0, 0, 0, 1), byrow=TRUE, ncol=3)

# adjacency matrix of (global) conditional independencies.
nodeStatsVAR1(A, P)
}
