\name{sparsifyVAR2}
\alias{sparsifyVAR2}
\title{
Function that determines the support of autoregression parameters of the VAR(2) model.
}
\description{
Function that determines the null and non-null elements of \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2}, the matrices of lag one and two (respectively) autoregression coefficients. 
}
\usage{
sparsifyVAR2(A1, A2, SigmaE, threshold=c("absValue", "localFDR", "top"), 
             absValueCut=c(0.25, 0.25), FDRcut=c(0.8, 0.8), top=c(10,10), 
             zerosA1=matrix(nrow=0, ncol=2), 
             zerosA2=matrix(nrow=0, ncol=2), 
             statistics=FALSE, verbose=TRUE)
}
\arguments{
\item{A1}{ A \code{matrix} \eqn{\mathbf{A}_1} of lag one autoregression parameters. }
\item{A2}{ A \code{matrix} \eqn{\mathbf{A}_2} of lag two autoregression parameters. }
\item{SigmaE}{ Covariance \code{matrix} of the errors (innovations). }
\item{threshold}{ A \code{character} signifying type of sparsification of both \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2} by thresholding. Must be one of: "absValue", "localFDR", or "top". }
\item{absValueCut}{ A \code{numeric} (of length two) giving the cut-offs for element selection based on absolute value thresholding. The elements are the cut-offs applied to \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2}, respectively. }
\item{FDRcut}{ A \code{numeric} (of length two) giving the cut-offs for element selection based on local false discovery rate (FDR) thresholding. The elements are the cut-offs applied to \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2}, respectively. }
\item{top}{ A \code{numeric} (of length two) specifying the numbers of elements of \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2}, respectively, which are to be selected. }
\item{zerosA1}{ A \code{matrix} with indices of entries of \eqn{\mathbf{A}_1} that are (prior to sparsification) known to be zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{A}_1}. The first column contains the row indices and the second the column indices. }
\item{zerosA2}{ A \code{matrix} with indices of entries of \eqn{\mathbf{A}_2} that are (prior to sparsification) known to be zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{A}_2}. The first column contains the row indices and the second the column indices. }
\item{statistics}{ A \code{logical} indicator: should test statistics be returned. This only applies when 
\code{threshold = "localFDR"}}
\item{verbose}{ A \code{logical} indicator: should intermediate output be printed on the screen? }
}
\value{
A \code{list}-object with slots:
\item{zerosA1}{ \code{Matrix} with indices of entries of \eqn{\mathbf{A}_1} that are identified to be null. It includes the elements of \eqn{\mathbf{A}_1} assumed to be zero prior to the sparsification as specified through the \code{zerosA1} option. }
\item{nonzerosA1}{ \code{Matrix} with indices of entries of \eqn{\mathbf{A}_1} that are identified to be non-null. }
\item{statisticsA1}{ \code{Matrix} with test statistics employed in the local FDR procedure. }
\item{zerosA2}{ \code{Matrix} with indices of entries of \eqn{\mathbf{A}_2} that are identified to be null. It includes the elements of \eqn{\mathbf{A}_2} assumed to be zero prior to the sparsification as specified through the \code{zerosA2} option. }
\item{nonzerosA2}{ \code{Matrix} with indices of entries of \eqn{\mathbf{A}_2} that are identified to be non-null. }
\item{statisticsA2}{ \code{Matrix} with test statistics employed in the local FDR procedure. }
The matrices \code{zerosA1}, \code{nonzerosA1}, \code{zerosA2} and \code{nonzerosA2} comprise two columns, each row corresponding to an entry of \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2}, respectively. The first column contains the row indices and the second the column indices.
}
\details{
When \code{threshold = "localFDR"} the function, following Lutkepohl (2005), divides the elements of (possibly regularized) input matrix \eqn{\mathbf{A}_1} (or \eqn{\mathbf{A}_2}) of lag one (or two) autoregression coefficients by (an approximation of) their standard errors. Subsequently, the support of the matrix \eqn{\mathbf{A}_1} (or \eqn{\mathbf{A}_2}) is determined by usage of the local FDR. In that case a mixture model is fitted to the nonredundant (standardized) elements of \eqn{\mathbf{A}_1} (or \eqn{\mathbf{A}_2}) by \code{\link[fdrtool:fdrtool]{fdrtool}}. The decision to retain elements is then based on the argument \code{FDRcut}. Elements with a posterior probability \eqn{>=q}  \code{FDRcut} (equalling 1 - local FDR) are retained. See Strimmer (2008) for \code{\link[rags2ridges:sparsify]{sparsify}}further details. Alternatively, the support of \eqn{\mathbf{A}_1} (or \eqn{\mathbf{A}_2})  is determined by simple thresholding on the absolute values of matrix entries (\code{threshold = "absValue"}). A third option (\code{threshold = "top"}) is to retain a prespecified number of matrix entries based on absolute values of the elements of \eqn{\mathbf{A}_1} (or \eqn{\mathbf{A}_2}). For example, one could wish to retain those entries representing the ten strongest cross-temporal coefficients. 
     
The argument \code{absValueCut} is only used when \code{threshold = "absValue"}. The argument \code{FDRcut} is only used when \code{threshold = "localFDR"}. The argument \code{top} is only used when \code{threshold = "top"}. 

When prior to the sparsification knowledge on the support of \eqn{\mathbf{A}_1} (or \eqn{\mathbf{A}_2}) is specified through the option \code{zerosA1} (or \code{zerosA2}), the corresponding elements of \eqn{\mathbf{A}_1} (or \eqn{\mathbf{A}_2}) are then not taken along in the local FDR procedure.
}
\references{
Lutkepohl, H. (2005), New Introduction to Multiple Time Series Analysis. Springer, Berlin.

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of network models from time-course omics data", \emph{submitted}. 

Strimmer, K. (2008), ``fdrtool: a versatile R package for estimating local and tail area-based false discovery rates'', \emph{Bioinformatics} 24(12): 1461-1462.

Van Wieringen, W.N., Peeters, C.F.W. (2016), ``Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data'', \emph{Computational Statistics and Data Analysis}, 103, 284-303.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>, Carel F.W. Peeters.
}
\seealso{
\code{\link{ridgeVAR2}}, \code{\link[rags2ridges:sparsify]{sparsify}}, \code{\link{sparsifyVAR1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A1 <- -createA(p, "clique", nCliques=1, nonzeroA=0.1)
A2 <- t(createA(p, "chain", nBands=1, nonzeroA=0.1))

# generate data
Y <- dataVAR2(n, T, A1, A2, SigmaE)

# fit VAR(1) model
VAR2hat <- ridgeVAR2(Y, 1, 1, 1)

# obtain support of adjacancy matrix
A1nullornot <- matrix(0, p, p)
A1nullornot[sparsifyVAR2(VAR2hat$A1, VAR2hat$A2, solve(VAR1hat$P), 
                         threshold="top", top=c(3,3))$nonzerosA1] <- 1

## plot non-null structure of A1 
edgeHeat(A1nullornot)
}

