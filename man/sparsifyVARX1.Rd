\name{sparsifyVARX1}
\alias{sparsifyVARX1}
\title{
Function that determines the support of (auto)regression parameters of the VARX(1) model.
}
\description{
Function that determines the support of \eqn{\mathbf{A}}, the matrix with autoregression coefficients, and \eqn{\mathbf{B}}, the matrix with regression coefficients of the time-varying covariates, of the VARX(1) model. 
}
\usage{
sparsifyVARX1(X, A, B, SigmaE, threshold=c("absValue", "localFDR", "top"), 
              absValueCut=rep(0.25, 2), FDRcut=rep(0.8, 2), top=rep(10, 2), 
              zerosA=matrix(nrow=0, ncol=2), zerosB=matrix(nrow=0, ncol=2), 
              statistics=FALSE, verbose=TRUE)
}
\arguments{
\item{X}{ Three-dimensional \code{array} containing the time-varying covariate data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{A}{ A \code{matrix} \eqn{\mathbf{A}} of autoregression parameters. }
\item{B}{ A \code{matrix} \eqn{\mathbf{B}} of regression parameters of the time-varying covariates stored in the array \code{X}. }
\item{SigmaE}{ Covariance \code{matrix} of the errors (innovations). }
\item{threshold}{ A \code{character} signifying type of sparsification of \eqn{\mathbf{A}} and \eqn{\mathbf{B}} by thresholding. Must be one of: "absValue", "localFDR", or "top". }
\item{absValueCut}{ A \code{numeric} of length two giving the cut-offs  for element selection based on absolute value thresholding. }
\item{FDRcut}{ A \code{numeric} of length two giving the cut-off (for \eqn{\mathbf{A}} and \eqn{\mathbf{B}}, respectively) for element selection based on local false discovery rate (FDR) thresholding. }
\item{top}{ A \code{numeric} of length two giving the number of elements of \eqn{\code{A}} and \eqn{\mathbf{B}}, respectively, which is to be selected, based on absolute value thresholding. }
\item{zerosA}{ \code{Matrix} with indices of entries of \eqn{\mathbf{A}} that are (prior to sparsification) known to be zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{A}}. The first column contains the row indices and the second the column indices. }
\item{zerosB}{ A \code{Matrix} with indices of entries of \eqn{\mathbf{B}} that are (prior to sparsification) known to be zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{B}}. The first column contains the row indices and the second the column indices. }
\item{statistics}{ A \code{Logical} indicator: should test statistics be returned. This only applies when 
\code{threshold = "localFDR"}}
\item{verbose}{ \code{Logical} indicator: should intermediate output be printed on the screen? }
}
\value{
A \code{list}-object with slots:
\item{zerosA}{ \code{Matrix} with indices of entries of \eqn{\mathbf{A}} that are identified to be null. It includes the elements of \eqn{\mathbf{A}} assumed to be zero prior to the sparsification as specified through the \code{zerosAknown} option. }
\item{nonzerosA}{ \code{Matrix} with indices of entries of \eqn{\mathbf{A}} that are identified to be non-null. }
\item{statisticsA}{ \code{Matrix} with test statistics employed in the local FDR procedure for the sparsification of \eqn{\mathbf{A}}. }
\item{zerosB}{ \code{Matrix} with indices of entries of \eqn{\mathbf{B}} that are identified to be null. It includes the elements of \eqn{\mathbf{B}} assumed to be zero prior to the sparsification as specified through the \code{zerosBknown} option. }
\item{nonzerosB}{ \code{Matrix} with indices of entries of \eqn{\mathbf{B}} that are identified to be non-null. }
\item{statisticsB}{ \code{Matrix} with test statistics employed in the local FDR procedure for the sparsification of \eqn{\mathbf{B}}. }
The matrices \code{zerosA}, \code{nonzerosA}, \code{zerosB} and \code{nonzerosB} comprise two columns, each row corresponding to an entry of either \eqn{\mathbf{A}} or \eqn{\mathbf{B}}. The first column contains the row indices and the second the column indices.
}
\details{
When \code{threshold = "localFDR"} the function, following Lutkepohl (2005), divides the elements of (possibly regularized) input matrix \eqn{\mathbf{A}} of autoregression coefficients by (approximations of) their standard errors. Subsequently, the support of the matrix \eqn{\mathbf{A}} is determined by usage of local FDR. In that case a mixture model is fitted to the nonredundant (standardized) elements of \eqn{\mathbf{A}} by \code{fdrtool}. The decision to retain elements is then based on the argument \code{FDRcut}. Elements with a posterior probability \eqn{>=q}  \code{FDRcut} (equalling 1 - local FDR) are retained. See Strimmer (2008) for further details. Alternatively, the support of \eqn{\mathbf{A}} is determined by simple thresholding on the absolute values of matrix entries (\code{threshold = "absValue"}). A third option (\code{threshold = "top"}) is to retain a prespecified number of matrix entries based on absolute values of the elements of \eqn{\mathbf{A}}. For example, one could wish to retain those
entries representing the ten strongest cross-temporal coefficients. 
     
The argument \code{absValueCut} is only used when \code{threshold = "absValue"}. The argument \code{FDRcut} is only used when \code{threshold = "localFDR"}. The argument \code{top} is only used when \code{threshold = "top"}. 

When prior to the sparsification knowledge on the support of \eqn{\mathbf{A}} and/or \eqn{\mathbf{B}} is specified through the options \code{zerosA} and/or \code{zerosB}, the corresponding elements of \eqn{\mathbf{A}} and  \eqn{\mathbf{B}} are then not taken along in the local FDR procedure.
}
\references{
Lutkepohl, H. (2005), \emph{New Introduction to Multiple Time Series Analysis}. Springer, Berlin.

Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', \emph{submitted}. 

Strimmer, K. (2008), ``fdrtool: a versatile R package for estimating local and tail area-based false discovery rates'', \emph{Bioinformatics} 24(12): 1461-1462.

Van Wieringen, W.N., Peeters, C.F.W. (2016), ``Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data'', \emph{Computational Statistics and Data Analysis}, 103, 284-303.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>, Carel F.W. Peeters.
}
\seealso{
\code{\link{ridgeVAR1}}, \code{\link[rags2ridges:sparsify]{sparsify}}, \code{\link[ragt2ridges:sparsifyVAR1]{sparsifyVAR1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
Ax     <- createA(p, "chain", nBands=1)

# generate time-varying covariates in accordance with VAR(1) process
X <- dataVAR1(n, T, Ax, SigmaE)

# set model parameters
B <- createA(p, "clique", nCliques=1)
A <- createA(p, "hub", nHubs=1)

# generate time-varying covariates in accordance with VAR(1) process
Y <- dataVARX1(X, A, B, SigmaE, lagX=0)

# fit VARX(1) model
VARX1hat <- ridgeVARX1(Y, X, 1, 1, 1, lagX=0)

## determine which elements of A are non-null
\dontrun{ Anullornot <- matrix(0, p, p) }
\dontrun{ Anullornot[sparsifyVARX1(X, VARX1hat$A, VARX1hat$B,
        solve(VARX1hat$P), threshold="localFDR")$nonzeros] <- 1 }
## REASON FOR NOT RUN:
## the employed local FDR approximation is only valid for reasonably sized 
## number of elements of A (say) at least p > 10 and, 
## consequently, a vector of 100 regression coefficients.

## plot non-null structure of A 
\dontrun{ edgeHeat(Anullornot) }
}

