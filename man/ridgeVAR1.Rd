\name{ridgeVAR1}
\alias{ridgeVAR1}
\title{
Ridge ML estimation of the VAR(1) model
}
\description{
Ridge penalized maximum likelihood estimation of the parameters of the VAR(1), first-order vector auto-regressive, model. The VAR(1) model explains the current vector of observations \eqn{\mathbf{Y}_{\ast,t+1}} by a linear combination of the previous observation vector: \eqn{\mathbf{Y}_{\ast,t+1} = \mathbf{A} \mathbf{Y}_{\ast,t} + \mathbf{\varepsilon}_{\ast,t+1}}, where \eqn{\mathbf{A}} is the autoregression coefficient matrix and \eqn{\mathbf{\varepsilon}_{\ast,t+1}} the vector of errors (or innovations). The VAR(1)-process is assumed to have mean zero. The experimental design is allowed to be unbalanced. 
}
\usage{
ridgeVAR1(Y, lambdaA=0, lambdaP=0, 
          targetA=matrix(0, dim(Y)[1], dim(Y)[1]), 
          targetP=matrix(0, dim(Y)[1], dim(Y)[1]), targetPtype="none",
          fitA="ml", zerosA=matrix(nrow=0, ncol=2), 
          zerosAfit="sparse", zerosP=matrix(nrow=0, ncol=2), 
          cliquesP=list(), separatorsP=list(), 
          unbalanced=matrix(nrow=0, ncol=2), diagP=FALSE, 
          efficient=TRUE, nInit=100, minSuccDiff=0.001)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{lambdaA}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of \eqn{\mathbf{A}}, the matrix with autoregression coefficients. }
\item{lambdaP}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of the inverse error covariance matrix (\eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}): the precision matrix of the errors. }
\item{targetA}{ Target \code{matrix} to which the matrix \eqn{\mathbf{A}} is to be shrunken. }
\item{targetP}{ Target \code{matrix} to which the in the inverse error covariance matrix, the precision matrix, is to be shrunken. }
\item{fitA}{ A \code{character}. If \code{fitA="ml"} the parameter \eqn{\mathbf{A}} is estimate by (penalized) maximum likelihood. If \code{fitA="ss"} the parameter \eqn{\mathbf{A}} is estimate by (penalized) sum of squares. The latter being much faster as it need not iterate. }
\item{targetPtype}{ A \code{character} indicating the type of target to be used for the precision matrix. When specified it overrules the \code{targetP}-option. See the \code{default.target}-function for the options. }
\item{zerosA}{ A \code{matrix} with indices of entries of \eqn{\mathbf{A}} that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{A}}. The first column contains the row indices and the second the column indices. }
\item{zerosAfit}{ A \code{character}, either "sparse" or "dense". With "sparse", the matrix \eqn{\mathbf{A}} is assumed to contain many zeros and a computational efficient implementation of its estimation is employed. If "dense", it is assumed that \eqn{\mathbf{A}} contains only few zeros and the estimation method is optimized computationally accordingly.  }
\item{zerosP}{ A \code{matrix}-object with indices of entries of the precision matrix that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of the adjacency matrix. The first column contains the row indices and the second the column indices. The specified graph should be undirected and decomposable. If not, it is symmetrized and triangulated (unless \code{cliquesP} and \code{seperatorsP} are supplied). Hence, the employed zero structure may differ from the input \code{zerosP}. }
\item{cliquesP}{ A \code{list}-object containing the node indices per clique as object from the \code{rip}-function. }
\item{separatorsP}{ A \code{list}-object containing the node indices per clique as object from the \code{rip}-function. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{diagP}{ A \code{logical}, indicates whether the inverse error covariance matrix is assumed to be diagonal. }
\item{efficient}{ A \code{logical}, affects estimation of \eqn{\mathbf{A}}. Details below. }
\item{nInit}{ Maximum number of iterations (positive \code{numeric} of length 1) to be used in maximum likelihood estimation. }
\item{minSuccDiff}{ Minimum distance (positive \code{numeric} of length 1) between estimates of two successive iterations to be achieved. }
}
\details{
The ridge ML estimator employs the following estimator of the variance of the VAR(1) process:
\deqn{ \frac{1}{n (\mathcal{T} - 1)} \sum_{i=1}^{n} \sum_{t=2}^{\mathcal{T}} \mathbf{Y}_{\ast,i,t} \mathbf{Y}_{\ast,i,t}^{\mathrm{T}}. }
This is used when \code{efficient=FALSE}. However, a more efficient estimator of this variance can be used
\deqn{ \frac{1}{n \mathcal{T}} \sum_{i=1}^{n} \sum_{t=1}^{\mathcal{T}} \mathbf{Y}_{\ast,i,t} \mathbf{Y}_{\ast,i,t}^{\mathrm{T}}, }
which is achieved by setting when \code{efficient=TRUE}. Both estimators are adjusted accordingly when dealing with an unbalanced design. 
}
\value{
A list-object with slots:
\item{A}{Ridge ML estimate of the matrix \eqn{\mathbf{A}}, the \code{matrix} with lag one auto-regressive coefficients.}
\item{P}{Ridge ML estimate of the inverse error covariance \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
\item{lambdaA}{Positive \code{numeric} of length one: ridge penalty used in the estimation of \eqn{\mathbf{A}}.}
\item{lambdaP}{Positive \code{numeric} of length one: ridge penalty used in the estimation of inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
}
\note{
When the target of the precision matrix is specified through the \code{targetPtype}-argument, the target is data-driven (for both \code{fitA="ss"} and \code{fitA="ml"}). In particular, it is updated at each iteration when \code{fitA="ml"}. 
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link{loglikLOOCVVAR1}}, \code{\link[rags2ridges:ridgeP]{ridgeP}}, \code{\link[rags2ridges:default.target]{default.target}}, \code{\link[rags2ridges:ridgePchordal]{ridgePchordal}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A      <- createA(p, "chain")

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

# fit VAR(1) model
ridgeVAR1(Y, 1, 1)$A
}

