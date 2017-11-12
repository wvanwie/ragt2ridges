\name{ridgeVAR2}
\alias{ridgeVAR2}
\title{
Ridge ML estimation of the VAR(2) model
}
\description{
Ridge penalized maximum likelihood estimation of the parameters of the VAR(2), second-order vector auto-regressive, model. The VAR(2) model explains the current vector of observations \eqn{\mathbf{Y}_{\ast,t+2}} by a linear combination of the previous two observation vectors: \eqn{\mathbf{Y}_{\ast,t+1} = \mathbf{A}_1 \mathbf{Y}_{\ast,t+1} + \mathbf{A}_2 \mathbf{Y}_{\ast,t} + \mathbf{\varepsilon}_{\ast,t+2}}, where \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2} are the lag one and two autoregression coefficient matrices, respectively, and \eqn{\mathbf{\varepsilon}_{\ast,t+2}} the vector of errors (or innovations). The VAR(2)-process is assumed to have mean zero. The experimental design is allowed to be unbalanced. 
}
\usage{
ridgeVAR2(Y, lambdaA1=-1, lambdaA2=-1, lambdaP=-1, 
          targetA1=matrix(0, dim(Y)[1], dim(Y)[1]), 
          targetA2=matrix(0, dim(Y)[1], dim(Y)[1]), 
          targetP=matrix(0, dim(Y)[1], dim(Y)[1]), 
          targetPtype="none", fitA12="ml", 
          zerosA1=matrix(nrow=0, ncol=2), 
          zerosA2=matrix(nrow=0, ncol=2), zerosA1fit="sparse", 
          zerosA2fit="sparse", zerosP=matrix(nrow=0, ncol=2), 
          cliquesP=list(), separatorsP=list(), 
          unbalanced=matrix(nrow=0, ncol=2), diagP=FALSE, 
          efficient=TRUE, nInit=100, minSuccDiff=0.001)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the response data. The first, second and third dimensions correspond to variates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{lambdaA1}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of \eqn{\mathbf{A}_1}, the matrix with the lag one autoregression coefficients. }
\item{lambdaA2}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of \eqn{\mathbf{A}_2}, the matrix with the lag two autoregression coefficients. }
\item{lambdaP}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of the inverse error covariance matrix (\eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}): the precision matrix of the errors. }
\item{targetA1}{ Target \code{matrix} to which the matrix \eqn{\mathbf{A}_1} is to be shrunken. }
\item{targetA2}{ Target \code{matrix} to which the matrix \eqn{\mathbf{A}_2} is to be shrunken. }
\item{targetP}{ Target \code{matrix} to which the in the inverse error covariance matrix, the precision matrix, is to be shrunken. }
\item{fitA12}{ A \code{character}. If \code{fitAB="ml"} the parameters \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2} are estimated by (penalized) maximum likelihood. If \code{fitAB="ss"} the parameters \eqn{\mathbf{A}_1}  and \eqn{\mathbf{A}_2} are estimated by (penalized) sum of squares. }
\item{targetPtype}{ A \code{character} indicating the type of target to be used for the precision matrix. When specified it overrules the \code{targetP}-option. See the \code{default.target}-function for the options. }
\item{zerosA1}{ A \code{matrix} with indices of entries of \eqn{\mathbf{A}_1} that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{A}_1}. The first column contains the row indices and the second the column indices. }
\item{zerosA2}{ A \code{matrix} with indices of entries of \eqn{\mathbf{A}_2} that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{A}_2}. The first column contains the row indices and the second the column indices. }
\item{zerosA1fit}{ A \code{character}, either "sparse" or "dense". With "sparse", the matrix \eqn{\mathbf{A}_1} is assumed to contain many zeros and a computational efficient implementation of its estimation is employed. If "dense", it is assumed that \eqn{\mathbf{A}_1} contains only few zeros and the estimation method is optimized computationally accordingly. }
\item{zerosA2fit}{ A \code{character}, either "sparse" or "dense". With "sparse", the matrix \eqn{\mathbf{A}_2} is assumed to contain many zeros and a computational efficient implementation of its estimation is employed. If "dense", it is assumed that \eqn{\mathbf{A}_2} contain only few zeros and the estimation method is optimized computationally accordingly. }
\item{zerosP}{ A \code{matrix}-object with indices of entries of the precision matrix that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of the adjacency matrix. The first column contains the row indices and the second the column indices. The specified graph should be undirected and decomposable. If not, it is symmetrized and triangulated (unless \code{cliquesP} and \code{seperatorsP} are supplied). Hence, the employed zero structure may differ from the input \code{zerosP}. }
\item{cliquesP}{ A \code{list}-object containing the node indices per clique as object from the \code{\link[gRbase:rip]{rip}}-function. }
\item{separatorsP}{ A \code{list}-object containing the node indices per clique as object from the \code{\link[gRbase:rip]{rip}}-function. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{diagP}{ \code{Logical}, indicates whether the inverse error covariance matrix is assumed to be diagonal. }
\item{efficient}{ \code{Logical}, affects estimation of \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2} directly. Details below. }
\item{nInit}{ Maximum number of iterations (positive \code{numeric} of length 1) to be used in maximum likelihood estimation. }
\item{minSuccDiff}{ Minimum distance (positive \code{numeric} of length 1) between estimates of two successive iterations to be achieved. }
}
\details{
If \code{diagP=TRUE}, no penalization to estimation of the covariance matrix is applied. Consequently, the arguments \code{lambdaP} and \code{targetP} are ignored (if supplied).

The ridge ML estimator employs the following estimator of the variance of the VARX(1) process:
\deqn{ \frac{1}{n (\mathcal{T} - 1)} \sum_{i=1}^{n} \sum_{t=2}^{\mathcal{T}} \mathbf{Y}_{\ast,i,t} \mathbf{Y}_{\ast,i,t}^{\top}. }
This is used when \code{efficient=FALSE}. However, a more efficient estimator of this variance can be used
\deqn{ \frac{1}{n \mathcal{T}} \sum_{i=1}^{n} \sum_{t=1}^{\mathcal{T}} \mathbf{Y}_{\ast,i,t} \mathbf{Y}_{\ast,i,t}^{\top}, }
which is achieved by setting when \code{efficient=TRUE}. Both estimators are adjusted accordingly when dealing with an unbalanced design. 
}
\value{
A list-object with slots:
\item{A1}{Ridge ML estimate of the matrix \eqn{\mathbf{A}_1}, the \code{matrix} with lag one auto-regressive coefficients.}
\item{A2}{Ridge ML estimate of the matrix \eqn{\mathbf{A}_2}, the \code{matrix} with lag two auto-regressive coefficients.}
\item{P}{Ridge ML estimate of the inverse error covariance \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
\item{lambdaA1}{Positive \code{numeric} of length one: ridge penalty used in the estimation of \eqn{\mathbf{A}_1}.}
\item{lambdaA2}{Positive \code{numeric} of length one: ridge penalty used in the estimation of \eqn{\mathbf{A}_2}.}
\item{lambdaP}{Positive \code{numeric} of length one: ridge penalty used in the estimation of inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
}
\note{
When the target of the precision matrix is specified through the \code{targetPtype}-argument, the target is data-driven and updated at each iteration when \code{fitAB="ml"}.

In case \eqn{\lambda_{a1} \not= \lambda_{a2}}, the ridge ML estimates (conditional on the current estimate of \eqn{\mathbf{\Omega}_{\varepsilon}}) of \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2} are approximations. Their explicit evaluation involves a Kronecker product times a vector. Its full expansion is memory inefficient, in particular for (say) \eqn{p_y > 50}. If \eqn{\lambda_{a1} = \lambda_{a2}} the estimates of \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2} can be evaluated through multiplications of matrices of dimensions \eqn{2 p_y \times 2 p_y} (instead of \eqn{2 p_y^2 \times  2 p_y^2}). The evaluation cannot be simplified computationally when \eqn{\lambda_{a1} \not= \lambda_{a1}}. To avoid the use of prohitedly large matrices an approximation is used. For the approximation it is assumed that the contemporaneous covariance between \eqn{\mathbf{Y}_t} and \eqn{\mathbf{X}_t} is zero (for lag one), after which the estimates can be evaluated using matrices of order \eqn{p_y+p_x}.
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link[rags2ridges:default.target]{default.target}}, \code{\link{loglikLOOCVVAR1}}, \code{\link[rags2ridges:ridgeP]{ridgeP}}, \code{\link[rags2ridges:ridgePchordal]{ridgePchordal}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A1     <- createA(p, "clique", nCliques=1)
A2     <- createA(p, "hub", nHubs=1)

# generate time-varying covariates in accordance with VAR(2) process
Y <- dataVAR2(n, T, A1, A2, SigmaE)

# fit VAR(2) model
ridgeVAR2(Y, 1, 1, 1)
}

