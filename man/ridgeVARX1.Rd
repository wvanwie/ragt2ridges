\name{ridgeVARX1}
\alias{ridgeVARX1}
\title{
Ridge ML estimation of the VARX(1) model
}
\description{
Ridge penalized maximum likelihood estimation of the parameters of the first-order vector auto-regressive model with time-varying covariates, in shorthand VARX(1) model. The VARX(1) model explains the current vector of observations \eqn{\mathbf{Y}_{\ast,t+1}} by a linear combination of the previous observation endogeneous vector and an exogeneous time-varying covariate: \eqn{\mathbf{Y}_{\ast,t+1} = \mathbf{A} \mathbf{Y}_{\ast,t} + \mathbf{B} \mathbf{X}_{\ast,t+1} + \mathbf{\varepsilon}_{\ast,t+1}}, where \eqn{\mathbf{A}} and \eqn{\mathbf{B}} are the lag one autoregression and time-varying regression coefficient matrix, respectively, and \eqn{\mathbf{\varepsilon}_{\ast,t+2}} the vector of errors (or innovations). The VARX(1)-process is assumed to have mean zero. The experimental design is allowed to be unbalanced. 
}
\usage{
ridgeVARX1(Y, X, lambdaA=-1, lambdaB=-1, lambdaP=-1, lagX,
           targetA=matrix(0, dim(Y)[1], dim(Y)[1]), 
           targetB=matrix(0, dim(Y)[1], dim(X)[1]), 
           targetP=matrix(0, dim(Y)[1], dim(Y)[1]), targetPtype="none",
           fitAB="ml", zerosA=matrix(nrow=0, ncol=2), 
           zerosB=matrix(nrow=0, ncol=2), zerosAfit="sparse", 
           zerosBfit="sparse", zerosP=matrix(nrow=0, ncol=2), cliquesP=list(), 
           separatorsP=list(), unbalanced=matrix(nrow=0, ncol=2), diagP=FALSE, 
           efficient=TRUE, nInit=100, minSuccDiff=0.001)
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the response data. The first, second and third dimensions correspond to variates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{X}{ Three-dimensional \code{array} containing the time-varying covariate data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{lambdaA}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of \eqn{\mathbf{A}}, the matrix with autoregression coefficients. }
\item{lambdaB}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of \eqn{\mathbf{B}}, the matrix with regression coefficients of the time-varying covariates stored in \code{array} \code{X}. }
\item{lambdaP}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of the inverse error covariance matrix (\eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}): the precision matrix of the errors. }
\item{lagX}{ \code{Integer}, either \code{0} or {1}, specifying whether \eqn{\mathbf{X}_t} or \eqn{\mathbf{X}_{t-1}} affects \eqn{\mathbf{Y}_t}, respectively. }
\item{targetA}{ Target \code{matrix} to which the matrix \eqn{\mathbf{A}} is to be shrunken. }
\item{targetB}{ Target \code{matrix} to which the matrix \eqn{\mathbf{B}} is to be shrunken. }
\item{targetP}{ Target \code{matrix} to which the in the inverse error covariance matrix, the precision matrix, is to be shrunken. }
\item{fitAB}{ A \code{character}. If \code{fitAB="ml"} the parameters \eqn{\mathbf{A}} and \eqn{\mathbf{B}} are estimated by (penalized) maximum likelihood. If \code{fitAB="ss"} the parameters \eqn{\mathbf{A}}  and \eqn{\mathbf{B}} are estimated by (penalized) sum of squares. }
\item{targetPtype}{ A \code{character} indicating the type of target to be used for the precision matrix. When specified it overrules the \code{targetP}-option. See the \code{default.target}-function for the options. }
\item{zerosA}{ A \code{matrix} with indices of entries of \eqn{\mathbf{A}} that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{A}}. The first column contains the row indices and the second the column indices. }
\item{zerosB}{ A \code{matrix} with indices of entries of \eqn{\mathbf{B}} that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of \eqn{\mathbf{B}}. The first column contains the row indices and the second the column indices. }
\item{zerosAfit}{ A \code{character}, either "sparse" or "dense". With "sparse", the matrix \eqn{\mathbf{A}} is assumed to contain many zeros and a computational efficient implementation of its estimation is employed. If "dense", it is assumed that \eqn{\mathbf{A}} contains only few zeros and the estimation method is optimized computationally accordingly. }
\item{zerosBfit}{ A \code{character}, either "sparse" or "dense". With "sparse", the matrix \eqn{\mathbf{B}} is assumed to contain many zeros and a computational efficient implementation of its estimation is employed. If "dense", it is assumed that \eqn{\mathbf{B}} contain only few zeros and the estimation method is optimized computationally accordingly. }
\item{zerosP}{ A \code{matrix}-object with indices of entries of the precision matrix that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of the adjacency matrix. The first column contains the row indices and the second the column indices. The specified graph should be undirected and decomposable. If not, it is symmetrized and triangulated (unless \code{cliquesP} and \code{seperatorsP} are supplied). Hence, the employed zero structure may differ from the input \code{zerosP}. }
\item{cliquesP}{ A \code{list}-object containing the node indices per clique as object from the \code{\link[gRbase:rip]{rip}}-function. }
\item{separatorsP}{ A \code{list}-object containing the node indices per clique as object from the \code{\link[gRbase:rip]{rip}}-function. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{diagP}{ A \code{logical}, indicates whether the inverse error covariance matrix is assumed to be diagonal. }
\item{efficient}{ A \code{logical}, affects estimation of \eqn{\mathbf{A}}. Details below. }
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
\item{A}{Ridge ML estimate of the matrix \eqn{\mathbf{A}}, the \code{matrix} with lag one autoregression coefficients.}
\item{B}{Ridge ML estimate of the matrix \eqn{\mathbf{B}}, the \code{matrix} with regression coefficients of the time-varying covariates.}
\item{P}{Ridge ML estimate of the inverse error covariance \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
\item{lambdaA}{Positive \code{numeric} of length one: ridge penalty used in the estimation of \eqn{\mathbf{A}}.}
\item{lambdaB}{Positive \code{numeric} of length one: ridge penalty used in the estimation of \eqn{\mathbf{B}}.}
\item{lambdaP}{Positive \code{numeric} of length one: ridge penalty used in the estimation of inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
}
\note{
When the target of the precision matrix is specified through the \code{targetPtype}-argument, the target is data-driven and updated at each iteration when \code{fitAB="ml"}.

In case \eqn{\lambda_a \not= \lambda_b}, the ridge ML estimates (conditional on the current estimate of \eqn{\mathbf{\Omega}_{\varepsilon}}) of \eqn{\mathbf{A}} and \eqn{\mathbf{B}} are approximations. Their explicit evaluation involves a Kronecker product times a vector. Its full expansion is memory inefficient, in particular for (say) \eqn{p_y > 50}. If \eqn{\lambda_a = \lambda_b} the estimates of \eqn{\mathbf{A}} and \eqn{\mathbf{B}} can be evaluated through multiplications of matrices of order \eqn{p_y+p_x} (instead of \eqn{p_y^2 + p_y p_x}). The evaluation cannot be simplified computationally 
when \eqn{\lambda_a \not= \lambda_b}. To avoid the use of matrices of order \eqn{p_y^2 + p_y p_x} an approximation is used. For the approximation it is assumed that the contemporaneous covariance between \eqn{\mathbf{Y}_t} and \eqn{\mathbf{X}_t} is zero (for lag one), after which the estimates can be evaluated using matrices of order \eqn{p_y+p_x}.
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
Ax     <- createA(p, "chain", nBands=1)

# generate time-varying covariates in accordance with VAR(1) process
X <- dataVAR1(n, T, Ax, SigmaE)

# set model parameters
B <- createA(p, "clique", nCliques=1)
A <- createA(p, "hub", nHubs=1)

# generate time-varying covariates in accordance with VAR(1) process
Y <- dataVARX1(X, A, B, SigmaE, lagX=0)

# fit VARX(1) model
ridgeVARX1(Y, X, 1, 1, 1, lagX=0)
}

