\name{ridgeVAR1fused}
\alias{ridgeVAR1fused}
\title{
Fused ridge ML estimation of multiple VAR(1) model
}
\description{
Ridge penalized maximum likelihood estimation of the parameters of the first-order Vector Auto-Regressive model, with a (possibly) unbalanced experimental set-up. The VAR(1)-process is assumed to have mean zero.
}
\usage{
ridgeVAR1fused(Y, id, lambdaA=0, lambdaF=0, lambdaP=0, 
               targetA=matrix(0, dim(Y)[1], dim(Y)[1]), 
               targetP=matrix(0, dim(Y)[1], dim(Y)[1]), 
               targetPtype="none", fitA="ml", 
               zerosA=matrix(nrow=0, ncol=2), zerosAfit="sparse", 
               zerosP=matrix(nrow=0, ncol=2), cliquesP=list(), 
               separatorsP=list(), unbalanced=matrix(nrow=0, ncol=2), 
               diagP=FALSE, efficient=TRUE, nInit=100, nInitA=5, 
               minSuccDiff=0.001, minSuccDiffA=0.001)
} 
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{id}{ A vector with group indices comprising of integers only. First group is represented by '0', the next by '1', and so on until the last. }
\item{lambdaA}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of the \eqn{\mathbf{A}_g}, the matrices with autoregression coefficients. }
\item{lambdaF}{ Fused ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of \eqn{\mathbf{A}_g}, the matrices with autoregression coefficients. }
\item{lambdaP}{ Ridge penalty parameter (positive \code{numeric} of length 1) to be used in the estimation of the inverse error covariance matrix (\eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}): the precision matrix of the errors. }
\item{targetA}{ Target \code{matrix} to which the matrix \eqn{\mathbf{A}} is to be shrunken. This target is shared among the groups (otherwise why fuse?). }
\item{targetP}{ Target \code{matrix} to which the in the inverse error covariance matrix, the precision matrix, is to be shrunken. }
\item{fitA}{ A \code{character}. If \code{fitA="ml"} the parameters \eqn{\mathbf{A}_g} ared estimated by (penalized) maximum likelihood. If \code{fitA="ss"} the parameter \eqn{\mathbf{A}} is estimate by (penalized) sum of squares. The latter is much faster. }
\item{targetPtype}{ A \code{character} indicating the type of target to be used for the precision matrix. When specified it overrules the \code{targetP}-option. See the \code{default.target}-function for the options. }
\item{zerosA}{ A \code{matrix} with indices of entries of \eqn{\mathbf{A}_g}'s that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of the \eqn{\mathbf{A}_s}. The first column contains the row indices and the second the column indices. The support is shared among the groups (otherwise why fuse?). }
\item{zerosAfit}{ A \code{character}, either "sparse" or "dense". With "sparse", the matrix \eqn{\mathbf{A}} is assumed to contain many zeros and a computational efficient implementation of its estimation is employed. If "dense", it is assumed that \eqn{\mathbf{A}_g}'s contain only few zeros and the estimation method is optimized computationally accordingly.  }
\item{zerosP}{ A \code{matrix}-object with indices of entries of the precision matrix that are constrained to zero. The matrix comprises two columns, each row corresponding to an entry of the adjacency matrix. The first column contains the row indices and the second the column indices. The specified graph should be undirected and decomposable. If not, it is symmetrized and triangulated (unless \code{cliquesP} and \code{seperatorsP} are supplied). Hence, the employed zero structure may differ from the input \code{zerosP}. }
\item{cliquesP}{ A \code{list}-object containing the node indices per clique as object from the \code{rip}-function. }
\item{separatorsP}{ A \code{list}-object containing the node indices per clique as object from the \code{rip}-function. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{diagP}{ A \code{logical}, indicates whether the inverse error covariance matrix is assumed to be diagonal. }
\item{efficient}{ A \code{logical}, affects estimation of the \eqn{\mathbf{A}_g}. Details below. }
\item{nInit}{ Maximum number of iterations (positive \code{numeric} of length 1) to be used in maximum likelihood estimation. }
\item{nInitA}{ Maximum number of iterations (positive \code{numeric} of length 1) to be used in fused estimation of the autoregression matrices \eqn{\mathbf{A}_g}, given the current estimate of \eqn{\mathbf{\Omega}_{\varepsilon}}. }
\item{minSuccDiff}{ Minimum distance (positive \code{numeric} of length 1) between estimates of two successive iterations to be achieved. }
\item{minSuccDiffA}{ Minimum distance (positive \code{numeric} of length 1) between the \eqn{\mathbf{A}_g} estimates of two successive fused estimation iterations to be achieved. }
}
\details{
If \code{diagP=TRUE}, no penalization to estimation of the covariance matrix is applied. Consequently, the arguments \code{lambdaP} and \code{targetP} are ignored (if supplied).

The ridge ML estimator employs the following estimator of the variance of the VAR(1) process:
\deqn{ \frac{1}{n (\mathcal{T} - 1)} \sum_{i=1}^{n} \sum_{t=2}^{\mathcal{T}} \mathbf{Y}_{\ast,i,t} \mathbf{Y}_{\ast,i,t}^{\mathrm{T}}. }
This is used when \code{efficient=FALSE}. However, a more efficient estimator of this variance can be used
\deqn{ \frac{1}{n \mathcal{T}} \sum_{i=1}^{n} \sum_{t=1}^{\mathcal{T}} \mathbf{Y}_{\ast,i,t} \mathbf{Y}_{\ast,i,t}^{\mathrm{T}}, }
which is achieved by setting when \code{efficient=TRUE}. Both estimators are adjusted accordingly when dealing with an unbalanced design. 
}
\value{
A list-object with slots:
\item{As}{Ridge ML estimates of the matrices \eqn{\mathbf{A}_g}, stacked and stored as a single rectangular \code{matrices}.}
\item{P}{Ridge ML estimate of the inverse error covariance \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
\item{lambdaA}{Positive \code{numeric} of length one: ridge penalty used in the estimation of the \eqn{\mathbf{A}_g}.}
\item{lambdaF}{Positive \code{numeric} of length one: fused ridge penalty used in the estimation of the \eqn{\mathbf{A}_g}.}
\item{lambdaP}{Positive \code{numeric} of length one: ridge penalty used in the estimation of inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
}
\note{
When the target of the precision matrix is specified through the \code{targetPtype}-argument, the target is data-driven and updated at each iteration when \code{fitA="ml"}. 
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link{loglikLOOCVVAR1}}, \code{\link[ragt2ridges:ridgeVAR1]{ridgeVAR1}}, \code{\link[rags2ridges:ridgePchordal]{ridgePchordal}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points, G=groups)
p <- 3; n <- 12; T <- 10; G <- 3

# set model parameters
SigmaE <- matrix(1/2, p, p)
diag(SigmaE) <- 1
A1 <- -createA(p, "clique", nCliques=1, nonzeroA=0.1)
A2 <- t(createA(p, "chain", nBands=1, nonzeroA=0.1))
A3 <- (A1 + A2) / 2

# generate data
Y1 <- dataVAR1(n/G, T, A1, SigmaE)
Y2 <- dataVAR1(n/G, T, A2, SigmaE)
Y3 <- dataVAR1(n/G, T, A3, SigmaE)
Y  <- abind::abind(Y1, Y2, Y3, along=3)
id <- c(rep(1, n/G), rep(2, n/G), rep(3, n/G))-1

VAR1hats <- ridgeVAR1fused(Y, id, lambdaA=1, lambdaF=1, lambdaP=1)
}

