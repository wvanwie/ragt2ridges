\name{loglikLOOCVcontourVAR1fused}
\alias{loglikLOOCVcontourVAR1fused}
\title{
Contourplot of LOOCV log-likelihood of multiple VAR(1) models
}
\description{
Evaluates the leave-one-out cross-validated log-likelihood of multiple jointly estimated VAR(1) models over a grid of the (fused) ridge penalty parameters (\eqn{\lambda_a} and \eqn{\lambda_{f}}) for the autoregression coefficient matrices \eqn{\mathbf{A}_g}, while keeping \eqn{\lambda_{\omega}}, the penalty parameter of the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}, fixed at a user-specified value. The result is plotted as a contour plot, which facilitates the choice of optimal penalty parameters. The function also works with a (possibly) unbalanced experimental set-up. The VAR(1)-processes are assumed to have mean zero.
}
\usage{
loglikLOOCVcontourVAR1fused(lambdaAgrid, lambdaFgrid, Y, id, 
                            lambdaP, figure=TRUE, verbose=TRUE, ...)
}
\arguments{
\item{lambdaAgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to \eqn{\lambda_a}, the ridge penalty parameter for the autoregression coefficient matrices \eqn{\mathbf{A}_g}. }
\item{lambdaFgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to \eqn{\lambda_{f}}, the fused ridge penalty parameter for the autoregression coefficient matrices \eqn{\mathbf{A}_g}. }
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{id}{ A vector with group indices comprising of integers only. First group is represented by '0', the next by '1', and so on until the last. }
\item{lambdaP}{ A \code{numeric} specifying the value at which \eqn{\lambda_{\omega}} should be fixed. } 
\item{figure}{ A \code{logical}, indicating whether the contour plot should be generated. }
\item{verbose}{ A \code{logical} indicator: should intermediate output be printed on the screen? }
\item{...}{  Other arguments to be passed on (indirectly) to \code{\link{ridgeVAR1fused}}. }
}
\value{
A \code{list}-object with slots:
\item{lambdaA}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_a} (the ridge penalty parameter for the autoregression coefficient matrices \eqn{\mathbf{A}_g}).}
\item{lambdaF}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_{f}} (the fused ridge penalty parameter for the autoregression coefficient matrices \eqn{\mathbf{A}_g}. }
\item{llLOOCV}{A \code{matrix} of leave-one-out cross-validated log-likelihoods. Rows and columns correspond to \eqn{\lambda_a} and \eqn{\lambda_{\omega}} values, respectively.}
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of network models from time-course omics data", \emph{submitted}. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>, Viktorian Miok.

}
\note{
Internally, this function calls the \code{\link{loglikLOOCVVAR1}}-function, which evaluates the minus (!) LOOCV log-likelihood (for practical reasons). For interpretation purposes
\code{\link{loglikLOOCVcontourVAR1}} provides the regular LOOCV log-likelihood (that is, without the minus).
}
\seealso{
\code{\link{loglikLOOCVcontourVAR1}}, \code{\link{loglikLOOCVcontourVARX1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A <- createA(p, "chain")

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

## plot contour of cross-validated likelihood
\dontrun{ lambdaAgrid <- seq(0.01, 1, length.out=20) }
\dontrun{ lambdaPgrid <- seq(0.01, 1000, length.out=20) }
\dontrun{ loglikLOOCVcontourVAR1(lambdaAgrid, lambdaPgrid, Y) }

## determine optimal values of the penalty parameters
\dontrun{optLambdas <- constrOptim(c(1,1), loglikLOOCVVAR1, gr=NULL, }
\dontrun{              ui=diag(2), ci=c(0,0), Y=Y, }
\dontrun{              control=list(reltol=0.01))$par }

## add point of optimum
\dontrun{ points(optLambdas[1], optLambdas[2], pch=20, cex=2, }
\dontrun{ col="red") }
}

