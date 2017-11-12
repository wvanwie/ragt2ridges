\name{loglikLOOCVcontourVARX1}
\alias{loglikLOOCVcontourVARX1}
\title{
Contourplot of the LOOCV log-likelihood of VARX(1) model
}
\description{
Evaluates the leave-one-out cross-validated log-likelihood of the VARX(1) model over a grid of the ridge penalty parameters (\eqn{\lambda_a} and \eqn{\lambda_b}) for the autoregression and time-varying covariate regression coefficient matrices \eqn{\mathbf{A}} and \eqn{\mathbf{B}}, respectively, while keeping \eqn{\lambda_{\omega}}, the penalty parameter of the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}, fixed at a user-specified value. The result is plotted as a contour plot, which facilitates the choice of optimal penalty parameters. The function also works with a (possibly) unbalanced experimental set-up. The VARX(1)-process is assumed to have mean zero.
}
\usage{
loglikLOOCVcontourVARX1(lambdaAgrid, lambdaBgrid, lambdaPgrid, Y, X,
                        lagX=0, figure=TRUE, verbose=TRUE, ...)
}
\arguments{
\item{lambdaAgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to the \eqn{\lambda_a} (the penalty parameter for the lag one autoregression coefficient matrix \eqn{\mathbf{A}}). }
\item{lambdaBgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to the \eqn{\lambda_b} (the penalty parameters for the regression coefficient matrix \eqn{\mathbf{B}} of the time-varying covariates. }
\item{lambdaPgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to the \eqn{\lambda_{\omega}} (the penalty parameters for the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}). }
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered variate-wise. }
\item{X}{ Three-dimensional \code{array} containing the time-varying covariates. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{lagX}{ An \code{integer}, either \code{0} or {1}, specifying whether \eqn{\mathbf{X}_t} or \eqn{\mathbf{X}_{t-1}} affects \eqn{\mathbf{Y}_t}, respectively. }
\item{figure}{ A \code{logical} indicating whether the contour plot should be generated. }
\item{verbose}{ A \code{logical} indicator: should intermediate output be printed on the screen? }
\item{...}{  Other arguments to be passed on (indirectly) to \code{\link{ridgeVARX1}}. }
}
\value{
A \code{list}-object with slots:
\item{lambdaA}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_a} (the penalty parameter for the autoregression coefficient matrix \eqn{\mathbf{A}}). }
\item{lambdaB}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_b} (the penalty parameter for time-varying covariate regression coefficient matrix \eqn{\mathbf{B}}. }
\item{lambdaP}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_{\omega}} (the penalty parameter for the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}.}
\item{llLOOCV}{A \code{matrix} of leave-one-out cross-validated log-likelihoods. Rows and columns correspond to \eqn{\lambda_a} and \eqn{\lambda_{b}} values, respectively. }
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of network models from time-course omics data", submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>, Viktorian Miok.
}
\note{
When \code{lambdaAgrid}, \code{lambdaBgrid} and \code{lambdaPgrid} are all vectors of length exceeding one, the contour is determined for the grid formed by the Cartesius product of \code{lambdaAgrid} and \code{lambdaBgrid} meanwhile restricting \code{lambdaPgrid} to its first element.

Internally, this function calls the \code{\link{loglikLOOCVVARX1}}-function, which evaluates the minus (!) LOOCV log-likelihood (for practical reasons). For interpretation purposes
\code{\link{loglikLOOCVcontourVARX1}} provides the regular LOOCV log-likelihood (that is, without the minus).
}
\seealso{
\code{\link{loglikLOOCVcontourVAR1}}, \code{\link{loglikLOOCVcontourVAR1fused}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
Ax     <- createA(3, "chain")

# generate time-varying covariate data
X <- dataVAR1(n, T, Ax, SigmaE)

# regression parameter matrices of VARX(1) model
A <- createA(p, topology="clique", nonzeroA=0.1, nClique=1)
B <- createA(p, topology="hub", nonzeroA=0.1, nHubs=1)

# generate data
Y <- dataVARX1(X, A, B, SigmaE, lagX=0)

## plot contour of cross-validated likelihood
\dontrun{ lambdaAgrid <- seq(0.01, 1, length.out=20) }
\dontrun{ lambdaBgrid <- seq(0.01, 1000, length.out=20) }
\dontrun{ lambdaPgrid <- seq(0.01, 1000, length.out=20) }
\dontrun{ loglikLOOCVcontourVARX1(lambdaAgrid, lambdaBgrid, lambdaPgrid, Y, X) }

## determine optimal values of the penalty parameters
\dontrun{optLambdas <- constrOptim(c(1,1,1), loglikLOOCVVARX1, gr=NULL, }
\dontrun{              ui=diag(3), ci=c(0,0,0), Y=Y, X=X,}
\dontrun{              control=list(reltol=0.01))$par }

## add point of optimum
\dontrun{ points(optLambdas[1], optLambdas[2], pch=20, cex=2, }
\dontrun{ col="red") }
}

