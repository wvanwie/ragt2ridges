\name{loglikLOOCVcontourVAR2}
\alias{loglikLOOCVcontourVAR2}
\title{
Contourplot of LOOCV log-likelihood of the VAR(2) model
}
\description{
Evaluates the leave-one-out cross-validated log-likelihood of the estimated VAR(2) model over a grid of the ridge penalty parameters (\eqn{\lambda_{a1}} and \eqn{\lambda_{a2}}) for the lag one and lag two autoregression coefficient matrices \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2} (respectively), while keeping \eqn{\lambda_{\omega}}, the penalty parameter of the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}, fixed at a user-specified value. The result is plotted as a contour plot, which facilitates the choice of optimal penalty parameters. The function also works with a (possibly) unbalanced experimental set-up. The VAR(2)-process are assumed to have mean zero.
}
\usage{
loglikLOOCVcontourVAR2(lambdaA1grid, lambdaA2grid, 
                       lambdaPgrid, Y, figure=TRUE, verbose=TRUE, ...)
}
\arguments{
\item{lambdaA1grid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to \eqn{\lambda_{a1}}, the ridge penalty parameter for the lag one autoregression coefficient matrix \eqn{\mathbf{A}_1}. }
\item{lambdaA2grid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to \eqn{\lambda_{a2}}, the ridge penalty parameter for the lag two autoregression coefficient matrix \eqn{\mathbf{A}_2}. }
\item{lambdaPgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to the \eqn{\lambda_{\omega}} (the penalty parameters for the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}). }
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{figure}{ A \code{logical}, indicating whether the contour plot should be generated. }
\item{verbose}{ A \code{logical} indicator: should intermediate output be printed on the screen? }
\item{...}{  Other arguments to be passed on (indirectly) to \code{\link{ridgeVAR2}}. }
}
\value{
A \code{list}-object with slots:
\item{lambdaA1}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_{a1}} (the ridge penalty parameter for the lag one autoregression coefficient matrix \eqn{\mathbf{A}_1}).}
\item{lambdaA2}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_{a2}} (the ridge penalty parameter for the lag one autoregression coefficient matrix \eqn{\mathbf{A}_2}).}
\item{lambdaP}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_{\omega}} (the penalty parameter for the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}).}
\item{llLOOCV}{A \code{matrix} of leave-one-out cross-validated log-likelihoods. Rows and columns correspond to \eqn{\lambda_{a1}} and \eqn{\lambda_{a2}} values, respectively.}
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>, Viktorian Miok.
}
\note{
When \code{lambdaA1grid}, \code{lambdaA2grid} and \code{lambdaPgrid} are all vectors of length exceeding one, the contour is determined for the grid formed by the Cartesius product of \code{lambdaA1grid} and \code{lambdaA2grid} meanwhile restricting \code{lambdaPgrid} to its first element.

Internally, this function calls the \code{\link{loglikLOOCVVAR2}}-function, which evaluates the minus (!) LOOCV log-likelihood (for practical reasons). For interpretation purposes
\code{\link{loglikLOOCVcontourVAR2}} provides the regular LOOCV log-likelihood (that is, without the minus).
}
\seealso{
\code{\link{loglikLOOCVcontourVAR1}}, \code{\link{loglikLOOCVcontourVARX1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A1 <- -createA(p, "clique", nCliques=1, nonzeroA=0.1)
A2 <- t(createA(p, "chain", nBands=1,   nonzeroA=0.1))

# generate data
Y <- dataVAR2(n, T, A1, A2, SigmaE)

## plot contour of cross-validated likelihood
\dontrun{ lambdaA1grid <- seq(0.01, 1, length.out=20) }
\dontrun{ lambdaA2grid <- seq(0.01, 1000, length.out=20) }
\dontrun{ lambdaPgrid  <- seq(0.01, 1000, length.out=20) }
\dontrun{ loglikLOOCVcontourVAR2(lambdaA1grid, lambdaA2grid, lambdaPgrid, Y) }

## determine optimal values of the penalty parameters
\dontrun{optLambdas <- constrOptim(c(1,1,1), loglikLOOCVVAR2, gr=NULL, }
\dontrun{              ui=diag(3), ci=c(0,0,0), Y=Y, }
\dontrun{              control=list(reltol=0.01))$par }

## add point of optimum
\dontrun{ points(optLambdas[1], optLambdas[2], pch=20, cex=2, }
\dontrun{ col="red") }
}

