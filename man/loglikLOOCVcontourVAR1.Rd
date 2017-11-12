\name{loglikLOOCVcontourVAR1}
\alias{loglikLOOCVcontourVAR1}
\title{
Contourplot of LOOCV log-likelihood of VAR(1) model
}
\description{
Evaluates the leave-one-out cross-validated log-likelihood of the VAR(1) model for a given grid of the ridge penalty parameters (\eqn{\lambda_a} and \eqn{\lambda_{\omega}}) for the autoregression coefficient matrix \eqn{\mathbf{A}} and the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}, respectively). The result is plotted as a contour plot, which facilitates the choice of optimal penalty parameters. The function also works with a (possibly) unbalanced experimental set-up. The VAR(1)-process is assumed to have mean zero.
}
\usage{
loglikLOOCVcontourVAR1(lambdaAgrid, lambdaPgrid, Y, figure=TRUE, 
                       verbose=TRUE, ...)
}
\arguments{
\item{lambdaAgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to the \eqn{\lambda_a} (the penalty parameter for the autoregression coefficient matrix \eqn{\mathbf{A}}). }
\item{lambdaPgrid}{ A \code{numeric} of length larger than one, comprising positive numbers only. It contains the grid points corresponding to the \eqn{\lambda_{\omega}} (the penalty parameters for the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}). }
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to centered covariate-wise. }
\item{figure}{ A \code{logical}, indicating whether the contour plot should be generated. }
\item{verbose}{ A \code{logical} indicator: should intermediate output be printed on the screen? }
\item{...}{  Other arguments to be passed on (indirectly) to \code{\link{ridgeVAR1}}. }
}
\value{
A \code{list}-object with slots:
\item{lambdaA}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_a} (the penalty parameter for the autoregression coefficient matrix \eqn{\mathbf{A}}).}
\item{lambdaP}{A \code{numeric} with the grid points corresponding to \eqn{\lambda_{\omega}} (the penalty parameter for the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}).}
\item{llLOOCV}{A \code{matrix} of leave-one-out cross-validated log-likelihoods. Rows and columns correspond to \eqn{\lambda_a} and \eqn{\lambda_{\omega}} values, respectively.}
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\note{
Internally, this function calls the \code{\link{loglikLOOCVVAR1}}-function, which evaluates the minus (!) LOOCV log-likelihood (for practical reasons). For interpretation purposes
\code{\link{loglikLOOCVcontourVAR1}} provides the regular LOOCV log-likelihood (that is, without the minus).
}
\seealso{
\code{\link{loglikLOOCVVAR1}}.
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

