\name{loglikLOOCVVAR1}
\alias{loglikLOOCVVAR1}
\title{
Leave-one-out (minus) cross-validated log-likelihood of VAR(1) model
}
\description{
Evaluation of the (minus) leave-one-out cross-validated log-likelihood of the VAR(1) model for given choices of the ridge penalty parameters (\eqn{\lambda_a} and \eqn{\lambda_{\omega}} for the autoregression coefficient matrix \eqn{\mathbf{A}} and the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}, respectively). The functions also works with a (possibly) unbalanced experimental set-up. The VAR(1)-process is assumed to have mean zero.
}
\usage{
loglikLOOCVVAR1(lambdas, Y, unbalanced=matrix(nrow=0, ncol=2), ...)
}
\arguments{
\item{lambdas}{ A \code{numeric} of length two, comprising positive values only. It contains the ridge penalty parameters to be used in the estimation of \eqn{\mathbf{A}} and the precision matrix of the errors, respectively. }
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{...}{  Other arguments to be passed to \code{ridgeVAR1}. }
}
\value{
A \code{numeric} of length one: the minus (!) LOOCV log-likelihood.
}
\note{
The minus LOOCV log-likelihood is returned as standard optimization procedures in R like \code{nlminb} and \code{constrOptim} minimize (rather then maximize). Hence, by providing the minus LOOCV log-likelihood the function \code{\link{loglikLOOCVVAR1}} can directly used by these optimization procedures.
}

\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data", \emph{Biometrical Journal}, 59(1), 172-191. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link[rags2ridges:ridgeP]{ridgeP}} and \code{\link{ridgeVAR1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A      <- createA(p, "chain")

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

## determine optimal values of the penalty parameters
\dontrun{optLambdas <- constrOptim(c(1,1), loglikLOOCVVAR1, gr=NULL, }
\dontrun{              ui=diag(2), ci=c(0,0), Y=Y, }
\dontrun{              control=list(reltol=0.01))$par }

# ridge ML estimation of the VAR(1) parameter estimates with 
# optimal penalty parameters
optLambdas <- c(0.1, 0.1)
ridgeVAR1(Y, optLambdas[1], optLambdas[2])$A
}

