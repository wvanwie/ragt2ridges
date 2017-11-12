\name{loglikLOOCVVAR2}
\alias{loglikLOOCVVAR2}
\title{
Leave-one-out (minus) cross-validated log-likelihood of VAR(2) model
}
\description{
Evaluation of the (minus) leave-one-out cross-validated log-likelihood of the VAR(2) model for given choices of the ridge penalty parameters (\eqn{\lambda_{a1}}, \eqn{\lambda_{a2}} and \eqn{\lambda_{\omega}} for the lag one autoregression coefficient matrix \eqn{\mathbf{A}_1}, lag two autoregression coefficient matrix \eqn{\mathbf{A}_2} of time-varying covariates, and the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}, respectively). The functions also works with a (possibly) unbalanced experimental set-up. The VAR(2)-process is assumed to have mean zero.
}
\usage{
loglikLOOCVVAR2(lambdas, Y, unbalanced=matrix(nrow=0, ncol=2), ...)
}
\arguments{
\item{lambdas}{ A \code{numeric} of length three, comprising positive values only. It contains the ridge penalty parameters to be used in the estimation of \eqn{\mathbf{A}_1}, \eqn{\mathbf{A}_2} and the precision matrix of the errors. }
\item{Y}{ Three-dimensional \code{array} containing the response data. The first, second and third dimensions correspond to variates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{...}{  Other arguments to be passed to the \code{\link{ridgeVAR2}}-function. }
}
\value{
A \code{numeric} of length one: the minus (!) LOOCV log-likelihood.
}
\note{
The minus LOOCV log-likelihood is returned as standard optimization procedures in R like \code{nlminb} and \code{constrOptim} minimize (rather then maximize). Hence, by providing the minus LOOCV log-likelihood the function \code{\link{loglikLOOCVVAR2}} can directly used by these optimization procedures.
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link[rags2ridges:ridgeP]{ridgeP}} and \code{\link{ridgeVAR2}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 12; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A1 <- createA(p, topology="clique", nonzeroA=0.1, nClique=1)
A2 <- createA(p, topology="hub", nonzeroA=0.1, nHubs=1)

# generate data
Y <- dataVAR2(n, T, A1, A2, SigmaE)

## determine optimal values of the penalty parameters
\dontrun{optLambdas <- constrOptim(c(1,1, 1), loglikLOOCVVAR2, gr=NULL, }
\dontrun{              ui=diag(3), ci=c(0,0,0), Y=Y,}
\dontrun{              control=list(reltol=0.01))$par }

# ridge ML estimation of the VAR(2) parameter estimates with 
# optimal penalty parameters
optLambdas <- c(0.1, 0.1, 0.1)
ridgeVAR2(Y, optLambdas[1], optLambdas[2], optLambdas[3])
}

