\name{loglikLOOCVVAR1fused}
\alias{loglikLOOCVVAR1fused}
\title{
Leave-one-out (minus) cross-validated log-likelihood of multiple VAR(1) models
}
\description{
Evaluation of the (minus) leave-one-out cross-validated log-likelihood of multiple VAR(1) models for given choices of the (fused) ridge penalty parameters (\eqn{\lambda_a}, \eqn{\lambda_f} and \eqn{\lambda_{\omega}} for the autoregression coefficient matrice \eqn{\mathbf{A}_g}'s and the inverse error covariance matrix \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}, respectively). The functions also works with a (possibly) unbalanced experimental set-up. The VAR(1)-processes are assumed to have mean zero.
}
\usage{
loglikLOOCVVAR1fused(lambdas, Y, id, unbalanced=matrix(nrow=0, ncol=2), ...)
}
\arguments{
\item{lambdas}{ A \code{numeric} of length three, comprising positive values only. It contains the ridge penalty parameters to be used in the estimation of \eqn{\mathbf{A}} and the precision matrix of the errors, respectively. }
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{id}{ A vector with group indices comprising of integers only. First group is represented by '0', the next by '1', and so on until the last. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{...}{  Other arguments to be passed to \code{ridgeVAR1}. }
}
\value{
A \code{numeric} of length one: the minus (!) LOOCV log-likelihood.
}
\note{
The minus LOOCV log-likelihood is returned as standard optimization procedures in R like \code{nlminb} and \code{constrOptim} minimize (rather then maximize). Hence, by providing the minus LOOCV log-likelihood the function \code{\link{loglikLOOCVVAR1fused}} can directly used by these optimization procedures.
}
\references{
Miok, V., Wilting, S.M., Van Wieringen, W.N. (2017), ``Ridge estimation of network models from time-course omics data'', submitted. 
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link[rags2ridges:ridgeP]{ridgeP}} and \code{\link{ridgeVAR1fused}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points, G=groups)
p <- 3; n <- 12; T <- 10; G <- 3

# set model parameters
SigmaE <- matrix(1/2, p, p)
diag(SigmaE) <- 1
A1 <- -createA(p, "clique", nCliques=1, nonzeroA=0.1)
A2 <- t(createA(p, "chain", nBands=1,   nonzeroA=0.1))
A3 <- (A1 + A2) / 2

# generate data
Y1 <- dataVAR1(n/G, T, A1, SigmaE)
Y2 <- dataVAR1(n/G, T, A2, SigmaE)
Y3 <- dataVAR1(n/G, T, A3, SigmaE)
Y  <- abind::abind(Y1, Y2, Y3, along=3)
id <- c(rep(1, n/G), rep(2, n/G), rep(3, n/G))-1

## determine optimal values of the penalty parameters
\dontrun{optLambdas <- constrOptim(c(1,1,1), loglikLOOCVVAR1fused, gr=NULL, }
\dontrun{              ui=diag(3), ci=c(0,0,0), Y=Y, id=id}
\dontrun{              control=list(reltol=0.01))$par }

# ridge ML estimation of the VAR(1) parameter estimates with 
# optimal penalty parameters
optLambdas <- c(0.1, 0.1, 0.1)
VAR1hats <- ridgeVAR1fused(Y, id, optLambdas[1], optLambdas[2], optLambdas[3])
}

