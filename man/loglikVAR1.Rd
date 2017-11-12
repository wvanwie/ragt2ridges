\name{loglikVAR1}
\alias{loglikVAR1}
\title{
Log-likelihood of the VAR(1) model.
}
\description{
Log-likelihood of the VAR(1) model specified by the supplied parameters
}
\usage{
loglikVAR1(Y, A, P, unbalanced=matrix(nrow=0, ncol=2))
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{A}{ A \code{matrix} \eqn{\mathbf{A}} of autoregression parameters. }
\item{P}{ Inverse error covariance \code{matrix} \eqn{\mathbf{\Omega}_{\varepsilon} (=\mathbf{\Sigma_{\varepsilon}^{-1}})}. }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
}
\value{
The log-likelihood of the VAR(1) model with supplied parameters.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}

\seealso{
\code{\link{ridgeVAR1}}.
}
\examples{
# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
SigmaE <- diag(p)/4
A      <- createA(p, "chain")

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

# center data
Y <- centerVAR1data(Y)

# fit VAR(1) model
VAR1hat <- ridgeVAR1(Y, 1, 1)

# evaluate the log-likelihood of this fit.
loglikVAR1(Y, VAR1hat$A, VAR1hat$P)
}

