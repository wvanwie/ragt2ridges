\name{evaluateVAR1fit}
\alias{evaluateVAR1fit}
\title{
Visualize the fit of a VAR(1) model
}
\description{
Simple plots for assessment of the fit of an estimated VAR(1) model. 
}
\usage{
evaluateVAR1fit(Y, A, SigmaE, unbalanced=NULL, diag=FALSE, 
                fileType="eps", dir=getwd())
}
\arguments{
\item{Y}{ Three-dimensional \code{array} containing the data. The first, second and third dimensions correspond to covariates, time and samples, respectively. The data are assumed to be centered covariate-wise. }
\item{A}{ A \code{matrix} \eqn{\mathbf{A}} of autoregression parameters. }
\item{SigmaE}{ Covariance \code{matrix} of the errors (innovations). }
\item{unbalanced}{ A \code{matrix} with two columns, indicating the unbalances in the design. Each row represents a missing design point in the (time x individual)-layout. The first and second column indicate the time and individual (respectively) specifics of the missing design point. }
\item{diag}{ A  \code{logical}, should the diagonal be included in the evaluation of the fit. }
\item{fileType}{ A \code{character} specifying the format in which figures should be save. Either '\code{pdf}' or '\code{eps}'. }
\item{dir}{ A \code{character} specifying the directory where plots should be saved. }
}
\value{
Plots are saved in the specified directory.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>
}
\seealso{
\code{\link[rags2ridges:evaluateSfit]{evaluateSfit}}, \code{\link{ridgeVAR1}}, \code{\link{dataVAR1}}.
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

# evaluate fit VAR(1) model
\dontrun{ evaluateVAR1fit(Y, VAR1hat$A, symm(VAR1hat$P)) }
}

