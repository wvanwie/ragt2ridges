\name{motifStatsVAR1}
\alias{motifStatsVAR1}
\title{
Network motif detection for the VAR(1) model.
}
\description{
Function that detects standard network motifs from the lag-one autoregression relationships as implied by the VAR(1) model.
}
\usage{
motifStatsVAR1(sparseA, verbose=TRUE)
}
\arguments{
\item{sparseA}{ A \code{matrix} \eqn{\mathbf{A}} of autoregression parameters, which is assumed to be sparse. }
\item{verbose}{ A \code{logical} specifying whether summary output should be displayed. } 
}
\details{
Six types of motifs are detected: self-regulators, feedback pairs, feedforward loops, feedback loops, bi-fans, and diamonds (cf. Alon, 2007 for details). A detected motif is reported by the paths that 
constitute them. In line with Alon (2007), who distinguishes subtypes of these motived based on the sign of the path's contribution, the latter is also reported. When \code{verbose=TRUE} the summary output is also visualized. See this plot for the definition of the motifs in terms of the VAR(1) time series chain graph. 
}
\value{
An object of class \code{list} with slots:
\item{selfregulators}{ A \code{list} of motifs specified as matrices with each row one of the motif's paths and in the last column the sign of the path's contribution. }
\item{feedbackpairs}{  A \code{list} of motifs specified as matrices with each row one of the motif's paths and in the last column the sign of the path's contribution. }
\item{feedforwardloops}{ A \code{list} of motifs specified as matrices with each row one of the motif's paths and in the last column the sign of the path's contribution. }
\item{feedbackloops}{ A \code{list} of motifs specified as matrices with each row one of the motif's paths and in the last column the sign of the path's contribution. }
\item{bifans}{ A \code{list} of motifs specified as matrices with each row one of the motif's paths and in the last column the sign of the path's contribution. }
\item{diamonds}{ A \code{list} of motifs specified as matrices with each row one of the motif's paths and in the last column the sign of the path's contribution. }
Future versions of this function may include additional slots reporting more motif types.
}
\references{
Alon, U. (2007), "Network motifs: theory and experimental approaches", Nature Reviews Genetics, 8, 450-461.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>.
}
\seealso{
\code{\link{ridgeVAR1}}, \code{\link{sparsifyVAR1}}, \code{\link{graphVAR1}}, \code{\link{nodeStatsVAR1}}
}
\examples{
# specify lag one autoregression model
sparseA <- matrix(runif(2500), ncol=50)
sparseA[sparseA < 0.9] <- 0

# find motifs 
motifList <- motifStatsVAR1(sparseA)
}
