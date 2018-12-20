\name{pruneMotifStats}
\alias{pruneMotifStats}
\title{
Network motif list subsetting.
}
\description{
Limit a list of network motifs as returned by the \code{motifStatsVAR1}-function to those involving a particular node.
}
\usage{
pruneMotifStats(motifList, id)
}
\arguments{
\item{motifList}{ A \code{list} with motifs as returned by the \code{motifStatsVAR1}-function. }
\item{id}{ A \code{integer} representing the node of interest. } 
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
Alon, U. (2007), ``Network motifs: theory and experimental approaches'', Nature Reviews Genetics, 8, 450-461.
}
\author{
Wessel N. van Wieringen <w.vanwieringen@vumc.nl>.
}
\seealso{
\code{\link{graphVAR1}}, \code{\link{motifStatsVAR1}}
}
\examples{
# specify lag one autoregression model
sparseA <- matrix(runif(2500), ncol=50)
sparseA[sparseA < 0.9] <- 0

# find motifs 
motifList <- motifStatsVAR1(sparseA)

# prune motif list
pruneMotifStats(motifList, 1)
}
