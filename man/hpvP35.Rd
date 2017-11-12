\name{hpvP53}
\alias{hpvP53}
\alias{hpvP53rna}
\alias{hpvP53mir}
\alias{hpvP53cn}
\alias{mir2rna}
\alias{cn2rna}
\docType{data}

\title{ Time-course P53 pathway data }

\description{
Time-course omics data of the P53-signalling pathway during HPV-induced transformation. Four cell lines transfected with the human papilloma virus (HPV) are cultured in petri dishes. At eight time points 
the cells from each cell line are molecularly interrogated genomically and transcriptomically (both mRNA and microRNA expression). From each molecular levels data DNA segments/mRNAs/microRNAs mapping to the P53-signalling pathway (as defined by KEGG, Ogata et al., 1999) are included. The included microRNAs are predicted to target the mRNAs of the pathway. Finally, two adjacency matrices are included. One describes how the DNA copy number features link to that of the mRNA expression data. The map is based on the feature's genomic location, while the other specifies which microRNAs (columns) target which mRNA (rows) according target predicton data bases. Confer Wilting et al. (2017) for preprocessing details.
}
\usage{
data(hpvP53)
}
\value{
\item{hpvP53rna}{ An object of class \code{\link[Biobase:ExpressionSet]{ExpressionSet}} with the time-course mRNA expression data of \eqn{32 (=4 \mbox{ cell lines } \times 8 \mbox{ time points})} observations. }
\item{hpvP53mir}{ An object of class \code{\link[Biobase:ExpressionSet]{ExpressionSet}} with the time-course microRNA expression data of \eqn{32 (=4 \mbox{ cell lines } \times 8 \mbox{ time points})} observations. }
\item{hpvP53cn}{ An object of class \code{\link[CGHbase:cghRaw]{cghRaw}} with the time-course DNA copy number data of \eqn{32 (=4 \mbox{ cell lines } \times 8 \mbox{ time points})} observations. }
\item{mir2rna}{ An object of class \code{matrix}. Adjancency matrix indicating which microRNAs (features in \code{hpvP53mir}) target which mRNAs (features in \code{hpvP53rna}). Rows = mRNAs; columns = microRNAs. }
\item{cn2rna}{ An object of class \code{matrix}. Adjancency matrix indicating which DNA copy number features (in \code{hpvP53cn}) map to which mRNAs (in \code{hpvP53rna}). }
}
\source{
Ogata, H., Goto, S., Sato, K., Fujibuchi, W., Bono, H., Kanehisa, M. (1999), KEGG: Kyoto Encyclopedia of Genes and Genomes. \emph{Nucleic Acids Research}, 27(1), 29-34.

Wilting, S.M., Miok, V., Jaspers, A., Boon, D., Hanne Sorgard, H., Lando, M., Snoek, B.C., Van Wieringen, W.N., Meijer, C.J.L.M., Lyng, H., Snijders, P.J.F., Steenbergen, R.D.M. (2016), "Aberrant methylation-mediated silencing of microRNAs contributes to HPV-induced anchorage independence", \emph{Oncotarget}, Vol. 7, No. 28, 43805-43819, DOI:10.18632/oncotarget.9698 .
}
\examples{
data(hpvP53)
}
\keyword{datasets}



