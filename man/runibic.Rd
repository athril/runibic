% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runibic.R
\docType{package}
\name{runibic}
\alias{runibic}
\alias{runibic-package}
\alias{BCUnibic}
\alias{BCUnibicD}
\alias{runibic}
\title{runibic: parallel row-based biclustering algorithm
for analysis of gene expression data in R}
\usage{
BCUnibic(x = NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0,
  useLegacy = FALSE)

BCUnibicD(x = NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0,
  useLegacy = FALSE)

runibic(x = NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0, useLegacy=FALSE)
}
\arguments{
\item{x}{numeric or integer matrix (depends on the function)}

\item{t}{consistency level of the block (0.5-1.0].}

\item{q}{a double value for quantile discretization}

\item{f}{filtering overlapping blocks (default 1 do not remove any blocks)}

\item{nbic}{maximum number of biclusters in output}

\item{div}{number of ranks for up(down)-regulated genes: default: 0==ncol(x)}

\item{useLegacy}{boolean value for using legacy parameter settings}
}
\value{
\code{\link[biclust]{Biclust}} object with detected biclusters
}
\description{
\code{\link{runibic}} is a package that contains much faster parallel version of one of the most accurate biclustering algorithms, UniBic.
The original method was reimplemented from C to C++11, OpenMP was added for parallelization.

If you use this package, please cite it as:
Patryk Orzechowski, Artur Pańszczyk, Xiuzhen Huang, Jason H Moore; 
"runibic: a Bioconductor package for parallel row-based biclustering of gene expression data";
Bioinformatics, 2018, bty512,
doi: https://doi.org/10.1093/bioinformatics/bty512

Each of the following functions \code{\link{BCUnibic}}, \code{\link{BCUnibicD}}, 
\code{\link{runibic}} perform biclustering
using UniBic biclustering algorithm. The major difference
between the functions is that \code{\link{BCUnibicD}} require a discretized matrix, 
whilst \code{\link{BCUnibic}} (or \code{\link{runibic}})
could be applied to numeric one.
}
\details{
For a given input matrix we first perform discretization and create index matrix using \code{\link{runiDiscretize}} function.
The discretization is performed taking into account quantiles of the data.
The resulting index matrix allows to detect order-preserving trends between each pair of the rows
irrespective to the order of columns.
After the ranking, the matrix is split by rows into subgroups based on the significance of the future biclusters.
In each of the chunks, we calculate pairwise calculations of Longest Common Subsequence LCS between all pairs of the rows.
LCS calculations are performed using dynamic programming and determine the longest order-preserving trend between each pair of the rows.
After partitioning the matrix strict order-preserving biclusters are determined and later expanded
to approximate-trend biclusters within \code{\link{cluster}} function.

This package provides 3 main functions:
\code{\link{runibic}} and \code{\link{BCUnibic}} perform UniBic biclustering algorithm on numeric data, whilst
\code{\link{BCUnibicD}} could be applied to integer ones. The latter two methods are compatible with \code{\link[biclust]{Biclust}} class.
}
\section{Functions}{
\itemize{
\item \code{BCUnibic}: \code{\link{BCUnibic}} performs biclustering using UniBic on numeric matrix.
It is intended to use as a method called from \code{\link[biclust]{biclust}}.

\item \code{BCUnibicD}: perform biclustering using UniBic on integer matrix.
It is intended to use as a method called from \code{\link[biclust]{biclust}}.

\item \code{runibic}: perform biclustering using UniBic on numeric matrix.
}}

\examples{
A <- matrix(replicate(100, rnorm(100)), nrow=100, byrow=TRUE)
runibic(A)
BCUnibic(A)
BCUnibic(A, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
B <- runiDiscretize(A)
runibic(B)
BCUnibicD(B, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
biclust::biclust(A, method=BCUnibic(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
biclust::biclust(B, method=BCUnibicD(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
}
\references{
Wang, Zhenjia, et al. "UniBic: Sequential row-based biclustering algorithm for analysis of gene expression data." Scientific reports 6 (2016): 23466.

Patryk Orzechowski, Artur Pańszczyk, Xiuzhen Huang, Jason H. Moore: "runibic: a Bioconductor package for parallel row-based biclustering of gene expression data", bioRxiv (2017): 210682, doi: https://doi.org/10.1101/210682
}
\seealso{
\code{\link{runiDiscretize}} \code{\link{set_runibic_params}} \code{\link{BCUnibic-class}} \code{\link{BCUnibicD-class}} \code{\link{unisort}}
}
\author{
Patryk Orzechowski \url{patryk.orzechowski@gmail.com}, Artur Pańszczyk \url{panszczyk.artur@gmail.com}
}
