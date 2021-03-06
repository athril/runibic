% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calculateLCS}
\alias{calculateLCS}
\title{Calculate all Longest Common Subsequences between a matrix.}
\usage{
calculateLCS(discreteInput, useFibHeap = TRUE)
}
\arguments{
\item{discreteInput}{is a input discrete matrix}

\item{useFibHeap}{boolean value for choosing which sorting method 
should be used in sorting of output}
}
\value{
a list with sorted values based on calculation of the length of LCS
between all pairs of rows
}
\description{
This function computes unique pairwise Longest Common Subsequences 
between each row of input matrix. The function outputs a list 
sorted by Longest Common Subsequences (LCS) length. The output list contains
the length of calculated LCS, indices, of the first and second rows
between which LCS was calculated.
The function uses two different sorting methods. The default one 
uses Fibonacci Heap used in original implementation of Unibic, 
the second one uses standard sorting algorithm from C++ STL.
}
\examples{
A <- matrix(c(4, 3, 1, 2, 5, 8, 6, 7), nrow=2, byrow=TRUE)
calculateLCS(A, TRUE)
}
\seealso{
\code{\link{runibic}} \code{\link{backtrackLCS}}  \code{\link{pairwiseLCS}}
}
