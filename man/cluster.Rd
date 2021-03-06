% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cluster}
\alias{cluster}
\title{Calculate biclusters from sorted list of LCS scores and row indices}
\usage{
cluster(discreteInput, discreteInputValues, scores, geneOne, geneTwo, rowNumber,
  colNumber)
}
\arguments{
\item{discreteInput}{an integer matrix with indices of sorted columns}

\item{discreteInputValues}{an integer matrix with discrete values}

\item{scores}{a numeric vector with LCS length}

\item{geneOne}{a numeric vector with first row indexes 
from pairwise LCS calculation}

\item{geneTwo}{a numeric vector with second row indexes 
from pairwise LCS calculation}

\item{rowNumber}{a int with number of rows in the input matrix}

\item{colNumber}{a int with number of columns in the input matrix}
}
\value{
a list with information of found biclusters
}
\description{
This function search for biclusters in the input matrix. 
The calculations are based on the integer matrix with indexes 
indicating positions of j-th smallest element in each row 
and the results from calculations of Longest Common Subsequence 
between all rows in the input matrix. The paramteres of this function can be
obtained from other functions provided by this package.
}
\examples{
A <- matrix( c(4,3,1,2,5,8,6,7,9,10,11,12),nrow=4,byrow=TRUE)
iA <- unisort(A)
lcsResults <- calculateLCS(A)
cluster(iA, A, lcsResults$lcslen, lcsResults$a, lcsResults$b, nrow(A), ncol(A))

}
\seealso{
\code{\link{runibic}} \code{\link{calculateLCS}} \code{\link{unisort}}
}
