% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{unisort}
\alias{unisort}
\title{Computing the indexes of j-th smallest values of each row}
\usage{
unisort(x)
}
\arguments{
\item{x}{a integer matrix}
}
\value{
a integer matrix with indexes indicating positions 
of j-th smallest element in each row
}
\description{
This function sorts separately each row of a integer matrix 
and returns a matrix in which the value in i-th row and j-th column 
represent the index of the j-th smallest value of the i-th row.
}
\examples{
A <- matrix(c(4, 3, 1, 2, 5, 8, 6, 7), nrow=2, byrow=TRUE)
unisort(A)
}
\seealso{
\code{\link{runibic}} \code{\link{calculateLCS}} \code{\link{runiDiscretize}}
}
