% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{runiDiscretize}
\alias{runiDiscretize}
\title{Discretize an input matrix}
\usage{
runiDiscretize(x)
}
\arguments{
\item{x}{a numeric matrix}
}
\value{
a discretized matrix containing integers only
}
\description{
This function discretizes the input matrix. 
\code{\link{runiDiscretize}} uses paramaters: 'div' and 'q', 
which are set by set_runibic_params function.
The funtion returns a discrete matrix with given number of ranks 
based on the parameter div. In contrast to biclust::discretize
the function takes into consideration the quantile parameter 'q'. 
When 'q' parameter is higher or equal 0.5 a simple discretization is used
with equal sizes of the levels using the quantiles. If 'q' parameter 
is lower than 0.5 we use up(down)-regulated discretization divided
into three parts.
}
\examples{
A <- replicate(10, rnorm(20))
runiDiscretize(A)

}
\seealso{
\code{\link{set_runibic_params}} \code{\link{calculateLCS}} \code{\link[biclust]{discretize}}
}
