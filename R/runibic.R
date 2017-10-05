#' Parallel row-based biclustering algorithm for analysis of gene expression data in R
#' @docType package
#' @author Patryk Orzechowski
#' @name runibic
#' @useDynLib runibic
#' @import Rcpp biclust testthat
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL




#' Parallel row-based biclustering algorithm for analysis of gene expression data in R
#' @param x integer matrix
#' @param t consistency level of the block (0.5-1.0].
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks, default 1(do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks as which we treat the up(down)-regulated value: default: 0==ncol(x)
#' @return Biclust object with detected biclusters
#'
#' @usage runibic_d(x, t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0)

runibic_d <- function(x, t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0) {

  MYCALL <- match.call()

  LCSRes = calculateLCS(unisort(x),TRUE)
  x=t(x);
  res = cluster(x, LCSRes$lcslen,LCSRes$a,LCSRes$b, nrow(x), ncol(x) )
  return(biclust::BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                                matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
                                res["info"]))
}

#' Parallel row-based biclustering algorithm for analysis of gene expression data in R
#' @param x numeric matrix
#' @param t consistency level of the block (0.5-1.0].
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks, default 1(do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks as which we treat the up(down)-regulated value: default: 0==ncol(x)
#' @return Biclust object with detected biclusters
#'
#' @usage runibic(x, t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0)
#' @examples 
#' A=matrix(replicate(10, rnorm(20)), nrow=10, byrow=TRUE)
#' runibic(A)

runibic <- function(x, t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0) {

  runibic_params(t,q,f,nbic,div)
  x_d <- discretize(x)
  return(runibic_d(x_d, t, q, f, nbic, div))

}