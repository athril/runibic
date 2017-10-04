#' Parallel row-based biclustering algorithm for analysis of gene expression data in R
#' @docType package
#' @author Patryk Orzechowski
#' @name runibic
#' @useDynLib runibic
#' @import Rcpp biclust testthat
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"


runibic_d <- function(x, t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0) {

  MYCALL <- match.call()

  LCSRes = calculateLCS(unisort(x))
  x=t(x);
  res = cluster(x, LCSRes$lcslen,LCSRes$a,LCSRes$b, nrow(x), ncol(x) )
  return(biclust::BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                                matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
                                res["info"]))
}

#' Parallel row-based biclustering algorithm for analysis of gene expression data in R
#'
#' @usage runibic <- function(x, t = 0.95, q = 100, f = 1, nbic = 100, div = 0)

runibic <- function(x, t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0) {

  runibic_params(t,q,f,nbic,div)
  x_d <- discretize(x)
  return(runibic_d(x_d, t, q, f, nbic, div))

}