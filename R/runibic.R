#' Parallel row-based biclustering algorithm for analysis of gene expression data in R
#' @docType package
#' @author Patryk Orzechowski, Artur Pańszczyk
#' @name runibic
#' @useDynLib runibic
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp evalCpp
#' @import testthat SummarizedExperiment
#' @importFrom biclust biclust bicluster
NULL



#' Class BCUnibic
#'
#' Class BCUnibic performs biclustering on continuous matrix using UniBic biclustering algorithm
#'
#' @name BCUnibic-class
#' @rdname BCUnibic-class
setClass(Class = "BCUnibic", contains = "BiclustMethod",
         prototype = prototype(biclustFunction = function(x, ...) {
  runibic(x, ...)
}))



#' Class BCUnibicD
#'
#' Class BCUnibicD performs biclustering discrete matrix using UniBic biclustering algorithm.
#'
#' @name BCUnibicD-class
#' @rdname BCUnibicD-class
setClass(Class = "BCUnibicD", contains = "BiclustMethod",
         prototype = prototype(biclustFunction = function(x, ...) {
  runibic_d(x, ...)
}))



#' BCUnibicD
#'
#' \code{BCUnibicD} performs Unibic for a discrete matrix.
#'
#' @aliases UnibicD biclust,matrix,BCUnibicD-method
#'
#' @describeIn Unibic biclustering algrorithm for discrete data
#'
#' @examples
#' BCUnibicD(discretize(replicate(10, rnorm(20))))
#' @usage \S4method{biclust}{matrix,BCUnibicD}(x=NULL, method=BCUnibicD(), t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0)
BCUnibicD <- function(x=NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {
  if (is.null(x)) 
    return(methods::new("BCUnibicD"))
  res<- biclust(x = x, method = BCUnibicD(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
  res@Parameters$Call = match.call()
  return (res);
}


#' BCUnibic
#'
#' @aliases Unibic biclust,matrix,BCUnibic-method
#'
#' @describeIn Unibic biclustering algrorithm
#' @param x numeric matrix
#' @param t consistency level of the block (0.5-1.0].
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks, default 1(do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks as which we treat the up(down)-regulated value: default: 0==ncol(x)
#' @return Biclust object with detected biclusters
#' @param method Unibic method for continuous data
#'
#' @examples
#' BCUnibic(replicate(10, rnorm(20)))
#' @usage \S4method{biclust}{matrix,BCUnibic}(x=NULL, method=BCUnibic(), t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0)
BCUnibic <- function(x=NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {
  if (is.null(x)) 
    return(methods::new("BCUnibic"))
  res<- biclust(x = x, method = BCUnibic(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
  res@Parameters$Call = match.call()
  return (res);
}


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
runibic_d <- function(x, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {

  MYCALL <- match.call()

  iX = unisort(x)
  LCSRes = calculateLCS(x,TRUE)
  res = cluster(iX,x, LCSRes$lcslen,LCSRes$a,LCSRes$b, nrow(x), ncol(x) )
  return(biclust::BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                                matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
                                res["info"]))
}



#' Parallel row-based biclustering algorithm for analysis of gene expression data in R
#' @param x numeric matrix or SummarizedExperiment
#' @param t consistency level of the block (0.5-1.0].
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks, default 1(do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks as which we treat the up(down)-regulated value: default: 0==ncol(x)
#' @return Biclust object with detected biclusters
#'
#' @usage runibic(x, t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0)
#' @examples 
#' A=matrix(replicate(100, rnorm(100)), nrow=100, byrow=TRUE)
#' runibic(A)
runibic <- function(x, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {
  if(inherits(x,"SummarizedExperiment")){
    return (runibic_se(x,t,q,f,nbic,div))
  }
  runibic_params(t,q,f,nbic,div)
  x_d <- discretize(x)
  return(runibic_d(x_d, t, q, f, nbic, div))
}

#' Parallel row-based biclustering algorithm for analysis of gene expression data in R
#' @param x SummarizedExperiment
#' @param t consistency level of the block (0.5-1.0].
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks, default 1(do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks as which we treat the up(down)-regulated value: default: 0==ncol(x)
#' @return List of Biclust objects with detected bicluster for each result from every assay within input
#'
#' @usage runibic_se(x, t = 0.95, q = 0.5, f = 1, nbic = 100, div = 0)
#' @examples 
#' A=matrix(replicate(100, rnorm(100)), nrow=100, byrow=TRUE)
#' se = SummarizedExperiment(assays=list(counts=A))
#' runibic_se(se)
runibic_se <- function(x, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {

  if(!inherits(x,"SummarizedExperiment")) stop("x must be a SummarizedExperiment")
  runibic_params(t,q,f,nbic,div)
  x_d <- lapply(assays(x), discretize)
  return(lapply(x_d, runibic_d, t, q, f, nbic, div))
}

