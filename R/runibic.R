#' runibic: parallel row-based biclustering algorithm
#' for analysis of gene expression data in R
#' @docType package
#' @author Patryk Orzechowski, Artur Pańszczyk
#' @name runibic
#' @useDynLib runibic
#' @importFrom Rcpp evalCpp
#' @import testthat SummarizedExperiment
#' @importFrom biclust biclust bicluster
#' @export runibic
#' @export BCUnibic
#' @export BCUnibicD
#' @export runiDiscretize
#' @export set_runibic_params
#'
NULL



#' Class BCUnibic
#'
#' Class \code{BCUnibic} defines UniBic biclustering algorithm
#' for numeric input
#'
#' @name BCUnibic-class
#' @rdname BCUnibic-class
setClass(Class = "BCUnibic", contains = "BiclustMethod", 
    prototype = prototype(biclustFunction = function(x, ...) {
    runibic(x, ...)
}))



#' Class BCUnibicD
#'
#' Class \code{BCUnibicD} defines UniBic biclustering algorithm
#' for discrete input
#'
#' @name BCUnibicD-class
#' @rdname BCUnibicD-class
setClass(Class = "BCUnibicD", contains = "BiclustMethod",
    prototype = prototype(biclustFunction = function(x, ...) {
    runibic_d(x, ...)
}))



#' BCUnibic
#'
#' Both \code{BCUnibic} and \code{BCUnibicD} function perform biclustering
#' using UniBic biclustering algorithm. The major difference 
#' between both functions is that BCUnibicD requires a discretized matrix,
#' whilst BCUnibic could be applied to numeric one.
#'
#' @param x numeric matrix
#' @param t consistency level of the block (0.5-1.0].
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks (default 1 do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks for up(down)-regulated genes: default: 0==ncol(x)
#' @return [biclust::Biclust](https://cran.r-project.org/web/packages/biclust/biclust.pdf) object with detected biclusters
#'
#' @usage BCUnibic(x, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
#'
#' @examples
#' A <- replicate(10, rnorm(20))
#' BCUnibic(A)
#' BCUnibic(A, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
#' biclust::biclust(A, method=BCUnibic(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
BCUnibic <- function(x=NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {
    if (is.null(x))
        return(methods::new("BCUnibic"))
    res <- biclust(x = x, method = BCUnibic(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
    res@Parameters$Call <- match.call()
    return (res);
}




#' BCUnibicD
#'
#' Both \code{BCUnibic} and \code{BCUnibicD} function perform biclustering
#' using UniBic biclustering algorithm. The major difference 
#' between both functions is that BCUnibicD requires a discretized matrix, 
#' whilst BCUnibic could be applied to numeric one.
#'
#' @param x a discrete matrix
#' @param t consistency level of the block (0.5-1.0]
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks (default 1 do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks for up(down)-regulated genes: default: 0==ncol(x)
#' @return [biclust::Biclust](https://cran.r-project.org/web/packages/biclust/biclust.pdf) object with detected biclusters
#'
#' @seealso BCUnibic biclustering algrorithm for discrete data
#'
#' @usage BCUnibicD(x, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
#' @examples
#' A <- runiDiscretize(replicate(10, rnorm(20)))
#' BCUnibicD(A)
#' BCUnibicD(A, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
#' biclust::biclust(A, method=BCUnibicD(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
BCUnibicD <- function(x=NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {
    if (is.null(x)) 
        return(methods::new("BCUnibicD"))
    res <- biclust(x = x, method = BCUnibicD(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
    res@Parameters$Call <- match.call()
    return (res);
}


#' runibic_d
#' @describeIn runibic biclustering algrorithm for discrete data
runibic_d <- function(x, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {
    MYCALL <- match.call()

    iX <- unisort(x)
    LCSRes <- calculateLCS(x, TRUE)
    res <- cluster(iX, x, LCSRes$lcslen, LCSRes$a, LCSRes$b, nrow(x), ncol(x) )
    return(biclust::BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
        matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
        res["info"]))
}


#' runibic
#' @param x numeric matrix or SummarizedExperiment
#' @param t consistency level of the block (0.5-1.0].
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks, default 1(do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks for up(down)-regulated genes: default: 0==ncol(x)
#' @return Biclust object with detected biclusters
#'
#' @usage runibic(x, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
#' @examples 
#' A=matrix(replicate(100, rnorm(100)), nrow=100, byrow=TRUE)
#' runibic(A)
runibic <- function(x, t = 0.95, q = 0, f = 1, nbic = 100, div = 0) {
    if(inherits(x,"SummarizedExperiment")){
        x_d <- lapply(assays(x), runiDiscretize)
        return (lapply(x_d, runibic_d, t, q, f, nbic, div))
    }
    set_runibic_params(t, q, f, nbic, div)
    x_d <- runiDiscretize(x)
    return(runibic_d(x_d, t, q, f, nbic, div))
}


