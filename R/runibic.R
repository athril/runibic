#' runibic: parallel row-based biclustering algorithm
#' for analysis of gene expression data in R
#' @docType package
#' @author Patryk Orzechowski \url{patryk.orzechowski@gmail.com}, Artur Pańszczyk \url{panszczyk.artur@gmail.com}
#' @name runibic
#' @useDynLib runibic
#' @importFrom Rcpp evalCpp
#' @import testthat SummarizedExperiment
#' @importFrom biclust biclust bicluster
#' @export runiDiscretize
#' @export set_runibic_params
#' @export runibic
#' @export BCUnibic
#' @export BCUnibicD
#' @description \code{\link{runibic}} is a package that contains much faster parallel version of one of the most accurate biclustering algorithms, UniBic.
#' The original method was reimplemented from C to C++11, OpenMP was added for parallelization.
#'
#' If you use this package, please cite it as:
#' Patryk Orzechowski, Artur Pańszczyk, Xiuzhen Huang, Jason H. Moore:
#' "runibic: a Bioconductor package for parallel row-based biclustering of gene expression data";
#' \emph{bioRxiv}, 2017,
#' doi: https://doi.org/10.1101/210682
#'
#' @details
#' For a given input matrix we first perform discretization and create index matrix using \code{\link{runiDiscretize}} function.
#' The discretization is performed taking into account quantiles of the data.
#' The resulting index matrix allows to detect order-preserving trends between each pair of the rows
#' irrespective to the order of columns.
#' After the ranking, the matrix is split by rows into subgroups based on the significance of the future biclusters.
#' In each of the chunks, we calculate pairwise calculations of Longest Common Subsequence LCS between all pairs of the rows.
#' LCS calculations are performed using dynamic programming and determine the longest order-preserving trend between each pair of the rows.
#' After partitioning the matrix strict order-preserving biclusters are determined and later expanded
#' to approximate-trend biclusters within \code{\link{cluster}} function.
#'
#' This package provides 3 main functions:
#' \code{\link{runibic}} and \code{\link{BCUnibic}} perform UniBic biclustering algorithm on numeric data, whilst
#' \code{\link{BCUnibicD}} could be applied to integer ones. The latter two methods are compatible with \code{\link[biclust]{Biclust}} class.
#'
#' @references Wang, Zhenjia, et al. "UniBic: Sequential row-based biclustering algorithm for analysis of gene expression data." Scientific reports 6 (2016): 23466.
#' @references Patryk Orzechowski, Artur Pańszczyk, Xiuzhen Huang, Jason H. Moore: "runibic: a Bioconductor package for parallel row-based biclustering of gene expression data", bioRxiv (2017): 210682, doi: https://doi.org/10.1101/210682
NULL


#' Class BCUnibic
#'
#' An S4 class to represent \code{\link{BCUnibic-class}} UniBic biclustering algorithm
#' for numeric input. The class is intended to use with 
#'
#' @name BCUnibic-class
#' @rdname BCUnibic-class
#' @seealso \code{\link{runibic}} 
setClass(Class = "BCUnibic", contains = "BiclustMethod", 
    prototype = prototype(biclustFunction = function(x, ...) {    
    runibic(x, ...)
}))


#' Class BCUnibicD
#'
#' An S4 class \code{\link{BCUnibicD-class}} defines UniBic biclustering algorithm
#' for discrete input.
#'
#' @name BCUnibicD-class
#' @rdname BCUnibicD-class
setClass(Class = "BCUnibicD", contains = "BiclustMethod",
    prototype = prototype(biclustFunction = function(x, ...) {
    runibic_d(x, ...)
}))


#' @describeIn runibic \code{\link{BCUnibic}} performs biclustering using UniBic on numeric matrix.
#' It is intended to use as a method called from \code{\link[biclust]{biclust}}.
BCUnibic <- function(x = NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0, useLegacy = FALSE) {
    if (is.null(x))
        return(methods::new("BCUnibic"))
    res <- biclust(x = x, method = BCUnibic(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0, useLegacy = FALSE)
    res@Parameters$Call <- match.call()
    return (res);
}


#' @describeIn runibic perform biclustering using UniBic on integer matrix.
#' It is intended to use as a method called from \code{\link[biclust]{biclust}}.
BCUnibicD <- function(x = NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0, useLegacy = FALSE) {
    if (is.null(x))
        return(methods::new("BCUnibicD"))
    res <- biclust(x = x, method = BCUnibicD(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0, useLegacy = FALSE)
    res@Parameters$Call <- match.call()
    return (res);
}


runibic_d <- function(x = NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0, useLegacy = FALSE) {
    MYCALL <- match.call()
    
    set_runibic_params(t, q, f, nbic, div, useLegacy)
    iX <- unisort(x)
    LCSRes <- calculateLCS(x, TRUE)
    res <- cluster(iX, x, LCSRes$lcslen, LCSRes$a, LCSRes$b, nrow(x), ncol(x) )
    return(biclust::BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
        matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
        res["info"]))
}
#' runibic
#'
#' Each of the following functions \code{\link{BCUnibic}}, \code{\link{BCUnibicD}}, 
#' \code{\link{runibic}} perform biclustering
#' using UniBic biclustering algorithm. The major difference
#' between the functions is that \code{\link{BCUnibicD}} require a discretized matrix, 
#' whilst \code{\link{BCUnibic}} (or \code{\link{runibic}})
#' could be applied to numeric one.
#'
#' @param x numeric or integer matrix (depends on the function)
#' @param t consistency level of the block (0.5-1.0].
#' @param q a double value for quantile discretization
#' @param f filtering overlapping blocks (default 1 do not remove any blocks)
#' @param nbic maximum number of biclusters in output
#' @param div number of ranks for up(down)-regulated genes: default: 0==ncol(x)
#' @param useLegacy boolean value for using legacy parameter settings
#' @return \code{\link[biclust]{Biclust}} object with detected biclusters
#'
#' @usage runibic(x = NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0, useLegacy=FALSE)
#' @seealso \code{\link{runiDiscretize}} \code{\link{set_runibic_params}} \code{\link{BCUnibic-class}} \code{\link{BCUnibicD-class}} \code{\link{unisort}}
#' @describeIn runibic perform biclustering using UniBic on numeric matrix.
#'
#' @examples
#' A <- matrix(replicate(100, rnorm(100)), nrow=100, byrow=TRUE)
#' runibic(A)
#' BCUnibic(A)
#' BCUnibic(A, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
#' B <- runiDiscretize(A)
#' runibic(B)
#' BCUnibicD(B, t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
#' biclust::biclust(A, method=BCUnibic(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
#' biclust::biclust(B, method=BCUnibicD(), t = 0.95, q = 0, f = 1, nbic = 100, div = 0)
runibic <- function(x = NULL, t = 0.95, q = 0, f = 1, nbic = 100, div = 0, useLegacy = FALSE) {
    if(inherits(x,"SummarizedExperiment")){
        x_d <- lapply(assays(x), runiDiscretize)
        return (lapply(x_d, runibic_d, t, q, f, nbic, div))
    }
    set_runibic_params(t, q, f, nbic, div, useLegacy)
    x_d <- runiDiscretize(x)
    return(runibic_d(x_d, t, q, f, nbic, div))
}


