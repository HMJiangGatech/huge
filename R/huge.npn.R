#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.npn(): nonparanormal(npn) transofmration                         #
#-----------------------------------------------------------------------#

#' Nonparanormal(npn) transformation
#'
#' Implements the Gausianization to help relax the assumption of normality.
#'
#' The nonparanormal extends Gaussian graphical models to semiparametric Gaussian copula models.Motivated by sparse additive models, the nonparanormal method estimates the Gasussian copula by marginally transforming the variables using smooth functions.Computationally, the estimation of a nonparanormal transformation is very efficient and only requires one pass of the data matrix.
#'
#' @param x The \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} dimensions
#' @param npn.func The transformation function used in the npn transformation. If \code{npn.func = "truncation"}, the truncated ECDF is applied. If \code{npn.func = "shrinkage"}, the shrunken ECDF is applied. The default is \code{"shrinkage"}. If \code{npn.func = "skeptic"}, the nonparanormal skeptic is applied.
#' @param npn.thresh The truncation threshold used in nonparanormal transformation, ONLY applicable when \code{npn.func = "truncation"}. The default value is \code{1/(4*(n^0.25)*} \code{sqrt(pi*log(n)))}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @return
#' \item{data}{
#' A \code{d} by \code{d} nonparanormal correlation matrix if \code{npn.func = "skeptic"}, and A \code{n} by \code{d} data matrix representing \code{n} observations in \code{d} transformed dimensions other wise.
#' }
#' @seealso \code{\link{huge}} and \code{\link{huge-package}}.
#' @examples
#' # generate nonparanormal data
#' L = huge.generator(graph = "cluster", g = 5)
#' L$data = L$data^5
#'
#' # transform the data using the shrunken ECDF
#' Q = huge.npn(L$data)
#'
#' # transform the non-Gaussian data using the truncated ECDF
#' Q = huge.npn(L$data, npn.func = "truncation")
#'
#' # transform the non-Gaussian data using the truncated ECDF
#' Q = huge.npn(L$data, npn.func = "skeptic")
#' @export
huge.npn = function(x, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE){
  gcinfo(FALSE)
  n = nrow(x)
    d = ncol(x)
    x.col = colnames(x)
    x.row = rownames(x)

    # Shrinkaage transformation
  if(npn.func == "shrinkage"){
    if(verbose) cat("Conducting the nonparanormal (npn) transformation via shrunkun ECDF....")

    x = qnorm(apply(x,2,rank)/(n+1))
    x = x/sd(x[,1])

    if(verbose) cat("done.\n")
    rm(n,d,verbose)
       gc()
       colnames(x) = x.col
    rownames(x) = x.row
  }

  # Truncation transformation
  if(npn.func == "truncation"){
    if(verbose) cat("Conducting nonparanormal (npn) transformation via truncated ECDF....")
    if(is.null(npn.thresh)) npn.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))

    x = qnorm(pmin(pmax(apply(x,2,rank)/n, npn.thresh), 1-npn.thresh))
      x = x/sd(x[,1])

      if(verbose) cat("done.\n")
      rm(n,d,npn.thresh,verbose)
       gc()
       colnames(x) = x.col
    rownames(x) = x.row
  }

  if(npn.func == "skeptic"){
    if(verbose) cat("Conducting nonparanormal (npn) transformation via skeptic....")
    x = 2*sin(pi/6*cor(x,method="spearman"))
      if(verbose) cat("done.\n")
      rm(n,d,verbose)
       gc()
       colnames(x) = x.col
    rownames(x) = x.col
  }

  return(x)
}
