#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# glasso(): The graphical lasso (glasso) using sparse matrix output     #
#-----------------------------------------------------------------------#

#' The graphical lasso (glasso) using sparse matrix output 
#' 
#' See more details in \code{\link{huge}}
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A sequence of decreasing positive numbers to control the regularization when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, or the thresholding in \code{method = "ct"}. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. When \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, use with care - it is better to supply a decreasing sequence values than a single (small) value.
#' @param nlambda The number of regularization/thresholding parameters. The default value is \code{30} for \code{method = "ct"} and \code{10} for \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}.
#' @param lambda.min.ratio If \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, it is the smallest value for \code{lambda}, as a fraction of the upperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, and 0.05 \code{method = "ct"}.
#' @param scr If \code{scr = TRUE}, the lossy screening rule is applied to preselect the neighborhood before the graph estimation. The default value is  \code{FALSE}. NOT applicable when \code{method = "ct"}, {"mb"}, or {"tiger"}.
#' @param cov.output If \code{cov.output = TRUE}, the output will include a path of estimated covariance matrices. ONLY applicable when \code{method = "glasso"}. Since the estimated covariance matrices are generally not sparse, please use it with care, or it may take much memory under high-dimensional setting. The default value is \code{FALSE}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @export
huge.glasso = function(x, lambda = NULL, lambda.min.ratio = NULL, nlambda = NULL, scr = NULL, cov.output = FALSE, verbose = TRUE){

  gcinfo(FALSE)
  n = nrow(x)
  d = ncol(x)
  cov.input = isSymmetric(x)
  if(cov.input)
  {
    if(verbose) cat("The input is identified as the covariance matrix.\n")
    S = x
  }
  else
  {
    x = scale(x)
    S = cor(x)
  }
  rm(x)
  gc()
  if(is.null(scr)) scr = FALSE
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda))
  {
    if(is.null(nlambda))
      nlambda = 10
    if(is.null(lambda.min.ratio))
      lambda.min.ratio = 0.1
    lambda.max = max(max(S-diag(d)),-min(S-diag(d)))
    lambda.min = lambda.min.ratio*lambda.max
    lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }

  fit = .Call("_huge_hugeglasso",S,lambda,scr,verbose,cov.output)

  fit$scr = scr
  fit$lambda = lambda
  fit$cov.input = cov.input
  fit$cov.output = cov.output

  rm(S)
  gc()
  if(verbose){
       cat("\nConducting the graphical lasso (glasso)....done.                                          \r")
       cat("\n")
      flush.console()
  }
  return(fit)
}
