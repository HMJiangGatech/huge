#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.gect(): graph estimation via correlation thresholding (ct)       #
#-----------------------------------------------------------------------#

#' Graph estimation via correlation thresholding (ct) 
#' 
#' See more details in \code{\link{huge}}
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A sequence of decreasing positive numbers to control the regularization when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, or the thresholding in \code{method = "ct"}. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. When \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, use with care - it is better to supply a decreasing sequence values than a single (small) value.
#' @param nlambda The number of regularization/thresholding parameters. The default value is \code{30} for \code{method = "ct"} and \code{10} for \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}.
#' @param lambda.min.ratio If \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, it is the smallest value for \code{lambda}, as a fraction of the upperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, and 0.05 \code{method = "ct"}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @export
huge.ct = function(x, nlambda = NULL, lambda.min.ratio = NULL, lambda = NULL, verbose = TRUE)
{
  gcinfo(FALSE)
  n = nrow(x);
  d = ncol(x);
  fit = list()
  fit$cov.input = isSymmetric(x);
  if(fit$cov.input)
  {
    if(verbose) cat("The input is identified as the covariance matrix.\n")
    S = cov2cor(x);
  }
  if(!fit$cov.input)
  {
    x = scale(x)
    S = cor(x)
  }

  rm(x)
  gc()
  diag(S) = 0
  S = abs(S)
    S.rank = order(S,decreasing = TRUE)
  gc()

  if(is.null(lambda))
  {
    if(is.null(nlambda))
      nlambda = 20
    if(is.null(lambda.min.ratio))
      lambda.min.ratio = 0.05

    density.max = lambda.min.ratio*d*(d-1)/2
    density.min = 1
    density.all = ceiling(seq(density.min,density.max,length = nlambda))*2
    fit$sparsity = density.all/d/(d-1)
    fit$lambda = S[S.rank[density.all]]
    rm(density.max,lambda.min.ratio,density.min,S)
    gc()

    fit$path = list()
    for(i in 1:nlambda)
    {
      fit$path[[i]] = Matrix(0,d,d)
      fit$path[[i]][S.rank[1:density.all[i]]] = 1
      if(verbose)
      {
          cat(paste(c("Conducting the graph estimation via correlation thresholding (ct) ....in progress:", floor(100*i/nlambda), "%"), collapse=""), "\r")
              flush.console()
            }
    }
    rm(density.all,nlambda,S.rank)
    gc()
  }

  if(!is.null(lambda))
  {
    nlambda = length(lambda)
    fit$path = list()
    fit$sparsity = rep(0,nlambda)
    for(i in 1:nlambda)
    {
      fit$path[[i]] = Matrix(0,d,d)
      fit$path[[i]][S > lambda[i]] = 1
      fit$sparsity[i] = sum(fit$path[[i]])/d/(d-1)
      if(verbose)
      {
          mes <- paste(c("Conducting the graph estimation via correlation thresholding (ct)....in progress:", floor(100*i/nlambda), "%"), collapse="")
          cat(mes, "\r")
              flush.console()
            }
    }
    fit$lambda = lambda
    rm(S,lambda)
    gc()
  }

  if(verbose)
  {
        cat("Conducting the graph estimation via correlation thresholding (ct)....done.             \r\n")
        flush.console()
    }
  gc()
  return(fit)
}
